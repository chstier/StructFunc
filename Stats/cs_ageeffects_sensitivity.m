%% Sensitivity of the MEG power age-effect model to donor characteristics
%  (1) Spearman correlations of VIP with donor age, sample size and
%      |donor age - MEG mean age|; all maps and PET maps; BH-FDR across 6 tests
%  (2) Refits without the evolutionary expansion map and without human maps
%      from cohorts with n < 15 (spin test as in the main analysis)
%  (3) Regression of donor age from the z-scored PET maps (Hansen et al., 2022);
%      variants: all maps with donor age; log10 sample size instead of age
%  All PLS models use the settings, seed (-> identical CV folds) and spin
%  permutations of the main analysis.
%  Requires: Pow_agemap_standalone.mat, sens_load.m, cvR2_pls.m
%  CS 2026
clear; clc;

baseDir = '/mnt/ceph-hdd/projects/scc_uknl_focke/projects/cstier/aging_multiv';
resFile = fullfile(baseDir, 'results', 'agemap_standalone', 'Pow_zY_agemap_standalone.mat');
outDir  = fullfile(baseDir, 'results', 'agemap_sensitivity');
outXls  = fullfile(outDir, 'sens_donor.xlsx');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if exist(outXls, 'file'), delete(outXls); end

megAge = 53;    % mean age of the MEG sample in years (not the centred age column!)
nMin   = 15;    % small-sample threshold
topK   = 10;    % top-K overlap of VIP rankings
assert(megAge > 18, 'megAge looks centred.');
if isempty(gcp('nocreate')), parpool; end

%% Data and main-analysis results
S = sens_load(resFile);
X = S.X;  Y = S.ageMap;  permid = S.permid;  M = S.meta;  labels = string(S.labels(:));
R = load(resFile, 'R2obs', 'R2null', 'vipObs', 'opt');
opt   = R.opt;  vipFull = R.vipObs(:);
nMaps = size(X, 2);  nSpin = numel(R.R2null);
pFull = (1 + sum(R.R2null >= R.R2obs)) / (nSpin + 1);

isPET = M.type == "PET";  isComp = M.type == "composite";
age   = double(M.age_mean);  nDon = double(M.n);
rankFull = tiedrank(-vipFull);                       % 1 = highest VIP

% reruns must reproduce the main model (same seed -> same folds)
rng(opt.seed);
[R2chk, vipChk] = cvR2_pls(1:nMaps, X, Y, opt.ncomp, opt.kfold, opt.nRep, 'zscore');
if abs(R2chk - R.R2obs) > 1e-6 || max(abs(vipChk(:) - vipFull)) > 1e-6
    warning('Main model not reproduced (R2 %.4f vs %.4f): check seed/settings.', R2chk, R.R2obs);
end

%% (1) Donor characteristics vs VIP 
pred = table(age, nDon, abs(age - megAge), 'VariableNames', {'donorAge', 'sampleSize', 'ageGap'});
grp  = {'all (non-composite)', ~isComp;  'PET only', isPET};
rows = {};
for g = 1:size(grp, 1)
    for i = 1:width(pred)
        x  = pred{:, i};  ok = grp{g, 2} & ~isnan(x);
        [rho, p] = corr(x(ok), vipFull(ok), 'type', 'Spearman');
        rows(end+1, :) = {string(grp{g,1}), string(pred.Properties.VariableNames{i}), ...
                          nnz(ok), rho, p}; %#ok<SAGROW>
    end
end
Tdon = cell2table(rows, 'VariableNames', {'group', 'predictor', 'nMaps', 'rho', 'p'});
Tdon.pFDR = bhFDR(Tdon.p);
disp(Tdon);

Supp = [M(:, {'label', 'label_data', 'type', 'n', 'n_male', 'age_mean', 'age_sd', ...
              'age_min', 'age_max', 'human', 'adult'}), ...
        table(vipFull, rankFull, 'VariableNames', {'VIP', 'rank'})];
Supp = sortrows(Supp, 'rank');

%% (2) Refits 
smallN = M.human & nDon < nMin & ~isComp;
sets = {'Without evolutionary expansion',               M.human
        sprintf('Without human cohorts n < %d', nMin),  ~smallN};
fprintf('\nMaps from cohorts n < %d (%d): %s\n', nMin, nnz(smallN), strjoin(labels(smallN), ', '));

Tref = table("Full model", nMaps, R.R2obs, prctile(R.R2null, 95), pFull, NaN, ...
             topList(vipFull, 1:nMaps, labels), "", 'VariableNames', ...
             {'Model', 'nMaps', 'R2', 'R2null95', 'pSpin', 'rhoVIP', 'Top5', 'Excluded'});
for s = 1:size(sets, 1)
    feats = find(sets{s, 2});
    rng(opt.seed);
    [R2, vip] = cvR2_pls(feats, X, Y, opt.ncomp, opt.kfold, opt.nRep, 'zscore');
    vip    = pickVIP(vip, feats, nMaps);
    R2null = spinNull(feats, X, Y, permid, nSpin, opt);
    Tref(end+1, :) = {string(sets{s,1}), numel(feats), R2, prctile(R2null, 95), ...
        (1 + sum(R2null >= R2)) / (nSpin + 1), corr(vip, vipFull(feats), 'type', 'Spearman'), ...
        topList(vip, feats, labels), strjoin(labels(~sets{s,2}), ', ')}; %#ok<SAGROW>
end
disp(Tref(:, 1:6));

%% (3) Regression of donor age (similar to Hansen et al., 2022)
Xz = zscore(X);
V = {'donor age, PET maps (primary)', age,         isPET & ~isnan(age), true
     'donor age, all maps',           age,         ~isnan(age),         false
     'log10 sample size, PET maps',   log10(nDon), isPET & nDon > 0,    false};
reg = struct([]);
for v = 1:size(V, 1)
    sel = V{v, 3};  xc = V{v, 2}(sel);  xc = xc - mean(xc);
    A   = Xz(:, sel);
    slope = (A * xc) / (xc' * xc);                   % OLS slope per parcel
    Xr  = Xz;  Xr(:, sel) = A - slope * xc';         % parcel mean across maps preserved

    % similarity of original and regressed maps
    rMap = diag(corr(A, Xr(:, sel)));
    C0 = corr(A);  C1 = corr(Xr(:, sel));  ut = triu(true(size(C0)), 1);
    selLab = labels(sel);  [rMin, iMin] = min(rMap);

    % PLS with identical settings and folds
    rng(opt.seed);
    [R2, vip] = cvR2_pls(1:nMaps, Xr, Y, opt.ncomp, opt.kfold, opt.nRep, 'zscore');
    vip = vip(:);  pSpin = NaN;
    if V{v, 4}
        R2null = spinNull(1:nMaps, Xr, Y, permid, nSpin, opt);
        pSpin  = (1 + sum(R2null >= R2)) / (nSpin + 1);
    end
    [~, o0] = sort(vipFull, 'descend');  [~, o1] = sort(vip, 'descend');

    reg(v).Variant      = string(V{v, 1});
    reg(v).nRegressed   = nnz(sel);
    reg(v).r_pooled     = corr(A(:), reshape(Xr(:, sel), [], 1));
    reg(v).r_map_min    = rMin;
    reg(v).map_min      = selLab(iMin);
    reg(v).r_map_median = median(rMap);
    reg(v).r_similarity = corr(C0(ut), C1(ut));
    reg(v).R2_orig      = R.R2obs;
    reg(v).R2_reg       = R2;
    reg(v).pSpin_orig   = pFull;
    reg(v).pSpin_reg    = pSpin;
    reg(v).rhoVIP_all   = corr(vipFull, vip, 'type', 'Spearman');
    reg(v).rhoVIP_PET   = corr(vipFull(isPET), vip(isPET), 'type', 'Spearman');
    reg(v).topK_overlap = numel(intersect(o0(1:topK), o1(1:topK)));
    reg(v).nVIP1_orig   = nnz(vipFull > 1);
    reg(v).nVIP1_reg    = nnz(vip > 1);
    reg(v).nVIP1_both   = nnz(vipFull > 1 & vip > 1);

    if v == 1                                        % detailed table, primary variant
        rFull = nan(nMaps, 1);  rFull(sel) = rMap;
        VIPtab = table(labels, M.type, age, sel, vipFull, rankFull, vip, tiedrank(-vip), ...
                       tiedrank(-vip) - rankFull, rFull, 'VariableNames', {'Map', 'Type', ...
                       'DonorAge', 'Regressed', 'VIP_orig', 'Rank_orig', 'VIP_reg', 'Rank_reg', ...
                       'RankChange', 'r_orig_vs_reg'});
        VIPtab = sortrows(VIPtab, 'Rank_orig');
    end
end
Treg = struct2table(reg);
disp(Treg);

%% save
writetable(Tdon,   outXls, 'Sheet', 'donor_stats');
writetable(Supp,   outXls, 'Sheet', 'donor_table');
writetable(Tref,   outXls, 'Sheet', 'refits');
writetable(Treg,   outXls, 'Sheet', 'regression');
writetable(VIPtab, outXls, 'Sheet', 'regression_VIP');
save(fullfile(outDir, 'sens_donor.mat'), 'Tdon', 'Supp', 'Tref', 'Treg', 'VIPtab', ...
     'megAge', 'nMin', 'opt', 'nSpin');

%% numbers for the manuscript
fprintf('\n(1) max |rho| = %.2f, min p = %.2f, min q = %.2f (%d tests)\n', ...
        max(abs(Tdon.rho)), min(Tdon.p), min(Tdon.pFDR), height(Tdon));
for s = 2:height(Tref)
    fprintf('(2) %s: %d maps, R2 = %.3f (null95 %.2f, p_spin = %.3f), rho = %.2f\n', ...
            Tref.Model(s), Tref.nMaps(s), Tref.R2(s), Tref.R2null95(s), Tref.pSpin(s), Tref.rhoVIP(s));
end
P = reg(1);
fprintf(['(3) pooled r = %.2f; per-map r >= %.2f (%s); similarity r = %.2f\n' ...
         '    R2 %.3f (p = %.3f) -> %.3f (p = %.3f); rho all %.2f, PET %.2f; top-%d %d/%d; VIP>1 %d -> %d (%d shared)\n'], ...
        P.r_pooled, P.r_map_min, P.map_min, P.r_similarity, P.R2_orig, P.pSpin_orig, P.R2_reg, ...
        P.pSpin_reg, P.rhoVIP_all, P.rhoVIP_PET, topK, P.topK_overlap, topK, P.nVIP1_orig, ...
        P.nVIP1_reg, P.nVIP1_both);
disp(VIPtab(abs(VIPtab.RankChange) >= 5, :));        % largest rank changes

%% =================== local functions ===================
function R2null = spinNull(feats, X, Y, permid, nSpin, opt)
% identical to the spin loop of the main analysis (same seed per spin)
R2null = zeros(nSpin, 1);
parfor k = 1:nSpin
    rng(opt.seed + k);
    R2null(k) = cvR2_pls(feats, X(permid(:, k), :), Y, opt.ncomp, opt.kfold, opt.nRep, 'zscore');
end
end

function v = pickVIP(vip, feats, nMaps)
% cvR2_pls may return VIP for all maps or only for the selected maps
vip = vip(:);
if numel(vip) == nMaps && numel(feats) < nMaps, v = vip(feats); else, v = vip; end
end

function s = topList(vip, feats, labels)
[~, o] = sort(vip(:), 'descend');
s = strjoin(labels(feats(o(1:min(5, numel(o))))), ', ');
end

function q = bhFDR(p)
% Benjamini-Hochberg adjusted p-values
p = p(:);  m = numel(p);
[ps, o] = sort(p);
qs = min(flipud(cummin(flipud(ps .* m ./ (1:m)'))), 1);
q = zeros(m, 1);  q(o) = qs;
end