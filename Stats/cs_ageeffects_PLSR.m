%% Which maps spatially align with age-related changes in MEG power / AC?
% Step 1: age effects per parcel x feature: y ~ age + sex + tiv
% Step 2: spatial correlation of each map with the age map (spin test, FDR + FWE)
% Step 3: PLS maps -> age map, identical to the main analysis (cvR2_pls.m)
% Step 4: uncertainty of Step 3, identical to the main analyses (pls_uncertainty.m) 
%         a) distance-dependent CV + spin null, b) parcel bootstrap of VIP
% Normalization of the age map (partial r): mean-centred across (training) parcels 
% within each frequency bin / lag, no variance scaling, for BOTH power and AC.
% Outputs: <prefix>_agemap_standalone.mat, <prefix>_agemap_for_R.mat,
%          <prefix>_agemap_PLS_VIP.xlsx, agemap_maplabels.csv
% CS 2026
clear; clc;

metric = 'AC';                                       % 'power' or 'AC'
normY  = 'zscore';                                 % 'absolute' or 'zscore'

baseDir      = '/mnt/ceph-hdd/projects/scc_uknl_focke/projects/cstier/aging_multiv';
outDir       = fullfile(baseDir, 'results', 'agemap_standalone');
mapFile      = 'allmaps_concat2_selected.mat';
designFile   = fullfile(baseDir, 'design_all_age_sex_icv.csv');
centroidFile = 'Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';  
permFile     = 'spin_permid_schaefer200.mat';           % same spins as all analyses
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Settings: identical to the main analysis
opt.ncomp = 10;  opt.kfold = 10;  opt.nRep = 50;  opt.seed = 1234;
opt.kPath = 1:10;                % components for the R2 curve
opt.trainFraction = 0.75;        % distance-dependent CV
opt.nBoot = 10000;               % parcel bootstrap (as main analysis)
opt.baselinetype = normY;        % same for power and AC (both targets are partial r)
nSpin     = 1000;
alphaFDR  = 0.05;
doQuadPLS = false;               % true: also run Steps 2-4 on the quadratic age map

if isempty(gcp('nocreate')), parpool; end

%% Data
Mf     = load(mapFile, 'allmaps', 'all_labels');
X      = double(Mf.allmaps);                            % parcels x maps
labels = string(Mf.all_labels{:,1});
[nParc, nMaps] = size(X);

switch metric                                          
    case 'power'
        S = load('MEGpowspc_sorted_Schaefer200.mat', 'power_sorted');   % UNADJUSTED power
        P = double(S.power_sorted);                     % subjects x parcels x freqs
        L = log(P ./ mean(P, 3));                       % log relative power (as main)
        xvals = 1:60;  xlab = 'Frequency (Hz)';
        prefix = 'Pow';                                 % baselinetype set above
    case 'AC'
        S = load('MEGdata_sorted_Schaefer200_2.mat', 'dat_sorted');       % UNADJUSTED AC
        ACdat=S.dat_sorted(:,:,84:123);
        L = double(ACdat);                        % subjects x parcels x lags (1-40)
        xvals = (1:40) * 1000/300;  xlab = 'Lag (ms)';
        prefix = 'AC';                                  % baselinetype set above
    otherwise
        error('metric must be ''power'' or ''AC''');
end
if strcmpi(normY, 'zscore'), prefix = [prefix '_zY']; end 

[~, nParcP, nF] = size(L);
assert(nParcP == nParc, 'Parcel number differs between maps and MEG data.');
assert(nF == numel(xvals), 'Number of MEG features (%d) differs from xvals (%d).', nF, numel(xvals));
Ymeg = reshape(L, size(L,1), nParc*nF);                 % subjects x features

D = readmatrix(designFile);                             % columns: age, sex, tiv
ok = all(~isnan(D(:,1:3)), 2) & all(isfinite(Ymeg), 2);
if any(~ok), fprintf('Excluding %d subjects with missing data.\n', nnz(~ok)); end
Ymeg = Ymeg(ok,:);  D = D(ok,:);  nSubj = nnz(ok);

ageMeanMEG = mean(D(:,1));
fprintf('Age column: mean = %.2f, SD = %.2f, range %.2f-%.2f (N = %d)\n', ...
        ageMeanMEG, std(D(:,1)), min(D(:,1)), max(D(:,1)), nSubj);
if abs(ageMeanMEG) < 1 || min(D(:,1)) < 0
    warning(['Age column looks centred/standardised - ageMeanMEG is not the sample mean age. ' ...
             'Results are unaffected, but set the age gap manually.']);
end

ageC = D(:,1) - mean(D(:,1));
T = table(ageC, ageC.^2, categorical(D(:,2)), D(:,3), ...
          'VariableNames', {'age','ageSq','sex','tiv'});

%% Step 1: linear + quadratic age effects, dBIC (fitlm per feature) 
[rLin, pLin, bicLin] = ageEffect(Ymeg, T, 'y ~ age + sex + tiv',         'age');
[rQua, pQua, bicQua] = ageEffect(Ymeg, T, 'y ~ age + ageSq + sex + tiv', 'ageSq');
dBIC = bicQua - bicLin;                                 % < 0 favours the quadratic model

ageMap  = reshape(rLin, nParc, nF);                     % parcels x freqs/lags
qAge    = reshape(bhFDR(pLin), nParc, nF);
quadMap = reshape(rQua, nParc, nF);
qQuad   = reshape(bhFDR(pQua), nParc, nF);
dBIC    = reshape(dBIC, nParc, nF);
ageMeanByFeat = mean(ageMap, 1);                        % global (cortex-wide) age effect per bin/lag,
                                                        % removed by centring
fprintf('Linear age effect FDR<.05: %.1f%% of features\n', 100*mean(qAge(:) < alphaFDR));
fprintf('Quadratic age effect FDR<.05: %.1f%%; dBIC < -2: %.1f%%\n', ...
        100*mean(qQuad(:) < alphaFDR), 100*mean(dBIC(:) < -2));

% sanity check against fitlm on a few random features
T  = table(ageC, categorical(D(:,2)), D(:,3), 'VariableNames', {'age','sex','tiv'});
rs = RandStream('mt19937ar', 'Seed', opt.seed);
for k = randperm(rs, nParc*nF, 5)
    T.y = Ymeg(:,k);  m = fitlm(T, 'y ~ age + sex + tiv');
    tt  = m.Coefficients{'age','tStat'};
    assert(abs(tt/sqrt(tt^2 + m.DFE) - rLin(k)) < 1e-8, 'OLS and fitlm disagree.');
end

%% Spins (fixed seed, saved and shared across analyses) 
if isfile(permFile)
    load(permFile, 'permid');
else
    rng(opt.seed);
    [~, ~, permid] = spin_test_old(rand(1,nParc), rand(1,nParc), 'parcellation_name', 'schaefer_200');
    save(permFile, 'permid');
end
nSpin = min(nSpin, size(permid, 2));  opt.nSpin = nSpin;

%% Steps 2-3: alignment of maps with the age map
fprintf('\n=== Linear age map (%s, Y: %s) ===\n', metric, opt.baselinetype);   
AL = runAlignment(X, ageMap, permid, nSpin, opt, alphaFDR);
AQ = [];
if doQuadPLS
    fprintf('\n=== Quadratic age map ===\n');
    AQ = runAlignment(X, quadMap, permid, nSpin, opt, alphaFDR);
end

%% Step 4: uncertainty (distance-dependent CV, parcel bootstrap) 
Tc     = readtable(centroidFile, 'VariableNamingRule', 'preserve');
coords = [Tc.R, Tc.A, Tc.S];
fprintf('\n=== Uncertainty, linear age map ===\n');
UL = pls_uncertainty(X, ageMap, coords, permid, opt);   % uses opt.baselinetype (check, see notes)
UQ = [];
if doQuadPLS, UQ = pls_uncertainty(X, quadMap, coords, permid, opt); end

%% Step 5: PLS regression coefficients + parcel-bootstrap CI 
fprintf('\n=== PLS coefficients (bootstrap), linear age map ===\n');
BL = plsBetaBoot(X, ageMap, opt.ncomp, opt.nBoot, opt.seed, opt.baselinetype);
BQ = [];
if doQuadPLS, BQ = plsBetaBoot(X, quadMap, opt.ncomp, opt.nBoot, opt.seed, opt.baselinetype); end

%% VIP table
vipTab = table(labels, AL.vip, AL.pVIP, AL.qVIP, ...
    UL.boot.ci(:,1), UL.boot.ci(:,2), UL.boot.rankCI(:,1), UL.boot.rankCI(:,2), UL.boot.pVIPgt1, ...
    UL.dist.vip, UL.dist.rankMed, UL.dist.rankP5, UL.dist.rankP95, UL.dist.pTop10, ...
    UL.dist.pVIP, UL.dist.qVIP, ...
    'VariableNames', {'Map','VIP','pSpin','qFDR', ...
    'VIP_boot_lo','VIP_boot_hi','Rank_boot_lo','Rank_boot_hi','P_VIPgt1_boot', ...
    'VIP_dist','Rank_dist_med','Rank_dist_p5','Rank_dist_p95','P_top10_dist', ...
    'pSpin_dist','qFDR_dist'});
vipTab = sortrows(vipTab, 'VIP', 'descend');
disp(vipTab(1:15,:));

%% quick-look figures
figure('Color', 'w', 'Position', [50 50 1400 800]);
subplot(2,2,1); imagesc(xvals, 1:nParc, ageMap); axis xy; colorbar;
clim([-1 1]*max(abs(ageMap(:)))); xlabel(xlab); ylabel('Parcel'); title('Age effect (partial r)');
subplot(2,2,2); imagesc(xvals, 1:nMaps, AL.rSp); axis xy; colorbar;
clim([-1 1]*max(abs(AL.rSp(:)))); set(gca, 'YTick', 1:nMaps, 'YTickLabel', labels, 'FontSize', 6);
xlabel(xlab); title('Map vs age map');
subplot(2,2,3); histogram(AL.R2null); xline(AL.R2cv, 'r', 'LineWidth', 2);
xlabel('CV R^2'); title(sprintf('Spin null, p = %.4f', AL.pR2));
subplot(2,2,4); plot(opt.kPath, AL.R2path, 'o-'); xline(opt.ncomp, '--');
xlabel('PLS components'); ylabel('CV R^2'); title('R^2 vs. components');
exportgraphics(gcf, fullfile(outDir, [prefix '_agemap_quicklook.png']), 'Resolution', 200);

figure('Color', 'w', 'Position', [50 50 1400 400]);                  % NEW
[~, o] = sort(AL.vip, 'descend');
plot([1:nMaps; 1:nMaps], UL.boot.ci(o,:)', 'k-'); hold on
plot(1:nMaps, AL.vip(o), 'o', 'MarkerFaceColor', 'r'); yline(1, '--');
set(gca, 'XTick', 1:nMaps, 'XTickLabel', labels(o), 'FontSize', 6); xtickangle(60);
ylabel('VIP (95% bootstrap CI)'); title(sprintf('%s age map: VIP with parcel-bootstrap CI', metric));
exportgraphics(gcf, fullfile(outDir, [prefix '_agemap_VIP_CI.png']), 'Resolution', 200);

%% save for MATLAB / sensitivity scripts
writetable(vipTab, fullfile(outDir, [prefix '_agemap_PLS_VIP.xlsx']));
rSp = AL.rSp;  pSp = AL.pSpin;  qSp = AL.qSpin;  pFWE = AL.pFWE;
R2obs = AL.R2cv;  R2null = AL.R2null;  pR2 = AL.pR2;
vipObs = AL.vip;  pVIP = AL.pVIP;  qVIP = AL.qVIP;
freqs = xvals;                                          % legacy name (lags in ms for AC)
save(fullfile(outDir, [prefix '_agemap_standalone.mat']), ...
     'metric', 'ageMap', 'qAge', 'quadMap', 'qQuad', 'dBIC', 'rSp', 'pSp', 'qSp', 'pFWE', ...
     'R2obs', 'R2null', 'pR2', 'vipObs', 'pVIP', 'qVIP', 'labels', 'freqs', 'xvals', ...
     'permid', 'nSubj', 'ageMeanMEG', 'opt', 'AL', 'AQ', 'UL', 'UQ', 'vipTab', 'BL', 'BQ', ...
     'ageMeanByFeat');                                 

%% export for R (readMat turns '_' into '.') 
ex = struct('freqs', xvals, 'nSubj', nSubj, 'ageMeanMEG', ageMeanMEG, ...
            'lin_r', ageMap, 'lin_q', qAge, 'quad_r', quadMap, 'quad_q', qQuad, 'dBIC', dBIC, ...
            'lin_mean', ageMeanByFeat);                 % CHANGED: mean partial r per bin/lag (panel A line)
fA = {'rSp','pSpin','qSpin','pFWE','vip','pVIP','qVIP','R2cv','pR2','R2null', ...
      'R2path','kPath','ncomp','strucY','strucX'};
fD = {'R2_pooled','R2_pooled_testmean','R2_split','R2_p5_p95','mean_medrank','pR2','pRank', ...
      'R2null','vip','pVIP','qVIP','rankMed','rankP5','rankP95','pTop10','pVIPgt1'};
fB = {'ci','rankCI','pVIPgt1','pTop10'};
ex = addPrefixed(ex, 'L',  AL,      fA);
ex = addPrefixed(ex, 'Ld', UL.dist, fD);               
ex = addPrefixed(ex, 'Lb', UL.boot, fB);                
ex = addPrefixed(ex, 'Lw', BL, {'beta','lo','hi','pSign'});            % NEW: betas
if ~isempty(BQ), ex = addPrefixed(ex, 'Qw', BQ, {'beta','lo','hi','pSign'}); end
if ~isempty(AQ), ex = addPrefixed(ex, 'Q',  AQ,      fA); end
if ~isempty(UQ), ex = addPrefixed(ex, 'Qd', UQ.dist, fD); ex = addPrefixed(ex, 'Qb', UQ.boot, fB); end
save(fullfile(outDir, [prefix '_agemap_for_R.mat']), '-struct', 'ex', '-v6');
writetable(table(labels, 'VariableNames', {'annotation'}), fullfile(outDir, 'agemap_maplabels.csv'));
fprintf('\nSaved to %s\n', outDir);

%% =================== local functions ===================
function A = runAlignment(X, Ymap, permid, nSpin, opt, alphaFDR)
% Spatial correlations (spin FDR + FWE) and PLS (CV R2, spin null, VIP, R2 path)
nMaps = size(X, 2);  feats = 1:nMaps;

% --- Step 2: map x feature correlations; spin p (two-sided) and max-statistic FWE ---
A.rSp   = corr(X, Ymap);
nExceed = zeros(size(A.rSp));  maxNull = zeros(nSpin, 1);
for k = 1:nSpin
    rk = abs(corr(X(permid(:,k),:), Ymap));
    nExceed    = nExceed + (rk >= abs(A.rSp));
    maxNull(k) = max(rk, [], 'all');
end
A.pSpin = (1 + nExceed) / (nSpin + 1);
A.qSpin = reshape(bhFDR(A.pSpin), size(A.pSpin));
A.pFWE  = reshape((1 + sum(maxNull >= abs(A.rSp(:))', 1)) / (nSpin + 1), size(A.rSp));
fprintf('Step 2: %d of %d correlations spin-FDR<.05, %d FWE<.05\n', ...
        nnz(A.qSpin < alphaFDR), numel(A.qSpin), nnz(A.pFWE < .05));

% --- Step 3: PLS, identical settings to the main analysis ---
rng(opt.seed);
[A.R2cv, vip] = cvR2_pls(feats, X, Ymap, opt.ncomp, opt.kfold, opt.nRep, opt.baselinetype);
A.vip = vip(:);

R2null = zeros(nSpin, 1);  vipNull = zeros(nSpin, nMaps);
parfor k = 1:nSpin
    rng(opt.seed + k);
    [r2k, vipk] = cvR2_pls(feats, X(permid(:,k),:), Ymap, opt.ncomp, opt.kfold, opt.nRep, opt.baselinetype);
    R2null(k) = r2k;  vipNull(k,:) = vipk(:)';
end
A.R2null = R2null;
A.pR2    = (1 + sum(R2null >= A.R2cv)) / (nSpin + 1);
A.pVIP   = (1 + sum(vipNull >= A.vip', 1))' / (nSpin + 1);
A.qVIP   = bhFDR(A.pVIP);

kPath = opt.kPath;  R2path = nan(numel(kPath), 1);
parfor i = 1:numel(kPath)
    rng(opt.seed);
    R2path(i) = cvR2_pls(feats, X, Ymap, kPath(i), opt.kfold, opt.nRep, opt.baselinetype);
end
A.R2path = R2path;  A.kPath = kPath(:);  A.ncomp = opt.ncomp;
if any(kPath == opt.ncomp) && abs(R2path(kPath == opt.ncomp) - A.R2cv) > 1e-10
    warning('R2 path at ncomp = %d differs from main R2 (check RNG use in cvR2_pls).', opt.ncomp);
end

% structure coefficients from the same Y normalization as the CV model
Xz = zscore(X);
switch lower(opt.baselinetype)
    case 'zscore',   Yn = zscore(Ymap);
    case 'absolute', Yn = Ymap - mean(Ymap, 1);
end
[~, ~, XS] = plsregress(Xz, Yn, opt.ncomp);
A.strucY = corr(Yn, XS);
A.strucX = corr(Xz, XS);
fprintf('Step 3: CV R2 = %.3f (%d comps), p_spin = %.4f, VIP > 1: %d maps\n', ...
        A.R2cv, opt.ncomp, A.pR2, nnz(A.vip > 1));
end

function ex = addPrefixed(ex, pre, A, f)
for i = 1:numel(f), ex.([pre '_' f{i}]) = A.(f{i}); end
end

function [r, p, bic] = ageEffect(Y, T, formula, term)
nFeat = size(Y, 2);
t = zeros(1, nFeat);  p = t;  dfe = t;  bic = t;
parfor k = 1:nFeat
    Tk = T;  Tk.y = Y(:,k);
    mdl    = fitlm(Tk, formula);
    t(k)   = mdl.Coefficients{term, 'tStat'};
    p(k)   = mdl.Coefficients{term, 'pValue'};
    dfe(k) = mdl.DFE;
    bic(k) = mdl.ModelCriterion.BIC;
end
r = t ./ sqrt(t.^2 + dfe);
end

function q = bhFDR(p)
p = p(:); q = nan(size(p)); ok = ~isnan(p); m = nnz(ok);
if m == 0, return; end
[ps, o] = sort(p(ok));
qs = min(flipud(cummin(flipud(ps .* m ./ (1:m)'))), 1);
tmp = nan(m,1); tmp(o) = qs; q(ok) = tmp;
end
