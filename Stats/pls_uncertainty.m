function U = pls_uncertainty(X, Y, coords, permid, opt)
% PLS_UNCERTAINTY  Robustness/uncertainty of a PLS model X (maps) -> Y (MEG).
% Applied identically to the main power, main AC and age-map analyses.
%  (a) Distance-dependent CV (distcv_pls.m) with spin null
%      (rows of X spun, spatial splits fixed) and the 5-95% range across
%      the overlapping seed splits (robustness range, not an independent CI)
%  (b) Parcel bootstrap (bootVIP_pls.m): 95% percentile CI of VIP and VIP rank
% opt: ncomp, baselinetype, trainFraction, nSpin, nBoot, seed

nMaps = size(X,2);  feats = 1:nMaps;  nSpin = opt.nSpin;

%% (a) distance-dependent CV + spin null
[R2, vip, mrank, D] = distcv_pls(feats, X, Y, opt.ncomp, coords, opt.trainFraction, opt.baselinetype);
D.R2_pooled = R2;  D.vip = vip;  D.mean_medrank = mrank;

R2null = zeros(nSpin,1);  rankNull = zeros(nSpin,1);  vipNull = zeros(nMaps, nSpin);
parfor k = 1:nSpin
    [r2k, vipk, rankk] = distcv_pls(feats, X(permid(:,k),:), Y, opt.ncomp, coords, ...
                                    opt.trainFraction, opt.baselinetype);
    R2null(k) = r2k;  vipNull(:,k) = vipk;  rankNull(k) = rankk;
end
D.R2null    = R2null;
D.rankNull  = rankNull;
D.pR2       = (1 + sum(R2null   >= D.R2_pooled))    / (nSpin + 1);   % one-sided
D.pRank     = (1 + sum(rankNull >= D.mean_medrank)) / (nSpin + 1);
D.pVIP      = (1 + sum(vipNull  >= D.vip, 2))       / (nSpin + 1);
D.qVIP      = bh_fdr(D.pVIP);
D.R2_p5_p95 = prctile(D.R2_split, [5 95]);
D.weights   = squeeze(mean(D.beta_split(:,2:end,:), 1));  % as average_weights in main scripts

rk = vip_ranks(D.vip_split);
D.rankMed = median(rk, 2);  D.rankP5 = prctile(rk, 5, 2);  D.rankP95 = prctile(rk, 95, 2);
D.pTop10  = mean(rk <= 10, 2);
D.pVIPgt1 = mean(D.vip_split > 1, 2);

%% (b) parcel bootstrap of VIP
B.vip     = bootVIP_pls(feats, X, Y, opt.ncomp, opt.baselinetype, opt.nBoot, opt.seed);
B.ci      = prctile(B.vip, [2.5 97.5], 2);
rk        = vip_ranks(B.vip);
B.rankCI  = prctile(rk, [2.5 97.5], 2);
B.pVIPgt1 = mean(B.vip > 1, 2);
B.pTop10  = mean(rk <= 10, 2);

U.dist = D;  U.boot = B;
fprintf('Distance CV: R2 = %.3f (splits 5-95%%: %.2f to %.2f), p_spin = %.4f; median rank %.1f of %d\n', ...
        D.R2_pooled, D.R2_p5_p95, D.pR2, D.mean_medrank, D.nTest);
fprintf('Bootstrap: %d maps with VIP > 1 in >= 95%% of resamples\n', nnz(B.pVIPgt1 >= .95));
end

function R = vip_ranks(vip)      % maps x samples -> ranks (1 = highest VIP)
R = nan(size(vip));
for s = 1:size(vip,2)
    [~, o] = sort(vip(:,s), 'descend');  R(o,s) = 1:size(vip,1);
end
end

function q = bh_fdr(p)           % Benjamini-Hochberg
p = p(:); m = numel(p); [ps, idx] = sort(p);
qs = min(1, flipud(cummin(flipud(ps .* m ./ (1:m)'))));
q = nan(m,1); q(idx) = qs;
end