%% This script runs partial least squares regression (PLSR) to predict parcel-wise
% average regional autocorrelation (AC) from neurochemical/structural information.
%
% Overview                                                                    
%   1. Main model   : random 10-fold cross-validation (CV) across parcels
%   2. Significance : spin test (spatial permutations of the maps), same CV
%   3. Robustness   : distance-dependent CV (train on near, test on distant parcels)
%   4. Stability    : parcel bootstrap of VIP scores and PLS coefficients
%   5. Summary table and figure
% Steps 3-4 use pls_uncertainty.m, identical to the power and age-map analyses.
% Required functions: cvR2_pls.m, distcv_pls.m, bootVIP_pls.m,
%                     pls_uncertainty.m, plsBetaBoot.m, spin_test_old.m
% CS 2026

%% Settings
kfold = 10;
ncomp = 10;          % number of components for partial least squares (PLS)
nRep  = 50;          % CV repetitions 
nSpin = 1000;
seed  = 1234;
baselinetype  = 'zscore';  % same as power
trainFraction = 0.75;     % distance CV: nearest 75% train, most distant 25% test   
nBoot         = 10000;    % parcel bootstrap resamples                            
centroidFile  = 'Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';  
permFile      = 'spin_permid_schaefer200.mat';   % same spins for all analyses    

%% Load neurochemical/structual maps
load('allmaps_concat2_selected.mat')

% Don't zscore globally, but do so in CV loop - keep raw maps
X_raw = double(allmaps);   % original predictors (200 x nMaps)
feats = 1:size(X_raw,2);   % use all maps
labels = string(all_labels{:,1});                                                 

%% Get MEG data: autocorrelation, adjusted for covariates
load('adjusted_MEGACall_sorted.mat')

% get average values (parcel x lag);
gACav = squeeze(mean(all_adj, 1));
Y_raw = gACav;                                   % (200 x 40)

%% Quick check on explained variance (full data, exploratory)
allmapsz = zscore(double(allmaps),0,1);
[xl,yl,XS,YS,beta,PCTVAR,mse,stat] = plsregress(allmapsz,Y_raw,30); % 30 components
figure
plot(cumsum(100*PCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in y'); % exploratory

%% 1. Main model: cross-validated PLS (R2, VIP, ranks, betas)
rng(seed);
[average_R2, average_vip, average_rank, det] = ...
    cvR2_pls(feats, X_raw, Y_raw, ncomp, kfold, nRep, baselinetype);

R2_all       = det.R2_rep;
medrank      = det.medrank_rep;
parcel_rank  = det.parcel_rank;
beta_all     = det.beta_all;                 % reps x (1+nMaps) x lags, intercept included
rmse_parcel  = det.rmse_parcel;
rmse_matrix  = det.rmse_matrix;
Y1           = det.Y_pred_last;              % predictions of the last repetition

average_parcelrank = mean(parcel_rank,2);
average_weights    = squeeze(mean(beta_all(:,2:end,:),1));

fprintf('CV R2 = %.3f, median rank = %.1f\n', average_R2, average_rank);

%% Save VIP ranking and prediction results
[VIPs,I] = sort(average_vip, 'descend');
maps_selected = all_labels(I,:);
maps_selected.scores = VIPs;
maps_selected.predvariable = repmat('AC_10ncomp_',length(VIPs),1);
writetable(maps_selected, 'maps_allselected_ACall_zscore_10ncomp_noleak.xlsx')

Y = Y_raw;
save('ACall_allselected_zscore_10ncomp_pls_noleak.mat', 'average_R2', 'average_rank', ...
    'average_vip', 'average_weights', 'beta_all', 'maps_selected', 'parcel_rank', 'rmse_matrix', 'rmse_parcel')

Y1z = det.Y_pred_norm_last;   % predicted, z-scored per lag (fold-specific training stats)
Yz  = det.Y_true_norm_last;   % observed,  z-scored per lag (same stats)
save('orig_pred_ACall_allselected_zscore_10ncomp_noleak.mat', 'Y', 'Y1', 'Yz', 'Y1z')

%% 2. Significance: spin test (same function and settings as observed model)
% Spins are generated once with a fixed seed and saved, so that all analyses   
% (power, AC, age maps, distance CV) use identical spatial permutations.
if isfile(permFile)                                                           
    load(permFile, 'permid');
else
    rng(seed);
    [~, ~, permid] = spin_test_old(rand(1,200), rand(1,200), 'parcellation_name','schaefer_200');
    save(permFile, 'permid');
end
nSpin = min(nSpin, size(permid,2));                                           

rsqsh        = zeros(nSpin,1);
spin_medrank = zeros(nSpin,1);
vip_spin     = zeros(nSpin, numel(feats));                                   
parfor k = 1:nSpin
    rng(seed + k);
    [r2k, vipk, rankk] = cvR2_pls(feats, X_raw(permid(:,k),:), Y_raw, ncomp, kfold, nRep, baselinetype);  
    rsqsh(k) = r2k;  spin_medrank(k) = rankk;  vip_spin(k,:) = vipk';         
end
results_single_spin.rsqsh_spin   = rsqsh;
results_single_spin.spin_medrank = spin_medrank;

spin_99prctile      = prctile(rsqsh, 99);
spin_99prctile_rank = prctile(spin_medrank, 99);
pval_single  = (1 + sum(rsqsh        >= average_R2))   / (nSpin + 1);   % one-sided
pval_medrank = (1 + sum(spin_medrank >= average_rank)) / (nSpin + 1);   % one-sided
pval_vip     = (1 + sum(vip_spin >= average_vip', 1))' / (nSpin + 1);   % one-sided, per map  
qval_vip     = bh_fdr(pval_vip);                                        % FDR across maps    
fprintf('spin p: R2 = %.4f, rank = %.4f; %d maps with VIP spin-FDR < .05\n', ...
        pval_single, pval_medrank, nnz(qval_vip < .05));                                   

% add significance test using spins
save('ACall_allselected_zscore_10ncomp_pls_noleak.mat', 'average_R2', 'average_rank', ...
    'average_vip', 'average_weights', 'beta_all', 'maps_selected', 'parcel_rank', 'rmse_matrix', 'rmse_parcel',...
    'results_single_spin', 'pval_single', 'pval_medrank', 'spin_99prctile_rank', ...
    'pval_vip', 'qval_vip')                                                              

%% 3.-4. Robustness (distance-dependent CV) and stability (parcel bootstrap)   
% Identical to the power and age-map analyses (pls_uncertainty.m):
%  (a) Distance-dependent CV: each of the 200 parcels serves once as seed; the
%      75% of parcels nearest to it train the model, the 25% most distant are
%      predicted. Significance: spin null (maps rotated, splits kept fixed).
%      The 200 splits overlap, so their spread is a robustness range, not a CI.
%  (b) Parcel bootstrap: parcels resampled with replacement and PLS refitted
%      -> 95% CI of VIP scores and ranks, P(VIP > 1), P(top 10).
Tc     = readtable(centroidFile, 'VariableNamingRule', 'preserve');
coords = [Tc.R, Tc.A, Tc.S];             % parcel centroids, same order as maps
assert(size(coords,1) == size(X_raw,1), 'Centroids and maps differ in parcel number.');

opt = struct('ncomp', ncomp, 'baselinetype', baselinetype, 'trainFraction', trainFraction, ...
             'nSpin', nSpin, 'nBoot', nBoot, 'seed', seed);
U = pls_uncertainty(X_raw, Y_raw, coords, permid, opt);

% PLS coefficients (heatmap) with parcel-bootstrap 95% CI, as in the age-map analysis
B = plsBetaBoot(X_raw, Y_raw, ncomp, nBoot, seed, baselinetype);

%% 5. Summary: VIP ranking with significance, robustness and stability 
vipTab = table(labels, average_vip, pval_vip, qval_vip, ...
    U.boot.ci(:,1), U.boot.ci(:,2), U.boot.pVIPgt1, U.boot.pTop10, ...
    U.dist.vip, U.dist.qVIP, U.dist.rankMed, U.dist.pTop10, ...
    'VariableNames', {'Map','VIP','pSpin','qFDR', ...
    'VIP_boot_lo','VIP_boot_hi','P_VIPgt1_boot','P_top10_boot', ...
    'VIP_dist','qFDR_dist','Rank_dist_med','P_top10_dist'});
vipTab = sortrows(vipTab, 'VIP', 'descend');
disp(vipTab(1:15,:));
writetable(vipTab, 'maps_allselected_ACall_zscore_10ncomp_robustness.xlsx')

fprintf('\nSummary (AC)\n');
fprintf('  Random CV  : R2 = %.3f (p_spin = %.4f), median rank %.1f of 200\n', ...
        average_R2, pval_single, average_rank);
fprintf('  Distance CV: R2 = %.3f (p_spin = %.4f; splits 5-95%%: %.2f to %.2f; test-mean baseline: %.3f)\n', ...
        U.dist.R2_pooled, U.dist.pR2, U.dist.R2_p5_p95, U.dist.R2_pooled_testmean);
fprintf('  Bootstrap  : %d maps with VIP > 1 in >= 95%% of resamples\n', nnz(U.boot.pVIPgt1 >= .95));

save('ACall_allselected_zscore_10ncomp_pls_noleak.mat', ...
     'U', 'B', 'vipTab', 'opt', 'permid', '-append')

% VIP with 95% bootstrap CI (supplementary figure)
figure('Color','w','Position',[50 50 1400 400]);
[~, o] = sort(average_vip, 'descend');
plot([1:numel(o); 1:numel(o)], U.boot.ci(o,:)', 'k-'); hold on
plot(1:numel(o), average_vip(o), 'o', 'MarkerFaceColor', 'r'); yline(1, '--');
set(gca, 'XTick', 1:numel(o), 'XTickLabel', labels(o), 'FontSize', 6); xtickangle(60);
ylabel('VIP (95% bootstrap CI)'); title('AC: VIP scores with parcel-bootstrap CI');
exportgraphics(gcf, 'AC_VIP_bootCI.png', 'Resolution', 200);

%% local functions
function q = bh_fdr(p)           % Benjamini-Hochberg FDR
p = p(:); m = numel(p); [ps, idx] = sort(p);
qs = min(1, flipud(cummin(flipud(ps .* m ./ (1:m)'))));
q = nan(m,1); q(idx) = qs;
end
