%% This script predicts individual deviations from the group average (MEG power)
%% Set paths
outputDir = '/home/path/results/individuals';

%% Load respective Neuromaps
load('allmaps_concat2_selected.mat')   % gives 'allmaps' (200 x 55)
allmaps = double(allmaps);
X_full = allmaps; % [200 x 55] parcels x maps

%% Get MEG data for the average (group-level)
load('adjusted_MEGpower_orig.mat') % gives 'all_adj' presumably [nSub x 200 x 60] or similar
% here adjusted for age, age2, sex, tiv
power_sorted = all_adj;
gpowav = squeeze(mean(power_sorted, 1));   % group-average across subjects [200 x 60]
Y_group_log = log(gpowav ./ repmat(mean(gpowav,2), 1, 60));  % [200 x 60]

%% Load individual power spectra (not covariate adjusted)
indiv = load('MEGpowspc_sorted_Schaefer200.mat'); % contains indiv.power_sorted [nSub x 200 x 60]
nSub = size(indiv.power_sorted,1);
nParc = size(indiv.power_sorted,2);
nFreq = size(indiv.power_sorted,3);

%% Build individual log-ratio spectra
all_log = zeros(nSub, nParc, nFreq);
for s = 1:nSub
    pow = squeeze(indiv.power_sorted(s,:,:));   % [200 x 60]
    all_log(s,:,:) = log( pow ./ repmat(mean(pow,2), 1, nFreq) ); % [200 x 60]
end

%% Compute DEV_f = (individual_log - group_mean_log) / group_sd_log
% This is group-based z of log-ratios - here using the corrected mean!!!!
group_mean = squeeze(mean(Y_group_log , 1));  % [200 x 60]
group_sd   = squeeze(std(Y_group_log , 0, 1));% [200 x 60]

% Avoid zeros in group_sd (if a parcel/freq has zero variance)
group_sd(group_sd == 0) = eps;

DEV_f = zeros(nSub, nParc, nFreq);
for s = 1:nSub
    DEV_f(s,:,:) = ( squeeze(all_log(s,:,:)) - group_mean ) ./ group_sd;
end
% Now DEV_f(s,parcel,freq) is the z-scored deviation from the group (log scale).

%% PLS prediction per subject (parcel-wise cross-validation)
% We standardize X (neuromaps) *within each CV training fold* to avoid leakage.
for s = 1:350   % change to 1:nSub for full run
    fprintf('Processing subject %d/%d...\n', s, nSub);

    % Parameters
    kfold = 10;
    ncomp = 10;
    X = X_full;                           % [200 x 55] parcels x maps
    Y_raw = squeeze(DEV_f(s,:,:));        % [200 x 60] target (do not zscore across parcels here)

    numPredictors = size(X,2);            % 55
    numResponses  = size(Y_raw,2);        % 60

    % Preallocate
    beta_all      = zeros(50, numPredictors+1, numResponses);
    all_vip       = zeros(numPredictors, 50);
    R2_all        = zeros(1, 50);
    rmse_parcel   = zeros(nParc, 50);
    rmse_matrix   = zeros(50, nParc, numResponses);
    parcel_rank   = zeros(nParc, 50);
    medrank       = zeros(1, 50);

    for rep = 1:50
        cv = cvpartition(nParc, "Kfold", kfold);
        Y1 = zeros(size(Y_raw));    % collect predictions for all parcels [200 x 60]
        b = zeros(kfold, numPredictors+1, numResponses);
        vipScore = zeros(numPredictors, kfold);

        for ifold = 1:kfold
            trIdx = cv.training(ifold);
            teIdx = cv.test(ifold);

            % Standardize X using training parcels ONLY (prevents leakage)
            Xtr = X(trIdx, :);      % training predictors [nTrain x 55]
            muX  = mean(Xtr, 1);
            sdX  = std(Xtr, 0, 1);
            sdX(sdX == 0) = 1;      % avoid div-by-zero for constant predictors

            Xtr_z = (Xtr - muX) ./ sdX;
            Xte_z = (X(teIdx, :) - muX) ./ sdX;

            % Corresponding Y training/test 
            Ytr = Y_raw(trIdx, :);
            % Fit PLS on Xtr_z, Ytr
            [XL, yl, XS, YS, beta, PCTVAR, mse, stats] = plsregress(Xtr_z, Ytr, ncomp);

            % store fold betas (beta size = (numPredictors+1) x numResponses)
            b(ifold,:,:) = beta;

            % Predict on test parcels (apply intercept + scaled Xte_z)
            Y1(teIdx, :) = [ones(sum(teIdx),1) Xte_z] * beta;

            % VIP: use stats from plsregress (works on scaled Xtr_z)
            W0 = stats.W ./ sqrt(sum(stats.W.^2,1));  % normalized weight
            p = size(XL,1);
            sumSq = sum(XS.^2,1) .* sum(yl.^2,1);
            vipScore(:, ifold) = sqrt( p * sum(sumSq .* (W0.^2), 2) ./ sum(sumSq, 2) );
        end % folds

        % store average fold betas and VIPs for this repetition
        beta_all(rep,:,:) = squeeze(mean(b, 1));
        all_vip(:,rep)    = mean(vipScore, 2);

        % Compute R² the way you requested (scalar across all parcels & freqs)
        mean_obs = mean(Y_raw(:));  % scalar
        r2 = 1 - ( mean( (Y_raw(:) - Y1(:)).^2 ) / mean( (Y_raw(:) - mean_obs).^2 ) );
        R2_all(rep) = r2;

        % RMSE per parcel (across frequencies)
        rmse_parcel(:,rep) = sqrt( mean( (Y1 - Y_raw).^2, 2 ) );
        rmse_matrix(rep,:,:) = sqrt( (Y1 - Y_raw).^2 );

        % Rank correlations per parcel across frequencies 
        r = corr(Y_raw', Y1');
        r_rank = tiedrank(r);
        parcel_rank(:,rep) = diag(r_rank);
        medrank(rep) = median(diag(r_rank));
    end % reps

    % Summary metrics
    average_R2         = mean(R2_all);
    average_rank       = mean(medrank);
    average_parcelrank = mean(parcel_rank, 2);
    average_vip        = mean(all_vip, 2);
    average_weights    = squeeze(mean(beta_all(:, 2:end, :), 1)); % drop intercept

    % Sort maps by VIP
    [VIPs, I] = sort(average_vip, 'descend');
    selected = all_labels(I,:);
    maps_selected = selected;
    maps_selected.scores = VIPs;
    maps_selected.predvariable = repmat("Pow_10ncomp", length(VIPs), 1);

    % Save outputs for subject s
    writetable(maps_selected, fullfile(outputDir, sprintf('maps_allselected_Powadjusted_zscore_10ncomp_noleakCorr_sub%03d.xlsx', s)));

    save(fullfile(outputDir, sprintf('Powadjusted_allselected_zscore_10ncomp_pls_noleakgCorr_sub%03d.mat', s)), ...
        'average_R2', 'average_rank', 'average_vip', 'average_weights', ...
        'beta_all', 'maps_selected', 'parcel_rank', 'rmse_matrix', 'rmse_parcel');

    save(fullfile(outputDir, sprintf('orig_pred_powspctrm_allselected_zscore_10ncomp_noleakgCorr_sub%03d.mat', s)), ...
        'Y_raw', 'Y1');

    % Spin test 
    [p_spin, r_dist, permid] = spin_test_old(rand(1,200), rand(1,200), 'parcellation_name', 'schaefer_200');
    cv = cvpartition(nParc, "Kfold", kfold);
    rsqsh = zeros(1, 1000);
    spin_parcel_rank = zeros(nParc, 1000);
    spin_medrank = zeros(1, 1000);

    for k = 1:1000
        X_perm = X_full(permid(:,k), :);
        Y1_perm = zeros(size(Y_raw));
        for ifold = 1:kfold
            trIdx = cv.training(ifold); teIdx = cv.test(ifold);
            % standardize permuted X by training parcels (again avoid leakage)
            Xtr_perm = X_perm(trIdx, :);
            muX = mean(Xtr_perm,1); sdX = std(Xtr_perm,0,1); sdX(sdX==0)=1;
            Xte_perm = (X_perm(teIdx,:) - muX) ./ sdX;
            Xtr_perm_z = (Xtr_perm - muX) ./ sdX;

            [xl, yl, XS, YS, beta, PCTVAR, mse, stat] = plsregress(Xtr_perm_z, Y_raw(trIdx,:), ncomp);
            Y1_perm(teIdx,:) = [ones(sum(teIdx),1) Xte_perm] * beta;
        end
        mean_obs = mean(Y_raw(:));
        rsqsh(k) = 1 - ( mean( (Y_raw(:) - Y1_perm(:)).^2 ) / mean( (Y_raw(:) - mean_obs).^2 ) );
        r = corr(Y_raw', Y1_perm');
        r_rank = tiedrank(r);
        spin_parcel_rank(:,k) = diag(r_rank);
        spin_medrank(k) = median(diag(r_rank));
    end

    results_single_spin.rsqsh_spin = rsqsh;
    spin_99prctile      = prctile(rsqsh, 99);
    spin_99prctile_rank = prctile(spin_medrank, 99);
    pval_single         = mean(abs(rsqsh) > abs(average_R2));
    pval_medrank        = mean(abs(spin_medrank) > abs(average_rank));

    save(fullfile(outputDir, sprintf('Powadjusted_allselected_zscore_10ncomp_pls_noleakgCorr_sub%03d.mat', s)), ...
        'results_single_spin', 'pval_single', 'pval_medrank', 'spin_99prctile_rank', '-append');

    clear beta_all all_vip R2_all rmse_parcel rmse_matrix parcel_rank medrank ...
          maps_selected results_single_spin rsqsh spin_parcel_rank spin_medrank
end

disp('Done.');
