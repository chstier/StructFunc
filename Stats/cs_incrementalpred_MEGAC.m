%% This script does forward selection of neurochemical/structural maps
% based on best CV R² at each step

% CS 2025

% Forward selection: start from best single map, add maps iteratively
% based on highest accuracy improvement

%% Load respective data
% Load neurochemical/structual maps
load('allmaps_concat2_selected.mat')

% No zscoring globally, but in CV loop - keep raw maps
X_raw = double(allmaps);   % original predictors (200 x nMaps)

% Get MEG autocorrelation data
% adjusted for age, sex, icv
load('adjusted_MEGACall_sorted.mat')

% get average values (parcel x freq)
gACav = squeeze(mean(all_adj, 1));
Y_raw = gACav;

Y = Y_raw;

% % Choose normalization method to compare (will be applied within-fold using training stats)
% baselinetype = 'absolute'; 


kfold = 10;
reps_cv = 50;
max_feats = size(X_raw,2);  % total number of maps

selected_maps = [];
remaining_maps = 1:max_feats;

cv_R2_list = zeros(1, max_feats);   % CV R^2 at each step
ncomp_list = zeros(1, max_feats);   % number of PLS components
added_maps = zeros(1, max_feats);   % which map was added
added_labels = strings(1, max_feats);  % optional
selected_vips = cell(1, max_feats);  % store CV-averaged VIPs for chosen set (optional)

fprintf('Starting greedy forward selection (accuracy-based)...\n');

for step = 1:max_feats
    best_r2 = -Inf;
    best_map = NaN;
    best_vip = [];

    % Try adding each remaining map to the selected ones
    for m = remaining_maps
        candidate = [selected_maps, m];
        ncomp = ceil(numel(candidate) / 5);  % increase components gradually

        % Now cvR2_pls returns both r2 and mean_vip (VIPs averaged across CV)
        % Pass X_raw and baselinetype
        [r2, mean_vip] = cvR2_pls_ACall(candidate, X_raw, Y_raw, ncomp, kfold, reps_cv);

        if r2 > best_r2
            best_r2 = r2;
            best_map = m;
            best_vip = mean_vip;   % store vip for the best candidate (no extra fit)
        end
    end

    % Update lists with the best candidate
    selected_maps = [selected_maps, best_map];
    remaining_maps(remaining_maps == best_map) = [];

    cv_R2_list(step) = best_r2;
    ncomp_list(step) = ceil(numel(selected_maps) / 5);
    added_maps(step) = best_map;
    selected_vips{step} = best_vip;   % VIPs (length = numel(selected_maps))

    if exist('all_labels', 'var')
        added_labels(step) = string(all_labels.annotation(best_map, :));
    end

    fprintf('Step %d: Added map %d ? CV R^2 = %.4f (ncomp = %d)\n', ...
        step, best_map, best_r2, ncomp_list(step));
    fprintf('Step %d: Remaining maps = %d\n', step, numel(remaining_maps));
end

%% Plot R^2 and number of components
figure;
yyaxis left
plot(1:max_feats, cv_R2_list, '-o', 'LineWidth', 2);
ylabel('Cross-validated R^2');
ylim([min(0, min(cv_R2_list)), max(cv_R2_list)*1.1]);

yyaxis right
plot(1:max_feats, ncomp_list, '--x', 'Color', [0.4 0.4 0.4]);
ylabel('# PLS Components');

xlabel('Number of Maps Added');
title('Forward Selection: CV R^2 vs. Added Maps');
grid on;

%% Save results
T = table(added_maps(:), ...
          added_labels(:), ...
          cv_R2_list(:), ...
          ncomp_list(:), ...
          'VariableNames', {'MapIndex', 'MapLabel', 'CV_R2', 'NumComponents'});

writetable(T, 'forward_selection_results.csv');
fprintf('Results saved to forward_selection_results.csv\n');