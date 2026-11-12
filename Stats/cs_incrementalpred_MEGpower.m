%% This script does forward selection of neurochemical/structural maps
% based on best CV R^2 at each step (cvR2_pls.m, same as main analyses)
% Forward selection: start from best single map, add maps iteratively
% based on highest accuracy improvement

% CS 2025

%% Settings
metric  = 'power'; 
kfold   = 10;
reps_cv = 50;
seed    = 1234;      % same CV partitions for every candidate

%% Load neurochemical/structual maps
load('allmaps_concat2_selected.mat')
X_raw = double(allmaps);   % original predictors (200 x nMaps), z-scored inside CV

%% Load MEG data
switch metric
    case 'power'   % adjusted for age, age2, sex, tiv
        load('adjusted_MEGpower_orig.mat')
        gpowav = squeeze(mean(all_adj, 1));                       % parcel x freq
        Y_raw  = log(gpowav ./ repmat(mean(gpowav,2),1,60));      % (200 x 60)
        baselinetype = 'zscore';
    otherwise
        error('metric must be ''power''');
end

%% Forward selection
max_feats = size(X_raw,2);  % total number of maps

selected_maps  = [];
remaining_maps = 1:max_feats;

cv_R2_list    = zeros(1, max_feats);    % CV R^2 at each step
ncomp_list    = zeros(1, max_feats);    % number of PLS components
added_maps    = zeros(1, max_feats);    % which map was added
added_labels  = strings(1, max_feats);
selected_vips = cell(1, max_feats);     % CV-averaged VIPs for chosen set

if isempty(gcp('nocreate')), parpool; end
fprintf('Starting greedy forward selection (%s, accuracy-based)...\n', metric);

for step = 1:max_feats
    ncomp = ceil((numel(selected_maps) + 1) / 5);   % increase components gradually
    nCand = numel(remaining_maps);
    r2_cand  = zeros(1, nCand); 
    vip_cand = cell(1, nCand);

    % Try adding each remaining map to the selected ones
    parfor c = 1:nCand
        candidate = [selected_maps, remaining_maps(c)];
        rng(seed);                                  % identical partitions for all candidates
        [r2_cand(c), vip_cand{c}] = cvR2_pls(candidate, X_raw, Y_raw, ncomp, kfold, reps_cv, baselinetype);
    end

    % Pick the best candidate
    [best_r2, c_best] = max(r2_cand);
    best_map = remaining_maps(c_best);

    selected_maps = [selected_maps, best_map];
    remaining_maps(remaining_maps == best_map) = [];

    cv_R2_list(step)    = best_r2;
    ncomp_list(step)    = ncomp;
    added_maps(step)    = best_map;
    selected_vips{step} = vip_cand{c_best};         % VIPs (length = numel(selected_maps))
    if exist('all_labels', 'var')
        added_labels(step) = string(all_labels.annotation(best_map, :));
    end

    fprintf('Step %d: Added map %d (%s) -> CV R^2 = %.4f (ncomp = %d), remaining = %d\n', ...
        step, best_map, added_labels(step), best_r2, ncomp, numel(remaining_maps));
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
title(sprintf('Greedy Forward Selection (%s): CV R^2 vs. Added Maps', metric));
grid on;

%% Save results
T = table(added_maps(:), added_labels(:), cv_R2_list(:), ncomp_list(:), ...
          'VariableNames', {'MapIndex', 'MapLabel', 'CV_R2', 'NumComponents'});
outFile = sprintf('greedy_forward_selection_results_%s.csv', metric);
writetable(T, outFile);
save(sprintf('greedy_forward_selection_%s.mat', metric), 'T', 'selected_vips', 'metric');
fprintf('Results saved to %s\n', outFile);
