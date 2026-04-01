%% This script loads all suitable maps from the ENIGMA toolbox
%% CS Dec 2024

%% Make a selected set of features
% exclude cognitive and func connectivity gradients

data_dir = '/home/path/data/';
maps = {'struc_func_hubs.mat', 'big_brain.mat', 'rec_maps_schaefer200.mat', 'surface_maps_schaefer200.mat'};

for i = 1:length(maps)
 file = fullfile(data_dir, maps{i});
 load(file)
end

allmaps = [fc_ctx_dc, sc_ctx_dc, bb_gradient, bb_moments, squeeze(neuromaps_receptors), neuromaps_mic];

receptor_labels = cellstr(reclabels);
micro_labels = cellstr(microsclabels);

all_labels = {'func_hubs_hcp', 'struct_hubs_hcp', 'bb_gradient', 'bb_moment1', 'bb_moment2', 'bb_moment3', ...
 'bb_moment4', 'bb_moment5', receptor_labels{:}, micro_labels{:}};

% get labels for selection of features
selected_labels = readtable('/home/path/maps_classes_selected.csv');

idx = find(ismember(string(all_labels), string(selected_labels.annotation)));
selected_maps = allmaps(:, idx);

% allmaps(:,66:76) = []; % exclude func gradients fMRI
% allmaps(:,37) = []; % exclude neurosynth
% all_labels(:,66:76) = []; % exclude func gradients fMRI
% all_labels(:,37) = []; % exclude neurosynth

% get labels for pca neurotransmitter
pca_labels = readtable('/home/path/maps_pca_transm.csv', 'ReadVariableNames', false);

% extract single neurotransmitter maps for pca
idx = find(ismember(string(all_labels), string(pca_labels.Var1)));
trans_maps = allmaps(:, idx);
trans_maps = double(trans_maps);
[coeff, score, ~, ~, explained] = pca(trans_maps);
neurotransmpc1 = score(:,1);  % First principal component scores for parcels
first_pc_loadings = coeff(:,1);

% Extend and annotate different classes of information
% summary metrics
summ = {'neurotransmpc1', 'abagen-genepc1', 'sydnor2021-SAaxis'}; % exclude neurosynth here
iid = find(ismember(string(all_labels), string(summ)));
summary_data = allmaps(:, iid);
summary_data = [neurotransmpc1, summary_data];   
sum_labels = summ;
save('maps_summary.mat', 'summary_data', 'summ') % save summary data

allmaps = [selected_maps, neurotransmpc1]; % attach pc1 of neurotransmitters and save
selected_labels.annotation{end+1,:} = 'neurotransmpc1';
selected_labels.class{end,:} = 'neurotransmitter - summary';
selected_labels.metric{end,:} = 'pc1 neurotransmitter';
selected_labels.tags{end,:} = 'ntpc1';
selected_labels.modality{end,:} = "";
all_labels = selected_labels;

writetable(all_labels,'allmaps_labels_full_selected.txt','WriteVariableNames', false)
save('allmaps_concat2_selected.mat', 'allmaps', 'all_labels')


% Microstructure
clear iid
iid = find(ismember(string(all_labels.class), string('microstructural')));
microstruct_data = allmaps(:, iid);
microstruct_labels = all_labels(iid,:);

% Metabolism
clear iid
iid = find(ismember(string(all_labels.class), string('metabolism')));
metab_data = allmaps(:, iid);
metab_labels = all_labels(iid,:);

% Cortical expension
clear iid
iid = find(ismember(string(all_labels.class), string('cortical expansion')));
corticalexp_data = allmaps(:, iid);
corticalexp_labels = all_labels(iid,:);

% Network metrics(structural and functional)
clear iid
iid = find(ismember(string(all_labels.class), string('network-metrics')));
network_data = allmaps(:, iid);
network_labels = all_labels(iid,:);

% Neurotransmitter
clear iid
iid = find(contains(string(all_labels.class), 'neurotransmitter'));
neurotr_data = allmaps(:, iid);
neurotr_labels = all_labels(iid,:);
neurotr_data(:,end) = [];
neurotr_labels(end,:) = [];

% save classes
save('maps_classes_selected.mat', 'summary_data', 'sum_labels', 'microstruct_data', 'microstruct_labels', 'metab_data', 'metab_labels', ...
 'corticalexp_data', 'corticalexp_labels', 'network_data', 'network_labels', 'neurotr_data', 'neurotr_labels')
