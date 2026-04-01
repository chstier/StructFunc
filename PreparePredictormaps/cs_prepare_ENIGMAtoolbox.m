%% This script loads all suitable maps from the ENIGMA toolbox
%% CS Nov 2024

% set paths
addpath(genpath('/home/path/'))

% load annotation keys
load('annot_keys_Schaefer200.mat')
% load original annotation key suma
load('Schaefer200_suma_reordered.mat')


%% BigBrains moments
% moment
load('/path/ENIGMA/enigmatoolbox/histology/bb_moments_schaefer_200.csv');
bb_moments = bb_moments_schaefer_200';

% create vector with NaNs
for p = 1:length(bb_moments_schaefer_200(:,1))
 plot_dat(:,p) = NaN(size(suma_all.annot_key2{1,2}));
end

for b = 1:length(bb_moments_schaefer_200(:,1))
 % map back to our SUMA/schafer order
 for i = 1:200
  v = bb_moments_schaefer_200(b,i);
  id = key_rev_cortical(i);
  plot_dat(id,b) = v;
 end
end

BB_labels = {'1', '2', '3', '4', '5'};

% load maps for plotting
suma_schaefer = load([pwd '/Schaefer200_suma_reordered.mat'], 'suma_all');

for p = 1:length(bb_moments_schaefer_200(:,1))
 
 file_name = cell2mat(BB_labels(p));
 data2plot = plot_dat(:,p);
 keydata = [data2plot suma_schaefer.suma_all.annot_key{1,1}];
                    vertex_data = zeros(2338,2);

                    for v = 1:length(suma_schaefer.suma_all.annot)
                        keynr = suma_schaefer.suma_all.annot(v);
                        if keynr == 0
                            vertex_data(v,1) = 0;
                        else
                            rownr = find(keydata(:,2) == keynr);
                            log_p_v = keydata(rownr,1);
                            vertex_data(v,1) = log_p_v;
                            logp = vertex_data(:,1);
                        end
                    end  

 % Plot
  opt=[];
  opt.title = ['BigBrain moments'];
  opt.output = [results_dir '/' file_name '_bbmoments.png'];
  opt.per_hemi=1;
  opt.per_cortex=1;
  opt.rot=[90 0 ; -90 0];
  opt.thresh=0; % keep threshold here for NaNs in the data
  opt.clim=[min(data2plot) max(data2plot)];
  opt.opathresh = 0.00000000001;
  opt.colormap='hot';
  opt.colorbar='hot';
  opt.scale=1;

  hFig = cs_nmri_plot_surface_suma(suma_all, logp, opt);
end

%% BigBrains gradients

load('/home/uni10/nmri/projects/cstier/aging_multiv/ENIGMA/enigmatoolbox/histology/bb_gradient_schaefer_200.csv');

% create vector with NaNs
for p = 1:length(bb_gradient_schaefer_200(:,1))
 plot_dat(:,p) = NaN(size(suma_all.annot_key2{1,2}));
end

for b = 1:length(bb_gradient_schaefer_200(:,1))
 % map back to our SUMA/schafer order
 for i = 1:200
  v = bb_gradient_schaefer_200(b,i);
  id = key_rev_cortical(i);
  plot_dat(id,b) = v;
 end
end

bb_gradient = bb_gradient_schaefer_200';
BB_labels = {''};

% load maps for plotting
suma_schaefer = load([pwd '/Schaefer200_suma_reordered.mat'], 'suma_all');

for p = 1:length(bb_gradient_schaefer_200(:,1))
 
 file_name = cell2mat(BB_labels(p));
 data2plot = plot_dat(:,p);
 keydata = [data2plot suma_schaefer.suma_all.annot_key{1,1}];
                    vertex_data = zeros(2338,2);

                    for v = 1:length(suma_schaefer.suma_all.annot)
                        keynr = suma_schaefer.suma_all.annot(v);
                        if keynr == 0
                            vertex_data(v,1) = 0;
                        else
                            rownr = find(keydata(:,2) == keynr);
                            log_p_v = keydata(rownr,1);
                            vertex_data(v,1) = log_p_v;
                            logp = vertex_data(:,1);
                        end
                    end  

 % Plot
  opt=[];
  opt.title = ['BigBrain gradients'];
  opt.output = [results_dir '/' file_name 'bbgradient.png'];
  opt.per_hemi=1;
  opt.per_cortex=1;
  opt.rot=[90 0 ; -90 0];
  opt.thresh=0; % keep threshold here for NaNs in the data
  opt.clim=[min(data2plot) max(data2plot)];
  opt.opathresh = 0.00000000001;
  opt.colormap='parula';
  opt.colorbar='parula';
  opt.scale=1;

  hFig = cs_nmri_plot_surface_suma(suma_all, logp, opt);
end

%% save BigBrain variables
save('big_brain.mat', 'bb_gradient', 'bb_moments')

%% Structual hubs
% Load cortico-cortical connectivity map
[sc_ctx, sc_ctx_labels, ~, ~] = load_sc('schaefer_200');
[~, ~, sc_sctx, sc_sctx_labels] = load_sc('schaefer_200');

% compute weighted degree centrality measures from connectivity
sc_ctx_dc = sum(sc_ctx)';
sc_sctx_dc = sum(sc_sctx, 2);
sc_sctx_dc([1 8]) = [];

% Plot using Schaefer Parcellation
% attach cortical + subcortical
sc_all = [sc_ctx_dc;sc_sctx_dc];

% create vector with NaNs
% load original annotation key suma
load('Schaefer200_suma_reordered.mat')
plot_dat = NaN(size(suma_all.annot_key2{1,2}));

% map back to our SUMA/schafer order
for i = 1:length(sc_all)
 v = sc_all(i);
 id = key_rev_all(i);
 plot_dat(id) = v;
end

results_dir = pwd;
 
% load maps for plotting
load([pwd '/conf/suma-all-fsaverage-10.mat'],'suma_all')
suma_schaefer = load([pwd '/Schaefer200_suma_reordered.mat'], 'suma_all');

file_name = 'struct';
 
keydata = [plot_dat suma_schaefer.suma_all.annot_key{1,1}];
                   vertex_data = zeros(2338,2);

                   for v = 1:length(suma_schaefer.suma_all.annot)
                       keynr = suma_schaefer.suma_all.annot(v);
                       if keynr == 0
                           vertex_data(v,1) = 0;
                       else
                           rownr = find(keydata(:,2) == keynr);
                           log_p_v = keydata(rownr,1);
                           vertex_data(v,1) = log_p_v;
                           logp = vertex_data(:,1);
                       end
                   end  

% Plot
 opt=[];
 opt.title = ['Cortico-cortical degree centrality (hubs)'];
 opt.output = [results_dir '/' file_name '_hubs.png'];
 opt.per_hemi=1;
 opt.per_cortex=1;
 opt.rot=[90 0 ; -90 0];
 opt.thresh=0; % keep threshold here for NaNs in the data
 opt.clim=[min(plot_dat) max(plot_dat)];
 opt.opathresh = 0.00000000001;
 opt.colormap='parula';
 opt.colorbar='parula';
 opt.scale=1;

 hFig = cs_nmri_plot_surface_suma(suma_all, logp, opt);

%% Functional hubs
% Load cortico-cortical connectivity map
[fc_ctx, fc_ctx_labels, ~, ~] = load_fc('schaefer_200');
[~, ~, fc_sctx, fc_sctx_labels] = load_fc('schaefer_200');

% compute weighted degree centrality measures from connectivity
fc_ctx_dc = sum(fc_ctx)';
fc_sctx_dc = sum(fc_sctx, 2);
fc_sctx_dc([1 8]) = [];

% Plot using Schaefer Parcellation
% attach cortical + subcortical
fc_all = [fc_ctx_dc;fc_sctx_dc];

% Plot using Schaefer Parcellation
% create vector with NaNs
% load original annotation key suma
load('Schaefer200_suma_reordered.mat')
plot_dat = NaN(size(suma_all.annot_key2{1,2}));

% map back to our SUMA/schafer order
for i = 1:length(fc_all)
 v = fc_all(i);
 id = key_rev_all(i);
 plot_dat(id) = v;
end

results_dir = pwd;
 
% load maps for plotting
load([pwd '/conf/suma-all-fsaverage-10.mat'],'suma_all')
suma_schaefer = load([pwd '/Schaefer200_suma_reordered.mat'], 'suma_all');

file_name = 'func';
 
keydata = [plot_dat suma_schaefer.suma_all.annot_key{1,1}];
                   vertex_data = zeros(2338,2);

                   for v = 1:length(suma_schaefer.suma_all.annot)
                       keynr = suma_schaefer.suma_all.annot(v);
                       if keynr == 0
                           vertex_data(v,1) = 0;
                       else
                           rownr = find(keydata(:,2) == keynr);
                           log_p_v = keydata(rownr,1);
                           vertex_data(v,1) = log_p_v;
                           logp = vertex_data(:,1);
                       end
                   end  

% Plot
 opt=[];
 opt.title = ['Cortico-cortical degree centrality (hubs)'];
 opt.output = [results_dir '/' file_name '_hubs.png'];
 opt.per_hemi=1;
 opt.per_cortex=1;
 opt.rot=[90 0 ; -90 0];
 opt.thresh=0; % keep threshold here for NaNs in the data
 opt.clim=[min(plot_dat) max(plot_dat)];
 opt.opathresh = 0.00000000001;
 opt.colormap='hot';
 opt.colorbar='hot';
 opt.scale=1;

 hFig = cs_nmri_plot_surface_suma(suma_all, logp, opt);
 
%% save structural and functional hub data
save('struc_func_hubs.mat', 'sc_ctx_dc', 'sc_sctx_dc', 'fc_ctx_dc', 'fc_sctx_dc')

%% summary statistics disorders
% Yet to come - needs remapping from Desikan-Kiliany to SUMA and
% Schaefer200

