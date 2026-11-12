function S = sens_load(resFile, mapFile)
% Load predictor maps, saved age map/results, spins and aligned metadata
if nargin < 2, mapFile = 'allmaps_concat2_selected.mat'; end
Mf       = load(mapFile, 'allmaps', 'all_labels');
S.X      = double(Mf.allmaps);                       % parcels x maps
S.labels = string(Mf.all_labels{:,1});
S.res    = load(resFile);
S.ageMap = S.res.ageMap;                             % parcels x freqs (or lags)

if isfield(S.res, 'labels')
    assert(isequal(string(S.res.labels(:)), S.labels(:)), 'Map order differs from saved results.');
end
if isfield(S.res, 'permid')
    S.permid = S.res.permid;
else
    warning('permid not stored in %s - regenerating spins.', resFile);
    n = size(S.X, 1);
    [~, ~, S.permid] = spin_test_old(rand(1,n), rand(1,n), 'parcellation_name', 'schaefer_200');
end

% align metadata to data label order: alias first, then normalised match
M     = map_metadata();
A     = map_label_alias();
key   = S.labels(:);
hasA  = isKey(A, cellstr(key));
key(hasA) = string(values(A, cellstr(key(hasA))));

norm = @(s) lower(regexprep(string(s), '[\s_\-]', ''));
[tf, idx] = ismember(norm(key), norm(M.label));
if ~all(tf), error('No metadata for: %s', strjoin(S.labels(~tf), ', ')); end
if numel(unique(idx)) ~= numel(idx)
    [~, i1] = unique(idx); dup = setdiff(1:numel(idx), i1);
    error('Duplicate metadata matches for: %s', strjoin(S.labels(dup), ', '));
end
unused = setdiff(1:height(M), idx);
if ~isempty(unused)
    warning('Metadata rows not used: %s', strjoin(M.label(unused), ', '));
end
S.meta = M(idx, :);
S.meta.label_data = S.labels(:);
S.short = S.meta.label;                              % Table 1 names for plots/tables
end