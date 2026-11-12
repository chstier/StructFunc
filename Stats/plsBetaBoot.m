function B = plsBetaBoot(X, Y, ncomp, nBoot, seed, baselinetype)
% PLSBETABOOT  PLS regression coefficients (maps x freqs/lags) with parcel bootstrap.
% The model is fitted on all parcels; parcels are then resampled with replacement
% and the model is refitted (nBoot times).
% Normalization identical to cvR2_pls.m: X standardized per map,
% Y centred ('absolute') or z-scored ('zscore') per frequency bin / lag.
% Note: parcels are treated as independent, so CIs are somewhat optimistic.
% returns struct B with fields
%   beta  : coefficients from the full-data model (maps x freqs/lags)
%   lo,hi : 2.5 and 97.5 percentiles across bootstrap resamples
%   pSign : bootstrap sign test, 2*min(P(beta>0), P(beta<0))

[n, p] = size(X);  q = size(Y, 2);
B.beta = plsBeta(X, Y, ncomp, baselinetype);
bb = zeros(nBoot, p*q, 'single');
parfor b = 1:nBoot
    s   = RandStream('mt19937ar', 'Seed', seed + b);
    idx = randi(s, n, n, 1);
    tmp = plsBeta(X(idx,:), Y(idx,:), ncomp, baselinetype);
    bb(b,:) = tmp(:)';
end
B.lo    = reshape(prctile(bb,  2.5, 1), p, q);
B.hi    = reshape(prctile(bb, 97.5, 1), p, q);
B.pSign = reshape(2*min(mean(bb > 0, 1), mean(bb < 0, 1)), p, q);
B.baselinetype = baselinetype;
end

function beta = plsBeta(X, Y, ncomp, baselinetype)
% standardize X - as in cvR2_pls
sigmaX = std(X, 0, 1);  sigmaX(sigmaX==0) = 1;
X = (X - mean(X,1)) ./ sigmaX;
% normalize Y according to baselinetype - as in cvR2_pls
switch lower(baselinetype)
    case 'zscore'
        sigmaY = std(Y, 0, 1);  sigmaY(sigmaY==0) = 1;
    case 'absolute'
        sigmaY = ones(1, size(Y,2));
    otherwise
        error('baselinetype must be ''zscore'' or ''absolute''.');
end
Y = (Y - mean(Y,1)) ./ sigmaY;
[~, ~, ~, ~, BETA] = plsregress(X, Y, ncomp);
beta = BETA(2:end, :);                                   % drop intercept row
end