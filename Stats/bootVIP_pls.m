function vipBoot = bootVIP_pls(feats, X_raw, Y_raw, ncomp, baselinetype, nBoot, seed)
    % Parcel bootstrap of VIP scores (95% CIs for the VIP ranking).
    % Parcels are resampled with replacement; X and Y are normalized within
    % each bootstrap sample as in cvR2_pls.m; PLS is refitted.
    % Note: parcels are treated as independent, so CIs are somewhat optimistic.
    % returns: numel(feats) x nBoot matrix of VIP scores

    X_raw = X_raw(:, feats);
    [nObs, p] = size(X_raw);
    rng(seed);
    idx = randi(nObs, nObs, nBoot);          % fixed resamples (reproducible)
    vipBoot = zeros(p, nBoot);

    parfor b = 1:nBoot
        Xb = X_raw(idx(:,b), :);
        Yb = Y_raw(idx(:,b), :);

        % standardize X - as in cvR2_pls
        sigmaX = std(Xb, 0, 1);  sigmaX(sigmaX==0) = 1;
        Xb = (Xb - mean(Xb,1)) ./ sigmaX;

        % normalize Y according to baselinetype - as in cvR2_pls
        switch lower(baselinetype)
            case 'zscore'
                stdFreq = std(Yb, 0, 1);  stdFreq(stdFreq==0) = 1;
            case 'absolute'
                stdFreq = ones(1, size(Yb,2));
            otherwise
                error('baselinetype must be ''zscore'' or ''absolute''.');
        end
        Yb = (Yb - mean(Yb,1)) ./ stdFreq;

        [~, yl, XS, ~, ~, ~, ~, stats] = plsregress(Xb, Yb, ncomp);

        % VIP - as in cvR2_pls
        W0 = stats.W ./ sqrt(sum(stats.W.^2, 1));
        sumSq = sum(XS.^2, 1) .* sum(yl.^2, 1);
        denom = sum(sumSq);
        if denom > 0
            vipBoot(:,b) = sqrt( p * ( sum( ( (W0.^2) .* repmat(sumSq, p, 1) ), 2 ) ) / denom );
        end
    end
end