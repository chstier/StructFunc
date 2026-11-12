function [R2_pooled, mean_vip, mean_rank, details] = distcv_pls(feats, X_raw, Y_raw, ncomp, coords, trainFraction, baselinetype)

    % Distance-dependent cross-validation of PLS regression.
    % Counterpart of cvR2_pls.m with identical standardization, R^2 and VIP.
    % Instead of random k-fold partitions, each parcel serves once as seed:
    % the nearest trainFraction of parcels (Euclidean centroid distance) are
    % used for training, the most distant parcels exclusively for testing.
    % R^2 = 1 - sum((Y - Yhat).^2) / sum((Y - mean_train).^2), pooled over all
    % splits, held-out parcels and frequencies/lags (original Y units).
    % feats, X_raw, Y_raw, ncomp, baselinetype : as in cvR2_pls.m
    % coords        : parcels x 3 centroid coordinates (same order as X, Y)
    % trainFraction : fraction of nearest parcels used for training (0.75)
    % returns:
    % R2_pooled : scalar, pooled distance-CV R^2
    % mean_vip  : vector (length = numel(feats)), average VIP across splits
    % mean_rank : median parcel identification rank within the test set
    %             (of nTest parcels), mean across splits
    % details   : (optional) struct with per-split results
    % Because splits overlap, their spread is a robustness range, not a CI.

    nObs   = size(X_raw,1);
    p      = numel(feats);
    nY     = size(Y_raw,2);
    nTrain = floor(trainFraction * nObs);
    nTest  = nObs - nTrain;
    D      = pdist2(coords, coords, 'euclidean');

    R2_split      = nan(nObs,1);
    rho_split     = nan(nObs,1);
    medrank_split = nan(nObs,1);
    vip_split     = zeros(p, nObs);
    beta_split    = zeros(nObs, p+1, nY);           % intercept included (as beta_all)
    predSum   = zeros(nObs, nY);
    predCount = zeros(nObs, 1);
    trainIdx  = cell(nObs,1);
    testIdx   = cell(nObs,1);
    SSE_total = 0;  SST_total = 0;  SSTtest_total = 0;

    for s = 1:nObs
        % nearest parcels train, most distant parcels test
        [~, order] = sort(D(s,:), 'ascend');
        train_idx = order(1:nTrain);
        test_idx  = order(nTrain+1:end);
        trainIdx{s} = train_idx;  testIdx{s} = test_idx;

        % standardize X (training stats only) - as in cvR2_pls
        X_train_raw = X_raw(train_idx, feats);
        X_test_raw  = X_raw(test_idx, feats);
        muX = mean(X_train_raw, 1);
        sigmaX = std(X_train_raw, 0, 1);
        sigmaX(sigmaX==0) = 1;
        Xtrain = (X_train_raw - muX) ./ sigmaX;
        Xtest  = (X_test_raw  - muX) ./ sigmaX;

        % normalize Y according to baselinetype (training parcels only) - as in cvR2_pls
        Ytrain_raw = Y_raw(train_idx, :);
        Ytest_raw  = Y_raw(test_idx, :);
        switch lower(baselinetype)
            case 'zscore'
                meanFreq = mean(Ytrain_raw, 1);
                stdFreq  = std(Ytrain_raw, 0, 1);
                stdFreq(stdFreq==0) = 1;
            case 'absolute'
                meanFreq = mean(Ytrain_raw, 1);
                stdFreq  = ones(1, nY);
            otherwise
                error('baselinetype must be ''zscore'' or ''absolute''.');
        end
        Ytrain = (Ytrain_raw - meanFreq) ./ stdFreq;

        % fit PLS on near parcels, predict distant parcels, inverse-transform
        [XL, yl, XS, YS, beta, PCTVAR, mse, stats] = plsregress(Xtrain, Ytrain, ncomp);
        beta_split(s,:,:) = beta;
        Ypred = ([ones(nTest,1), Xtest] * beta) .* stdFreq + meanFreq;

        % prediction error vs. training-mean baseline (main) and test-mean baseline (stricter)
        SSE     = sum((Ytest_raw - Ypred).^2, 'all');
        SST     = sum((Ytest_raw - meanFreq).^2, 'all');
        SSTtest = sum((Ytest_raw - mean(Ytest_raw,1)).^2, 'all');
        R2_split(s)   = 1 - SSE/SST;
        SSE_total     = SSE_total + SSE;
        SST_total     = SST_total + SST;
        SSTtest_total = SSTtest_total + SSTtest;

        rho_split(s) = corr(Ytest_raw(:), Ypred(:), 'Type', 'Spearman', 'Rows', 'complete');

        % median identification rank within the test set (as in cvR2_pls)
        r = corr(Ytest_raw', Ypred', 'Rows', 'pairwise');
        r_rank = zeros(nTest);
        for k0 = 1:nTest
            r_rank(k0,:) = tiedrank(r(k0,:));
        end
        medrank_split(s) = median(diag(r_rank), 'omitnan');

        % VIP from the training model (as in cvR2_pls)
        W0 = stats.W ./ sqrt(sum(stats.W.^2, 1));
        sumSq = sum(XS.^2, 1) .* sum(yl.^2, 1);
        denom = sum(sumSq);
        if denom > 0
            vip_split(:,s) = sqrt( p * ( sum( ( (W0.^2) .* repmat(sumSq, p, 1) ), 2 ) ) / denom );
        end

        % accumulate predictions (parcels appear in several test sets)
        predSum(test_idx,:) = predSum(test_idx,:) + Ypred;
        predCount(test_idx) = predCount(test_idx) + 1;
    end

    R2_pooled = 1 - SSE_total / SST_total;
    mean_vip  = mean(vip_split, 2);
    mean_rank = mean(medrank_split, 'omitnan');

    if nargout > 3
        details = struct( ...
            'R2_split', R2_split, ...
            'R2_pooled_testmean', 1 - SSE_total / SSTtest_total, ...
            'mean_R2_split', mean(R2_split, 'omitnan'), ...
            'median_R2_split', median(R2_split, 'omitnan'), ...
            'rho_split', rho_split, 'mean_rho', mean(rho_split, 'omitnan'), ...
            'medrank_split', medrank_split, 'vip_split', vip_split, ...
            'beta_split', beta_split, ...
            'Ypred_average', predSum ./ predCount, ...   % NaN = never tested
            'predictionCount', predCount, ...
            'nTrain', nTrain, 'nTest', nTest);
        details.trainIdx = trainIdx;                     % cell arrays added separately
        details.testIdx  = testIdx;
    end
end