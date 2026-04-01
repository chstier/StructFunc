function [mean_r2, mean_vip] = cvR2_pls_ACall(feats, X_raw, Y_raw, ncomp, kfold, reps_cv)

    % Computes cross-validated R^2 using global MSE with per-fold standardization
    % feats : vector of selected feature indices (columns of X)
    % X, Y  : full data matrices (observations x features) and (observations x freqs)
    % ncomp : number of pls components
    % returns: mean R^2 across reps
    % mean_r2 : scalar, mean CV R^2 across reps
    % mean_vip: vector (length = numel(feats)), average VIP across folds & reps

    nObs = size(X_raw,1);
    p = numel(feats);
    r2_vals = zeros(1, reps_cv);

    % preallocate vip accumulator: reps_cv * kfold columns
    vip_accum = zeros(p, reps_cv * kfold);
    vip_col = 0;

       for rep = 1:reps_cv
        cv = cvpartition(nObs, 'KFold', kfold);
        Y_pred = zeros(size(Y_raw));

        for kf = 1:cv.NumTestSets
            train_idx = cv.training(kf);
            test_idx = cv.test(kf);

            % standardize X (training stats only)
            X_train_raw = X_raw(train_idx, feats);
            X_test_raw  = X_raw(test_idx, feats);

            muX = mean(X_train_raw, 1);
            sigmaX = std(X_train_raw, 0, 1);
            sigmaX(sigmaX==0) = 1;  % avoid division by zero

            Xtrain = (X_train_raw - muX) ./ sigmaX;
            Xtest  = (X_test_raw  - muX) ./ sigmaX;

            % normalize Y according to baselinetype(per-frequency) 
            % using training parcels only
            Ytrain_raw = Y_raw(train_idx, :);   % parcels x freqs
            Ytest_raw  = Y_raw(test_idx, :);

%             meanVals = repmat(mean(Ytrain_raw, 1), size(Ytrain_raw,1), 1);

%             switch lower(baselinetype)
%                 case 'absolute'
%                     Ytrain = Ytrain_raw - meanVals;
%                     Ytest_norm  = Ytest_raw  - repmat(mean(Ytrain_raw,1), size(Ytest_raw,1), 1);
% 
%                 otherwise
%                     error('Only ''absolute'' baseline implemented in this function for now.');
%             end
            Ytrain = Ytrain_raw;
            Ytest_norm = Ytest_raw;

            % Fit PLS on training set
            [XL, yl, XS, YS, beta, PCTVAR, mse, stats] = plsregress(Xtrain, Ytrain, ncomp);

            % Predict on test set (in normalized Y-space) and inverse-transform
            nTest = sum(test_idx);
            Ypred_norm = [ones(nTest,1), Xtest] * beta;
            Ypred = Ypred_norm;

%             % inverse absolute baseline
%             meanFreq_train = mean(Ytrain_raw,1);   % 1 x freqs
%             Ypred = Ypred_norm + repmat(meanFreq_train, nTest, 1);

            % store predictions
            Y_pred(test_idx, :) = Ypred;

            % Compute VIP from the training model (no leakage) 
            W0 = stats.W ./ sqrt(sum(stats.W.^2, 1));   % p x ncomp
            sumSq = sum(XS.^2, 1) .* sum(yl.^2, 1);     % 1 x ncomp
            denom = sum(sumSq);
            vip_col = vip_col + 1;
            if denom == 0
                vip_accum(:, vip_col) = zeros(p,1);
            else
                vip_fold = sqrt( p * ( sum( ( (W0.^2) .* repmat(sumSq, p, 1) ), 2 ) ) / denom );
                vip_accum(:, vip_col) = vip_fold;
            end
        end % folds

        % compute R^2 for this rep (in original Y units)
        mean_obs = mean(Y_raw,1);  % use raw Y (original units) for R2 denominator
        r2 = 1-(mean((Y_raw-Y_pred).^2)/(mean((Y_raw-mean_obs).^2)));
        r2_vals(rep) = r2;
        
%         mean_obs = mean(Y, 1);  % 1 x freqs
%         denom_all = mean((Y - repmat(mean_obs, size(Y,1), 1)).^2, 'all');
%         numer_all = mean((Y - Y_pred).^2, 'all');
%         r2_vals(rep) = 1 - numer_all / denom_all;
    end % reps

    mean_r2 = mean(r2_vals);

    if vip_col > 0
        mean_vip = mean(vip_accum(:, 1:vip_col), 2);  % average VIPs across all folds & reps
    else
        mean_vip = zeros(p,1);
    end
end