function [mean_r2, mean_vip, mean_rank, details] = cvR2_pls(feats, X_raw, Y_raw, ncomp, kfold, reps_cv, baselinetype)

    % Computes cross-validated prediction R^2 (1 - NMSE; Scheinost et al., 2019)
    % with per-fold standardization.
    % R^2 = 1 - sum((Y - Yhat).^2) / sum((Y - mean_train).^2), pooled over all
    % held-out parcels and frequencies (original Y units) within a repetition,
    % then averaged over repetitions.
    % feats : vector of selected feature indices (columns of X)
    % X, Y  : full data matrices (observations x features) and (observations x freqs)
    % ncomp : number of pls components
    % returns:
    % mean_r2  : scalar, mean CV R^2 across reps
    % mean_vip : vector (length = numel(feats)), average VIP across folds & reps
    % mean_rank: (optional) median parcel identification rank, mean across reps
    % details  : (optional) struct with per-rep results (only computed if requested):
    %            R2_rep, medrank_rep, parcel_rank, beta_all, rmse_parcel,
    %            rmse_matrix, Y_pred_last

    nObs = size(X_raw,1);
    p = numel(feats);
    r2_vals = zeros(1, reps_cv);
    medrank_vals = zeros(1, reps_cv);                           
    wantDetails  = nargout > 3;                                
    if wantDetails                                             
        nY = size(Y_raw,2);                                    
        parcel_rank = zeros(nObs, reps_cv);                   
        beta_all    = zeros(reps_cv, p+1, nY);                 
        rmse_parcel = zeros(nObs, reps_cv);                    
        rmse_matrix = zeros(reps_cv, nObs, nY);                
    end

    % preallocate vip accumulator: reps_cv * kfold columns
    vip_accum = zeros(p, reps_cv * kfold);
    vip_col = 0;

       for rep = 1:reps_cv
        cv = cvpartition(nObs, 'KFold', kfold);
        Y_pred = zeros(size(Y_raw));
        SSE_rep = 0;
        SST_rep = 0;
        Y_pred_norm = zeros(size(Y_raw));                               
        Y_true_norm = zeros(size(Y_raw));                              
        if wantDetails, b_rep = nan(kfold, p+1, size(Y_raw,2)); end   

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

            switch lower(baselinetype)
                case 'zscore'
                    % compute mean & std across parcels from the training set, per frequency
                    meanFreq = mean(Ytrain_raw, 1);   % 1 x freqs
                    stdFreq  = std(Ytrain_raw, 0, 1); % 1 x freqs
                    stdFreq(stdFreq==0) = 1;

                case 'absolute'                                         
                    % only subtract the training mean per column 
                    meanFreq = mean(Ytrain_raw, 1);   % 1 x lags        
                    stdFreq  = ones(1, size(Ytrain_raw,2));  % no scaling

                otherwise
                    error('baselinetype must be ''zscore'' or ''absolute''.');
            end

            % normalize training Y
            Ytrain = (Ytrain_raw - repmat(meanFreq, size(Ytrain_raw,1), 1)) ./ repmat(stdFreq, size(Ytrain_raw,1), 1);

            % normalize test Y with training stats
            Ytest_norm = (Ytest_raw - repmat(meanFreq, size(Ytest_raw,1), 1)) ./ repmat(stdFreq, size(Ytest_raw,1), 1);
                

            % Fit PLS on training set
            [XL, yl, XS, YS, beta, PCTVAR, mse, stats] = plsregress(Xtrain, Ytrain, ncomp);
            if wantDetails, b_rep(kf,:,:) = beta; end          

            % Predict on test set (in normalized Y-space) and inverse-transform
            nTest = sum(test_idx);
            Ypred_norm = [ones(nTest,1), Xtest] * beta;

            % inverse zscore using training mean/std
            Ypred = Ypred_norm .* repmat(stdFreq, nTest, 1) + repmat(meanFreq, nTest, 1);

            % store predictions
            Y_pred(test_idx, :) = Ypred;
            Y_pred_norm(test_idx, :) = Ypred_norm;                         
            Y_true_norm(test_idx, :) = Ytest_norm;                         

            % accumulate prediction error and baseline error (training mean)
            SSE_rep = SSE_rep + sum((Ytest_raw - Ypred).^2, 'all');
            SST_rep = SST_rep + sum((Ytest_raw - meanFreq).^2, 'all');

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

        % prediction R^2 for this rep: 1 - NMSE, pooled over parcels and freqs
        r2_vals(rep) = 1 - SSE_rep / SST_rep;

        % median identification rank (as in main script)       
        r = corr(Y_raw', Y_pred');                              
        r_rank = zeros(nObs);                                  
        for k0 = 1:nObs                                  
            r_rank(k0,:) = tiedrank(r(k0,:));                 
        end                                                  
        medrank_vals(rep) = median(diag(r_rank));              

        if wantDetails                                        
            parcel_rank(:,rep)   = diag(r_rank);               
            beta_all(rep,:,:)    = squeeze(mean(b_rep, 1));     
            rmse_parcel(:,rep)   = sqrt(mean((Y_pred - Y_raw).^2, 2));  
            rmse_matrix(rep,:,:) = sqrt((Y_pred - Y_raw).^2);    
        end
    end % reps

    mean_r2   = mean(r2_vals);
    mean_rank = mean(medrank_vals);                             

    if vip_col > 0
        mean_vip = mean(vip_accum(:, 1:vip_col), 2);  % average VIPs across all folds & reps
    else
        mean_vip = zeros(p,1);
    end

    if wantDetails                                              
        details = struct('R2_rep', r2_vals, 'medrank_rep', medrank_vals, ...
            'parcel_rank', parcel_rank, 'beta_all', beta_all, ...
            'rmse_parcel', rmse_parcel, 'rmse_matrix', rmse_matrix, ...
            'Y_pred_last', Y_pred, ...
            'Y_pred_norm_last', Y_pred_norm, 'Y_true_norm_last', Y_true_norm, ...
            'vip_all', vip_accum(:, 1:vip_col));                        
    end
end
