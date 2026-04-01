%% This script runs partial least squares regression to predict parcel-wise 
% average regional autocorrelation (MEG) from neurochemical/structural information 

% CS 2024/25

%% Load neurochemical/structual maps
load('allmaps_concat2_selected.mat')

% Don't zscore globally, but do so in CV loop - keep raw maps
X_raw = double(allmaps);   % original predictors (200 x nMaps)

%% Get MEG data
% get autocorrelation data
% adjusted for age, sex, icv
load('adjusted_MEGACall_sorted.mat')

% get average values (parcel x freq)
gACav = squeeze(mean(all_adj, 1));
Y_raw = gACav;

% Choose normalization method (across parcels, will be applied within-fold using training stats)
baselinetype = 'absolute'; 

%% Now do a quick check on predictive performance (on full data for exploratory PCTVAR)
% do z-scoring beforehand
allmapsz=zscore(double(allmaps),0,1);

% Use unscaled X_raw and Y_raw for this exploratory PLS (plsregress will center internally)
[xl,yl,XS,YS,beta,PCTVAR,mse,stat]=plsregress(allmapsz,Y_raw,30); % assume 30 components
figure
plot(cumsum(100*PCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in y'); % exploratory

%% Predict MEG power, compute ranks using ncomp
kfold=10;
ncomp=10; %number of components for partial least squares (PLS)

% We'll use X_raw and Y_raw and perform per-fold standardization inside the loop
Y1=zeros(size(Y_raw));

for rep=1:50

    cv = cvpartition(200,"Kfold",kfold);
    Y1=zeros(size(Y_raw));
    vipScore = nan(size(X_raw,2),kfold);
    b = nan(kfold, size(X_raw,2)+1, size(Y_raw,2)); % store betas per fold

    for ifold=1:kfold
        tr = cv.training(ifold);
        te = cv.test(ifold);

        %% STANDARDIZE X using training data stats
        muX = mean(X_raw(tr,:),1);
        sigmaX = std(X_raw(tr,:),0,1);
        sigmaX(sigmaX==0) = 1;  % avoid division by zero
        Xtrain = (X_raw(tr,:) - muX) ./ sigmaX;
        Xtest  = (X_raw(te,:)  - muX) ./ sigmaX;

        %% NORMALIZE Y using training data stats (according to baselinetype)
        Ytrain_raw = Y_raw(tr,:);   % training parcels x lag
        Ytest_raw  = Y_raw(te,:);   % test parcels x lag

        % Compute the mean across parcels (rows) for each lag (columns)
        meanVals = repmat(mean(Ytrain_raw, 1), size(Ytrain_raw,1), 1);
        
        switch baselinetype
            case 'absolute'
                Ytrain = Ytrain_raw - meanVals;
                Ytest  = Ytest_raw  - repmat(mean(Ytrain_raw,1), size(Ytest_raw,1), 1);
            otherwise
                error('Unsupported method for baseline normalization: %s', baselinetype);
        end

        %% Fit PLS on training set (note: Ytrain used here)
        [XL,yl,XS,YS,beta,PCTVAR,mse,stats] = plsregress(Xtrain, Ytrain, ncomp);
        b(ifold,:,:) = beta; % store betas for this fold

        % predict on test set; predictions are in normalized Y-space (if Y was normalized)
        Ypred_norm = [ones(sum(te),1), Xtest] * beta;

        % If Y was normalized, bring predictions back into original Y units before storing
        switch baselinetype
%             case 'db'
%                 % Ypred_norm already in dB units because we created Ytrain that way.
%                 Y1(te,:) = Ypred_norm;
            otherwise
                % For absolute/relative/relchange/zscore, convert back to original Y_raw units where needed
                if strcmp(baselinetype,'zscore')
                    % inverse zscore: Y = Ypred_norm .* stdFreq_train + meanFreq_train
                    meanFreq_train = mean(Ytrain_raw,1);
                    stdFreq_train  = std(Ytrain_raw,0,1);
                    stdFreq_train(stdFreq_train==0)=1;
                    Y1(te,:) = Ypred_norm .* repmat(stdFreq_train, sum(te),1) + repmat(meanFreq_train, sum(te),1);
                elseif strcmp(baselinetype,'absolute')
                    meanFreq_train = mean(Ytrain_raw,1);
                    Y1(te,:) = Ypred_norm + repmat(meanFreq_train, sum(te),1);
                elseif strcmp(baselinetype,'relative')
                    meanFreq_train = mean(Ytrain_raw,1);
                    Y1(te,:) = (Ypred_norm + 1) .* repmat(meanFreq_train, sum(te),1);
                elseif strcmp(baselinetype,'relchange')
                    meanFreq_train = mean(Ytrain_raw,1);
                    Y1(te,:) = (Ypred_norm .* repmat(mean(Ytrain_raw,1), size(Ypred_norm,1),1)) + repmat(mean(Ytrain_raw,1), size(Ypred_norm,1),1);
                else
                    % fallback: treat as no inverse transform (should not happen)
                    Y1(te,:) = Ypred_norm;
                end
        end

        % normalized PLS weights for each fold
        W0 = stats.W ./ sqrt(sum(stats.W.^2,1));

        % VIP scores for ncomp and fold
        p = size(XL,1);
        sumSq = sum(XS.^2,1).*sum(yl.^2,1);
        vipScore(:, ifold) = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));       
    end % ifold

    beta_all(rep,:,:) =  squeeze(mean(b, 1)); % intercept is there!
    all_vip(:,rep) = mean(vipScore,2);

    % compute metrics across folds
    mean_obs = mean(Y_raw);  % use raw Y (original units) for R2 denominator
    r2 = 1-(mean((Y_raw-Y1).^2)/(mean((Y_raw-mean_obs).^2)));
    R2_all(rep) = r2;

    % Compute RMS error per parcel for this repetition
    rmse_parcel(:,rep) = sqrt(mean((Y1 - Y_raw).^2, 2))';  % size: [200 x 1]
    % RMS for each frequency bin and parcel
    rmse_matrix(rep,:,:) = sqrt((Y1 - Y_raw).^2);

    r=corr(Y_raw',Y1'); % compute r

    % median rank
    for k0=1:200
        r_rank(k0,:)=tiedrank(r(k0,:));
    end
    parcel_rank(:,rep) = diag(r_rank);
    medrank(rep)=(median(diag(r_rank)));

end % rep

average_R2 = mean(R2_all);
average_rank = mean(medrank);
average_parcelrank = mean(parcel_rank,2);
average_vip = mean(all_vip,2);
average_weights = squeeze(mean(beta_all(:,2:end,:),1));


% load labels for the maps
[VIPs,I] = sort(average_vip, 'descend');
selected = all_labels(I,:);
maps_selected = selected;

% save important maps
maps_selected.scores = VIPs;
maps_selected.predvariable = repmat('Pow_10ncomp',length(VIPs),1);

writetable(maps_selected, 'maps_allselected_ACall_abs_10ncomp_noleak.xlsx')

% save prediction results (save Y_raw as Y to keep naming similar)
Y = Y_raw;
save('ACall_allselected_abs_10ncomp_pls_noleak.mat', 'average_R2', 'average_rank', ...
    'average_vip', 'average_weights', 'beta_all', 'maps_selected', 'parcel_rank', 'rmse_matrix', 'rmse_parcel')

save('orig_pred_ACall_allselected_abs_10ncomp_noleak.mat', 'Y', 'Y1')

%% Now create 0 distribution using spin-tests (updated to avoid leakage)
% do spin test to get a null distribution
[p_spin, r_dist, permid] = spin_test_old(rand(1,200),rand(1,200), 'parcellation_name','schaefer_200');

kfold=10;
cv = cvpartition(200,"Kfold",kfold);
Y1=zeros(size(Y_raw));

for k=1:1000
  % Permute rows of X_raw according to spin permutation
  Xperm = X_raw(permid(:,k),:);

  % Reset Y1 for this spin
  Y1=zeros(size(Y_raw));

  for ifold=1:kfold
      tr = cv.training(ifold);
      te = cv.test(ifold);

      % standardize permuted X using training data from permuted X
      muX = mean(Xperm(tr,:),1);
      sigmaX = std(Xperm(tr,:),0,1);
      sigmaX(sigmaX==0)=1;
      Xtrain = (Xperm(tr,:) - muX) ./ sigmaX;
      Xtest  = (Xperm(te,:)  - muX) ./ sigmaX;

      % standardize Y using training data (same as before)
      Ytrain_raw = Y_raw(tr,:);
      Ytest_raw  = Y_raw(te,:);
      switch baselinetype
          case 'absolute'
              meanVals = repmat(mean(Ytrain_raw,1), size(Ytrain_raw,1),1);
              Ytrain = Ytrain_raw - meanVals;
              Ytest  = Ytest_raw  - repmat(mean(Ytrain_raw,1), size(Ytest_raw,1),1);
          case 'zscore'
              meanFreq = repmat(mean(Ytrain_raw,1), size(Ytrain_raw,1), 1);
              stdFreq  = repmat(std(Ytrain_raw,0,1), size(Ytrain_raw,1), 1);
              stdFreq(stdFreq==0)=1;
              Ytrain = (Ytrain_raw - meanFreq) ./ stdFreq;
              meanFreq_test = repmat(mean(Ytrain_raw,1), size(Ytest_raw,1), 1);
              stdFreq_test  = repmat(std(Ytrain_raw,0,1), size(Ytest_raw,1), 1);
              stdFreq_test(stdFreq_test==0)=1;
              Ytest = (Ytest_raw - meanFreq_test) ./ stdFreq_test;
          otherwise
              % for spins, implement same handling as above if needed (here we keep zscore typical)
              meanFreq = repmat(mean(Ytrain_raw,1), size(Ytrain_raw,1), 1);
              stdFreq  = repmat(std(Ytrain_raw,0,1), size(Ytrain_raw,1), 1);
              stdFreq(stdFreq==0)=1;
              Ytrain = (Ytrain_raw - meanFreq) ./ stdFreq;
              meanFreq_test = repmat(mean(Ytrain_raw,1), size(Ytest_raw,1), 1);
              stdFreq_test  = repmat(std(Ytrain_raw,0,1), size(Ytest_raw,1), 1);
              stdFreq_test(stdFreq_test==0)=1;
              Ytest = (Ytest_raw - meanFreq_test) ./ stdFreq_test;
      end

      % build PLS on training set (perm)
      [xl,yl,XS,YS,beta,PCTVAR,mse,stat]=plsregress(Xtrain,Ytrain,ncomp);

      % predict power in test parcels (back-transform to original Y units if needed)
      Ypred_norm = [ones(sum(te),1), Xtest] * beta;
      % inverse transform
      meanFreq_train = mean(Ytrain_raw,1);
      stdFreq_train  = std(Ytrain_raw,0,1);
      stdFreq_train(stdFreq_train==0)=1;
%     Y1(te,:) = Ypred_norm .* repmat(stdFreq_train, sum(te),1) + repmat(meanFreq_train, sum(te),1);
      Y1(te,:) = Ypred_norm + repmat(mean(Ytrain_raw,1), sum(te),1);

  end

  mean_obs = mean(Y_raw);
  rsqsh(k) = 1-(mean((Y_raw-Y1).^2)/(mean((Y_raw-mean_obs).^2)));

  r=corr(Y_raw',Y1');
  for k0=1:200
     r_rank(k0,:)=tiedrank(r(k0,:));
  end
  spin_parcel_rank(:,k) = diag(r_rank);
  spin_medrank(k)=(median(diag(r_rank)));
end

results_single_spin.rsqsh_spin = rsqsh;
clear r rsqsh

spin_99prctile = prctile(results_single_spin.rsqsh_spin, 99);
spin_99prctile_rank = prctile(spin_medrank, 99);
pval_single=length(find(abs(results_single_spin.rsqsh_spin)>abs(average_R2)))/1000;
pval_medrank=length(find(abs(spin_medrank)>abs(average_rank)))/1000;

% add significance test using spins
save('ACall_allselected_abs_10ncomp_pls_noleak.mat', 'average_R2', 'average_rank', ...
    'average_vip', 'average_weights', 'beta_all', 'maps_selected', 'parcel_rank', 'rmse_matrix', 'rmse_parcel',...
    'results_single_spin', 'pval_single', 'pval_medrank', 'spin_99prctile_rank')

