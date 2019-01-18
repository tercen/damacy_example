function [S_scores, Histograms, DAMACY_output] = DAMACY_main_script(S_scaled, paired_data, DAMACY_settings, pc_used)
% This functions performs the DAMACY algortihm on the already pre-processed
% data.
%
% Input:
% struct S with i being the number of measurements and with
% the following fields:
% Data      The preprocessed data with m x n with m being the number of cells
% measured in that measurement and n being the number of surface proteins
% measured
% Labels
% ID
% paired_data true (1) or false (2)
% parameters (struct with all the parameters
% pc_used
%
% Output
% Struct S with i being the number of measurements with the
% following fields:
% Data      The scores of the PCA model with m x o with m being
% the number of cells measured in that measurement and 0 being
% the number of PCs chosen.
% Labels
% ID
% Histograms with i being the number of measurements by binsize
% default ix500x500 matrix
% struct with the DAMACY outputs
%
% Written by G.H. Tinnevelt on 26-2-2015 at Radboud University Nijmegen
% Last Edited on 04-03-2015 byS. van Staveren at UMC Utrecht in
% collaboration with G.H. Tinnevelt
% Last edited on 2-4-2015 G.H. Tinnevelt able to use other PCs
% Last edited on 12-6-2015 G.H. Tinnevelt, Pre-Processing and FLOOD
% removed, added new plot functions and gating.
% Last edited on 19-6-2015 G.H. Tinnevelt, transformed into function!
%% extract information and set parameters

Labels = vertcat(S_scaled.Labels); % 0 is control, 1 is diseased
if isnumeric(S_scaled(1).ID)
    ID = vertcat(S_scaled.ID);
else
    ID = {S_scaled.ID};
end
    if isfield(S_scaled(1), 'train')
        trainset = logical(vertcat(S_scaled.train));
    else
        trainset = true(length(S_scaled),1);
    end
    
    if isempty(paired_data)
        paired_data = 0;
    end
    
    if isempty(DAMACY_settings)
        minmaxcutoff = 0.01; %part not taken for determining the edges, standard: 0.01  = 1%
        edgesfactor = 1; %Increase the edgesfactor.
        smooth_factor = 5; %smooth_factor  = lambda (article), 1-5 normal, 0 for no smoothing.
        number_bins = 250000;
    else
        minmaxcutoff = DAMACY_settings.minmaxcutoff; %part not taken for determining the edges, standard: 0.01  = 1%
        edgesfactor = DAMACY_settings.edgesfactor; %Increase the edgesfactor.
        smooth_factor = DAMACY_settings.smooth_factor; %smooth_factor  = lambda (article), 1-5 normal, 0 for no smoothing.
        number_bins = DAMACY_settings.number_bins;
    end
    
    %% Create PCA model on whole data with number of PCs given
    % Time: 3 seconds
    
    if isempty(pc_used)
        pc_used = input('Which PCs should be used for the base model? '); % input in brackets [] bijv. [1 3] voor PC 1 en 3.
    end
    N_PCs = length(pc_used);
    [X] = blockscale_struct(S_scaled);
    [~, LDS, VAR] = pcafunction(X);
    
    % for L1=1:(min(length(VAR),7))
    %     if sum(L1 == pc_used) == 1
    %         disp([ '  Used | '  ' | Expained Variance - PC ' num2str(L1) ' = ' num2str(VAR(L1)) ' %']);
    %     else
    %         disp([ 'Unused | '  ' | Expained Variance - PC ' num2str(L1) ' = ' num2str(VAR(L1)) ' %']);
    %     end
    % end
    
    %% Project data
    % Time 0.3
    
    S_scores = S_scaled;
    for l1 = 1:length(S_scaled)
        S_scores(l1).Data = S_scaled(l1).Data*LDS(:,pc_used);
    end
    
    %% calculated individual variance explained by PCA
%     ind_VAR = individual_VAR(S_scaled, S_scores, LDS, pc_used);
    
    %% Create binedges
    X = vertcat(S_scores(trainset).Data);
    if N_PCs == 1
        binsize = 1000;
    else
        binsize = repmat(round(number_bins^(1/N_PCs)),1,N_PCs);
    end
    
    edges = binedges(X, binsize, minmaxcutoff, edgesfactor);
    
    %% Create Histograms
    % Time: 90 seconds
    
    Histograms = zeros([length(S_scores) binsize]);
    for l1 = 1:length(S_scores)
        Histograms(l1,:,:,:,:,:,:) = NDhist(S_scores(l1).Data, edges, smooth_factor);
    end
    
    %% remove empty bins (DAMACY)
    %Note: Empty bins should not be taken into account to circumvent seperation
    %based on empty vs not empty. The empty bins are not 0, but have a very low
    %value due to the smoothing.
    X = reshape(Histograms, size(Histograms,1), prod(binsize));
    ind = std(X) > 10^-6;
    X = X(:,ind);
    
    %% create class vector Y and train and test set (DAMACY)
    %Time: 0 seconds
    
    Y2 = Labels == 0;
    Y1 = ~Y2;
    Y = Y1 - Y2;
    
    %% OPLS (DAMACY)
    
    if paired_data == 1;
        [b_acc,n_LV] = DAMACY_top_leave_pair_out_cv(X,ID,Y,15);% determine the latent variables
    else
%         [b_acc, n_LV,~,y_hat_cross] = Damacy_top_cv(X(trainset,:), Y(trainset), 20);
    [b_acc, n_LV] = Damacy_top_cv_acc_lvonly(X(trainset,:), Y(trainset), 15);

    end
    DAMACY_output.accuracy = b_acc;
    
    if isfield(S_scaled(1), 'train')
        mean_Y = mean(Y(trainset));
        mean_X = mean(X(trainset,:));
        Y_train = (Y(trainset) - repmat(mean_Y,size(Y(trainset),1),1));%./(repmat(std_Y, size(Ytrain,1),1));
        X_train = (X(trainset,:) - repmat(mean_X, size(X(trainset,:),1),1));%./(repmat(std_X, size(Xtrain,1),1));
        X_test = (X(~trainset,:) - repmat(mean_X, size(X(~trainset,:),1),1));%./(repmat(std_X, size(Xtest,1),1));
        
        [w,~,p,    q,~,P_o,W_o] = OPLS(X_train,Y_train,n_LV); %2 class
        
        [~,~,yhat] = OPLSpred(X_test,P_o,W_o,w,q,n_LV); %2 class
        
%         SSres_top = sum((X_test - (X_test*w)*p').^2,2);
%         ind_VAR_DAMACY = sum(ind_VAR(~trainset,:),2).*(1- SSres_top./sum(X_test.^2,2));
        Y_pred1 = yhat >  - mean_Y;
        Y_pred2 = yhat <  - mean_Y;
        Y_pred = Y_pred1 - Y_pred2;
        Y_test = Y(~trainset);
        accf = sum(Y_pred == Y_test)/length(Y_test);
        DAMACY_output.weights = w;
        DAMACY_output.spec = sum(Y_pred((Y_test == -1)) == Y_test(Y_test == -1))/length(Y_test(Y_test == -1));
        DAMACY_output.sens = sum(Y_pred((Y_test == 1)) == Y_test(Y_test == 1))/length(Y_test(Y_test == 1));
        DAMACY_output.accuracy = accf;
        DAMACY_output.base_accuracy = b_acc;
    else
%         w = mean(squeeze(w_cross(:,:,n_LV)),2);
%         yhat = y_hat_cross(:,n_LV);
       
    end
    
    %% DAMACY output
    DAMACY_output.indices = ind;
    DAMACY_output.yhat = yhat;
    DAMACY_output.edges = edges;
    DAMACY_output.pc_used = pc_used;
    DAMACY_output.PCAloading = LDS;
    DAMACY_output.PCAvariance = VAR;
%     DAMACY_output.ind_PCAvariance = ind_VAR;
    DAMACY_output.OPLS_model.n_LV = n_LV;
%     DAMACY_output.ind_DAMACYvariance = ind_VAR_DAMACY;
%     DAMACY_output.OPLS_model.SSresiduals = SSres_top;
    
    
