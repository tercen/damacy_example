function [results_CV_DAMACY, DAMACY_output, S_scores, Histograms] = MFC_crossvalidate_DAMACY_numcells(S,trainset,testset, pre_process_modes, pc_used,paired_data, DAMACY_settings)
%%
% Function used to crossvalidate the DAMACY model
%
% Input:
% S                 struct with the fields raw data, labels, ID
% trainset          i x n_iterations matrix with i being the number of samples and
%                       n_iterations being the number of iterations, contains a 1 if belongs to
%                       trainset and 0 if belong to testset.
% testset           i x n_iterations matrix (reverse of trainset)
% pre_process_modes 3 x 1 cell contains the different preprocess_modes,
%                   like log transformation, centring, scaling
%
% pc_used           2x1 vector containing which PC are used for example
%                   pc_used = [1 2] for the first two PCs
%
% Output:
% Results struct containing the following fields
% facc          scalar containing the final accuracy
% fsens         scalar containing the final sensitivity
% fspec         scalar containing the final specificity
% yhat          i x 1 vector containing the prediction for y
% w_cross       containing all the crossvalidated weights
% ind_cross     containing all the indices of the weights
%
%
% Written by G.H. Tinnevelt at Radboud University on 16 February 2016
%% create training and test set
if isfield(S(1), 'train')
    S_train = S;
    trainset2 = logical(vertcat(S.train));
    S = S(trainset2); 
end
n_iterations = size(trainset,2); % number of iterations
n_repeats = sum(testset(1,:),2);
n_folds = n_iterations/n_repeats;
acc = zeros(1,n_iterations);
spec = zeros(1,n_iterations);
sens = zeros(1,n_iterations);
n_LV = zeros(n_iterations,1);
yhat = zeros(length(S),n_repeats);
Labels = vertcat(S.Labels); % 0 is control, 1 is diseased
w_cross = cell(n_iterations,1);
for l1 = 1:n_iterations
    [S] = trainsetinStruct(S,trainset(:,l1), testset(:,l1)); %Add train field to struct
    [S_scaled] = Pre_process_MFC_takingintoaccountthenumberofcells(S, paired_data, pre_process_modes);
    [~, ~, DAMACY_output] = DAMACY_main_script(S_scaled, paired_data,DAMACY_settings, pc_used);
    acc(l1) = DAMACY_output.accuracy*sum(testset(:,l1));
    n_LV(l1) = DAMACY_output.OPLS_model.n_LV;
    yhat(testset(:,l1),ceil(l1/n_folds)) = DAMACY_output.yhat;
%     ind_var_DAMACY(testset(:,l1)) = DAMACY_output.ind_DAMACYvariance;
%     SS_res(testset(:,l1)) =  DAMACY_output.OPLS_model.SSresiduals;
    spec(l1) = DAMACY_output.spec*sum(testset(Labels == 0,l1));
    sens(l1) = DAMACY_output.sens*sum(testset(Labels == 1,l1));
%     w_cross(l1) = {DAMACY_output.weights};
end
facc = sum(acc)/length(S);
fspec = sum(spec)/length(S(Labels == 0));
fsens = sum(sens)/length(S(Labels == 1));

results_CV_DAMACY.facc = facc;
results_CV_DAMACY.n_LV = n_LV;
results_CV_DAMACY.fspec = fspec;
results_CV_DAMACY.fsens = fsens;
results_CV_DAMACY.yhat = yhat;
% results_CV_DAMACY.ind_var_DAMACY = ind_var_DAMACY;
% results_CV_DAMACY.SSresiduals = SS_res;
% results_CV_DAMACY.w_cross = w_cross;

%% final model
if exist('S_train', 'var') == 1
    S = S_train; 
else
    S = rmfield(S,'train');
end

[S_scaled] = Pre_process_MFC_takingintoaccountthenumberofcells(S, 0, pre_process_modes);
[S_scores, Histograms, DAMACY_output] = DAMACY_main_wo_cv(S_scaled, 0, DAMACY_settings, pc_used, mode(n_LV));

if  exist('S_train', 'var') == 1
    yhat_instruct = zeros(length(S),1); 
    trainset = logical(vertcat(S.train)); 
    yhat_instruct(trainset) = yhat; 
    yhat_instruct(~trainset) = DAMACY_output.yhat; 
    DAMACY_output.yhat = yhat_instruct; 
else
    DAMACY_output.yhat = yhat;
end
