function [b_acc,n_LV,w_cross,y_hat_cross] = DAMACY_top_leave_pair_out_cv(X,ID,Y,numFolds)
%%
% Function used to perform leave pair out crossvalidation to decide the
% number of latent variables based on the maximum latent variables allowed.
%
% Input:
% X         Unfolded Histograms (m samples by n bins)
% ID        m by 1 double containing the patients ID (if you do not have
% double IDs the function works as leave one out)
% Y         m by 1 double containing the dummy label (-1 or 1)
% numfolds  the maximum number of latent variables used
%
% Output:
% b_acc     Maximum accuracy of crossvalidation
% n_LV      Optimal number of latent variables
% w_cross   All the weights of crossvalidation
% t_cross   All the scores of crossvalidation
% T_o_cross All the orthogonal scores of crossvalidation
% W_o_cross All the orthogonal weights of crossvalidation
% y_hat_cross All the yhat 
%
% last edited by Marietta Kokla, 20-05-2015
% changes applied : add more inputs(binsize) and more
%  outputs(,w_cross,t_cross,T_o_cross,W_o_cross,y_hat_cross)
% 
%% 
ff = unique(ID);
w_cross = 0;
y_hat_cross = 0; 
% w_cross = zeros(size(X,2),length(ff), numFolds);
% y_hat_cross = zeros(length(Y),length(ff), numFolds); 
acc = zeros(numFolds, length(ff));
for l1 = 1:length(ff)
    if isa(ID, 'cell') == 1
        tmp_id = strcmp(ID, ff(l1));
    elseif isa(ID, 'double') == 1
        tmp_id = ID == ff(l1);
    end
    Xtest = X(tmp_id,:);
    Xtrain = X(~tmp_id,:);
    Ytest = Y(tmp_id,:);
    Ytrain = Y(~tmp_id,:);
    [MC_Ytrain, mean_Y] = mncn(Ytrain);
    [MC_train, mean_X] = mncn(Xtrain);
    MC_test = Xtest - repmat(mean_X, size(Xtest,1),1);
    
   
    [W, Q, W_o, P_o] = OPLS_large(MC_train,MC_Ytrain,numFolds);
    YHAT = OPLSpred_large(MC_test,P_o,W_o,W,Q);
    for i_fold = 1:numFolds+1
            Y_pred1 = YHAT(:,i_fold) > 0 - mean_Y;
            Y_pred2 = YHAT(:,i_fold) < 0 - mean_Y;
            Y_pred = Y_pred1 - Y_pred2;
            acc(i_fold, l1) = sum(Y_pred == Ytest)/length(Ytest);
    end
end
[b_acc, n_LV] = max(mean(acc(:,:), 2));
end






