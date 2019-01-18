function [b_acc, n_LV] = Damacy_top_cv_acc_lvonly(X, Y, numFolds)
%%
%%
%Calculation accuracy based on 7 fold crossvalidation
%
% Input:
% X:        The selected unfolded histograms (already selected on
%           variables with a certain standard deviation)
% Labels:   Dummy vector containing 0 for control class and +1 for
%           challenged class.
%
% Output:
% b_acc:            Best accuracy among the 7 fold cross validation
% n_LV:             Best number of latent variables
%
% Written by G.H. Tinnevelt on 26-2-2015 at Raboud University Nijmegen
% Edited a lot! by G.H. Tinnevelt and M. Kokla on 21-5-2015 at Raboud University Nijmegen
%% create class vector Y
%Time: 0 seconds

X = X - repmat(mean(X), size(X,1), 1);
mean_Y = mean(Y);

%% OPLS-DA cross validation

idx1 = find(Y == -1);
idx2 = find(Y == 1);
R1 = randperm(length(idx1));
R2 = randperm(length(idx2));
R1 = idx1(R1);
R2 = idx2(R2);
use_7fold = 0;
if length(Y) > 50 % use 7 crossvalidation if the classes are well distributed
    cv = min([sum(Y == -1), sum(Y == 1), 7]);
    if cv == 7
        use_7fold = 1;
    end
else %if the number of subjects is lower than 50, always use leave one out cv
    cv = min([sum(Y == -1), sum(Y == 1)]);
end
if use_7fold == 1 % 7 fold cross validation
    m = ones(cv+1,1);
    n = ones(cv+1,1);
    acc = zeros(numFolds, cv, 20);
    for l1 = 1:cv
        m(l1) = round(length(R1)/cv*l1);
        n(l1) = round(length(R2)/cv*l1);
    end
    m(l1+1) = m(l1) + 1;
    n(l1+1) = n(l1) + 1;
    
    for l3 = 1:20
        R1 = randperm(length(idx1));
        R2 = randperm(length(idx2));
        R1 = idx1(R1);
        R2 = idx2(R2);
        for l1 = 1:cv
            R_test = [R1(m(l1):m(l1+1)-1); R2(n(l1):n(l1+1)-1)];
            R_train = setdiff([R1; R2], R_test);
            X_test = X(R_test,:);
            Y_test = Y(R_test,:);
            Y_train = Y(R_train,:) - mean_Y;
            X_train = X(R_train,:);
            [W, Q, W_o, P_o] = OPLS_large(X_train,Y_train,numFolds);
            YHAT = OPLSpred_large(X_test,P_o,W_o,W,Q);
            for l2 = 1:numFolds+1
                Y_pred1 = YHAT(:,l2) > 0 - mean_Y;
                Y_pred2 = YHAT(:,l2) < 0 - mean_Y;
                Y_pred = Y_pred1 - Y_pred2;
                acc(l2, l1, l3) = sum(Y_pred == Y_test)/length(Y_test);
            end
        end
    end
    [b_acc, n_LV] = max(mean(acc(:,:), 2));
else % leave part out cross validation (if equal class, leave one of every class out)
    m = ones(cv+1,1);
    n = ones(cv+1,1);
    acc = zeros(numFolds+1, cv);
    for l1 = 1:cv
        m(l1) = round(length(R1)/cv*l1);
        n(l1) = round(length(R2)/cv*l1);
    end
    m(l1+1) = m(l1) + 1;
    n(l1+1) = n(l1) + 1;
    R1 = randperm(length(idx1));
    R2 = randperm(length(idx2));
    R1 = idx1(R1);
    R2 = idx2(R2);
    
    for l1 = 1:cv
        R_test = [R1(m(l1):m(l1+1)-1); R2(n(l1):n(l1+1)-1)];
        R_train = setdiff([R1; R2], R_test);
        X_test = X(R_test,:);
        Y_test = Y(R_test,:);
        Y_train = Y(R_train,:) - mean_Y;
        X_train = X(R_train,:);
        [W, Q, W_o, P_o] = OPLS_large(X_train,Y_train,numFolds);
        YHAT = OPLSpred_large(X_test,P_o,W_o,W,Q);
        for l2 = 1:numFolds+1
            Y_pred1 = YHAT(:,l2) > 0 - mean_Y;
            Y_pred2 = YHAT(:,l2) < 0 - mean_Y;
            Y_pred = Y_pred1 - Y_pred2;
            acc(l2, l1) = sum(Y_pred == Y_test)/length(Y_test);
        end
    end
end
[b_acc, n_LV] = max(mean(acc(:,:), 2));
n_LV = n_LV - 1; 
end
