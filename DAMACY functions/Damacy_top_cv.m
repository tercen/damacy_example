function [b_acc, n_LV,w_cross,t_cross,T_o_cross,W_o_cross,y_hat_cross] = Damacy_top_cv(X, Y, numFolds)
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
    w_cross = zeros(size(X,2),cv, numFolds,20);
    t_cross = cell(cv, numFolds,20);
    T_o_cross = cell(cv, numFolds,20);
    W_o_cross = cell(size(X,2), numFolds,20);
    y_hat_cross = zeros(size(X,1), numFolds, 20);
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
        for l2 = 1:numFolds
            for l1 = 1:cv
                R_test = [R1(m(l1):m(l1+1)-1); R2(n(l1):n(l1+1)-1)];
                R_train = setdiff([R1; R2], R_test);
                X_test = X(R_test,:);
                Y_test = Y(R_test,:);
                Y_train = Y(R_train,:) - mean_Y;
                X_train = X(R_train,:);
                [w,t,~,    q,T_o,P_o,W_o] = OPLS(X_train, Y_train,l2);
                [~,~,yhat] = OPLSpred(X_test,P_o,W_o,w,q,l2);
                Y_pred1 = yhat > 0 - mean_Y;
                Y_pred2 = yhat < 0 - mean_Y;
                Y_pred = Y_pred1 - Y_pred2;
                acc(l2,l1, l3) = sum(Y_pred == Y_test)/length(Y_test);
                w_cross(:,l1,l2, l3) = w;
                t_cross{l1,l2, l3} = t;
                T_o_cross{l1,l2, l3} = T_o;
                W_o_cross{l1,l2, l3} = W_o;
                y_hat_cross(R_test,l2,l3) = yhat;
            end
        end
    end
    [b_acc, n_LV] = max(mean(acc(:,:), 2));
else % leave part out cross validation (if equal class, leave one of every class out)
    m = ones(cv+1,1);
    n = ones(cv+1,1);
    acc = zeros(numFolds, cv);
    w_cross = zeros(size(X,2),cv, numFolds);
    t_cross = cell(cv, numFolds);
    T_o_cross = cell(cv, numFolds);
    W_o_cross = cell(size(X,2), numFolds);
    y_hat_cross = zeros(size(X,1), numFolds);
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
    for l2 = 1:numFolds
        for l1 = 1:cv
            R_test = [R1(m(l1):m(l1+1)-1); R2(n(l1):n(l1+1)-1)];
            R_train = setdiff([R1; R2], R_test);
            X_test = X(R_test,:);
            Y_test = Y(R_test,:);
            Y_train = Y(R_train,:) - mean_Y;
            X_train = X(R_train,:);
            [w,t,~,    q,~,P_o,W_o] = OPLS(X_train, Y_train,l2);
            [~,Tnew_o,yhat] = OPLSpred(X_test,P_o,W_o,w,q,l2);
            Y_pred1 = yhat > 0 - mean_Y;
            Y_pred2 = yhat < 0 - mean_Y;
            Y_pred = Y_pred1 - Y_pred2;
            acc(l2, l1) = sum(Y_pred == Y_test)/length(Y_test);
            w_cross(:,l1,l2) = w;
            t_cross{l1,l2} = t;
            for l3 = 1:length(R_test)
                T_o_cross{R_test(l3),l2} = Tnew_o(l3,:);
            end
            W_o_cross{l1,l2} = W_o;
            y_hat_cross(R_test,l2) = yhat;
        end
    end
    [b_acc, n_LV] = max(mean(acc(:,:), 2));
end
end