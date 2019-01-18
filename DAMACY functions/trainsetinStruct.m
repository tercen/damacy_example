function [S] = trainsetinStruct(S,trainset, testset)
%%
% Adds a field to the struct S called S.train which contains a 1 if it
% belongs to the test set or 0 if it belongs to the testset
%%

for l1 = 1:length(S)
    if trainset(l1) == 1
        S(l1).train = 1;
    elseif testset(l1) == 1
        S(l1).train = 0;
    else
        S(l1).train = 1;
        warning('train and test set are inconsistent. Sample is set to belong to train set')
    end
end
end