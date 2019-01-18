function [S_mean,median_all] = unpaired_centering_numcells(S, center_mode)
%%
% Function that performs centering based on either all controls or on all
% measurements.
%
% Input:
% S:            A nx1 struct, where n is the number of individuals and
%               containing the folowwing fields: Data, Labels.
% center_mode:  Type of centering that needs to applied:
%               1 mean center over control measurements
%               2 median center over control measurements
%               3 mean center over all measurements
%               4 median center over all measurements
%               5 mean center for each individuals
%               6 median center for each individuals
%
% Output:
% S_mean:       A nx1 struct, where n is the number of individuals and
%               containing the folowwing fields: Data, Labels and ID.
%               Where the Data is centered.
%
% Written by G.H. Tinnevelt on 23-2-2015 at Radboud University Nijmegen
% (changed introducing the mean/median center for each individuals 28/04/2015)
%%
S_mean = S;
Labels = vertcat(S.Labels); % 0 is control, 1 is diseased
if isfield(S(1), 'train')
    trainset = logical(vertcat(S.train));
else
    trainset = true(length(S),1);
end
if center_mode == 1
    mean_control = mean_of_means(S(logical((Labels == 0).*trainset)));
    for l1 = 1:length(S)
        S_mean(l1).Data = S(l1).Data - repmat(mean_control,size(S(l1).Data,1),1);
    end
elseif center_mode == 2
    median_control = mean_of_medians(S(logical((Labels == 0).*trainset)));
    for l1 = 1:length(S)
        S_mean(l1).Data = S(l1).Data - repmat(median_control,size(S(l1).Data,1),1);
    end
    
elseif center_mode == 3 % mean center over all individuals
    mean_all = mean_of_means(S(trainset));
    for l1 = 1:length(S)
        S_mean(l1).Data = S(l1).Data - repmat(mean_all,size(S(l1).Data,1),1);
    end
elseif center_mode == 4 % median center over all individuals
    median_all = mean_of_medians(S(trainset));
    for l1 = 1:length(S)
        S_mean(l1).Data = S(l1).Data - repmat(median_all,size(S(l1).Data,1),1);
    end
    
elseif center_mode == 5 % mean for each individuals
    
    for l1 = 1:length(S)
        S_mean(l1).Data = S(l1).Data - repmat(mean(S(l1).Data),size(S(l1).Data,1),1);
    end
    
elseif center_mode == 6 % median for each individuals
    
    for l1 = 1:length(S)
        S_mean(l1).Data = S(l1).Data - repmat(median(S(l1).Data),size(S(l1).Data,1),1);
    end
end
