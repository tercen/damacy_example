function [S_mean] = paired_centering(S, center_mode)
%%
% Function that performs paired centering based on either both measurements
% or only on the control method.
% Input:
% S:            A nx1 struct, where n is the number of individuals and
%               containing the folowwing fields: Data, Labels and ID.
%               !Labels of the control needs to be 0!
% center_mode:  Type of centering that needs to applied:
%               1 mean center over control measurements
%               2 median center over control measurements
%               3 mean center per individual
%               4 median center per individual

% Output:
% S_mean:       A nx1 struct, where n is the number of individuals and
%               containing the folowwing fields: Data, Labels and ID.
%               Where the Data is centered.
%
% Written by G.H. Tinnevelt on 4-3-2015 at Radboud University Nijmegen
%%
S_mean = S;
id = vertcat(S.ID);
Labels = vertcat(S.Labels); % 0 is control, 1 is challenged
ff = unique(id);
for l1 = 1:length(ff)
    idxid = find(id == ff(l1) == 1);
    if center_mode == 1 %mean center over control measurement(s)
        mean_control = mean_of_means(S(logical(idxid(Labels(idxid) == 0))));
        for l2 = 1:length(idxid)
            S_mean(idxid(l2)).Data = S(idxid(l2)).Data - repmat(mean_control,size(S(idxid(l2)).Data,1),1);
        end
    elseif center_mode == 2 %median center over control measurement(s)
        median_control = mean_of_medians(S(logical(idxid(Labels(idxid) == 0))));
        for l2 = 1:length(idxid)
            S_mean(idxid(l2)).Data = S(idxid(l2)).Data - repmat(median_control,size(S(idxid(l2)).Data,1),1);
        end
        
    elseif center_mode == 3 % mean center over all measurements
        mean_all = mean_of_means(S(idxid));
        for l2 = 1:length(idxid)
            S_mean(idxid(l2)).Data = S(idxid(l2)).Data - repmat(mean_all,size(S(idxid(l2)).Data,1),1);
        end
    elseif center_mode == 4 % median center over all measurements
        median_all = mean_of_medians(S(idxid));
        for l2 = 1:length(idxid)
            S_mean(idxid(l2)).Data = S(idxid(l2)).Data - repmat(median_all,size(S(idxid(l2)).Data,1),1);
        end
    end
end
