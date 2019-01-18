function [S_scores, Histograms, DAMACY_output] = DAMACY_base(S_scaled,paired_data, DAMACY_settings, pc_used) 

%%
% Function creates the DAMACY base model
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
%
% Written by G.H. Tinnevelt on 10-2-2015 at Radboud University Nijmegen
%% extract information and set parameters
ID = vertcat(S_scaled.ID);
if isfield(S_scaled(1), 'train')
    trainset = logical(vertcat(S_scaled.train));
else
    trainset = true(length(S_scaled),1);
end

if isempty(paired_data)
    paired_data = 0;
end

if isempty(DAMACY_settings)
    minmaxcutoff = 0.001; %part not taken for determining the edges, standard: 0.01  = 1%
    edgesfactor = 1; %Increase the edgesfactor.
    smooth_factor = 5; %smooth_factor  = lambda (article), 1-5 normal, 0 for no smoothing.
    number_bins = 250000;
    
    DAMACY_settings.minmaxcutoff = minmaxcutoff; %part not taken for determining the edges, standard: 0.01  = 1%
    DAMACY_settings.edgesfactor = edgesfactor; %Increase the edgesfactor.
    DAMACY_settings.smooth_factor = smooth_factor; %smooth_factor  = lambda (article), 1-5 normal, 0 for no smoothing.
    DAMACY_settings.number_bins = number_bins;
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
[~, LDS, VARp] = pcafunction(X);

%% Project data
% Time 0.3

S_scores = S_scaled;
for l1 = 1:length(S_scaled)
    S_scores(l1).Data = S_scaled(l1).Data*LDS(:,pc_used);
end

%% calculated individual variance explained by PCA
ind_VAR = individual_VAR(S_scaled, S_scores, LDS, pc_used);

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

%% DAMACY output
DAMACY_output.edges = edges;
DAMACY_output.pc_used = pc_used;
DAMACY_output.PCAloading = LDS;
DAMACY_output.PCAvariance = VARp;
DAMACY_output.ind_PCAvariance = ind_VAR;
DAMACY_output.DAMACY_settings = DAMACY_settings; 
