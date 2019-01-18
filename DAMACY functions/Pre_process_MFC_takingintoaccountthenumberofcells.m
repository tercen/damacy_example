function [S_scaled, pre_process_modes] = Pre_process_MFC_takingintoaccountthenumberofcells(S, paired_data, pre_process_modes)
%% Pre-Procesisng of FLOW cytometry data
% Function that pre-processes flow cytometry S structs.
% The steps include:
% Log transformation
% Centering
% Scaling

% Note: paired centring/scaling does not take into account number of cells

% Written by G.H. Tinnevelt on 12-6-2015 at Radboud University Nijmegen
%% extract information

Labels = vertcat(S.Labels); % 0 is control, 1 is diseased

%% log transformation
% Note: Logicle scaling is used from the SPADE algorithm.
if isempty(pre_process_modes)
    log_mode = input('What way of center method should be used? \n 1 common log transformation \n 2 arcsinh transformation \n 3 subtract min value + common log transformation \n arcsinh tranformation cofactor 150 \n log_mode: ');
else
    log_mode = pre_process_modes(1);
end
S_scaled = S;
if log_mode == 1
    for l1 = 1:length(S_scaled)
        tmp_data = S(l1).Data;
        tmp_data(tmp_data < 1) = 1;
        S_scaled(l1).Data = log10(tmp_data);
    end
elseif log_mode == 2
    X = vertcat(S.Data);
    cofactor = zeros(size(X,2));
    for l2 = 1:size(X,2)
        cofactor(l2) = abs(prctile(X(X(:,l2) < 0,l2),1));
    end
    cofactor(isnan(cofactor)) = 1;
    for l1 = 1:length(S_scaled)
        S_scaled(l1).Data = flow_arcsinh(S(l1).Data',cofactor)';
    end
elseif log_mode == 3
    min_value = prctile(vertcat(S.Data),0.1);
    for l1 = 1:length(S_scaled)
        tmp_data = S(l1).Data;
        tmp_data = tmp_data - repmat(min_value-1,size(tmp_data,1),1);
        tmp_data(tmp_data < 1) = 1;
        S_scaled(l1).Data = log10(tmp_data);
    end
elseif log_mode == 4
    for l1 = 1:length(S_scaled)
        S_scaled(l1).Data = flow_arcsinh(S(l1).Data',150)';
    end
elseif log_mode == 5
    X = vertcat(S.Data);
    Xneg = X(any(X < 0,2),:);
    min_value = median(Xneg) - 2*1.4826*mad(Xneg);
    for l1 = 1:length(S_scaled)
        tmp_data = S(l1).Data;
        tmp_data(any(tmp_data < min_value,2),:) = [];
        S_scaled(l1).Data = flow_arcsinh(tmp_data',150)';
    end
else
    warning('Data is not log transformed!')
end

%% centring

if paired_data == 1
    if isempty(pre_process_modes)
        center_mode = input('What way of center method should be used? \n 1 mean center over control measurements \n 2 median center over control measurements \n 3 mean center per individual \n 4 median center per individual \n 5 mean center over all indviduals \n 6 median center over all individuals \n center_mode: ');
    else
        center_mode = pre_process_modes(2);
    end
    if center_mode == 5 || center_mode == 6
        S_scaled = unpaired_centering_numcells(S_scaled,center_mode -2);
    elseif sum(center_mode == 1:4) == 1
        S_scaled = paired_centering(S_scaled, center_mode);
    else
        S_scaled = S_scaled;
        warning('Data is not centered!')
        
    end
else
    if isempty(pre_process_modes)
        center_mode = input('What way of center method should be used? \n 1 mean center over all controls \n 2 median center over all controls \n 3 mean center over all individuals \n 4 median center over all individuals \n 5 mean center for each individuals \n 6 median center for each individuals \n center_mode: ');
    else
        center_mode = pre_process_modes(2);
    end
    
    if sum(center_mode == 1:6) == 1
        S_scaled = unpaired_centering_numcells(S_scaled, center_mode);
    else
        warning('Data is not centered!')
    end
end

%% scaling
if isfield(S(1), 'train')
    trainset = logical(vertcat(S.train));
else
    trainset = true(length(S),1);
end

if isempty(pre_process_modes)
    auto_mode = input('What way of autoscaling method should be used? \n 1 scale over whole set \n 2 scale over control \n3 scaling for each individual \n 4 no scaling \n auto_mode: ');
else
    auto_mode = pre_process_modes(3);
end
if auto_mode == 1
    std_whole = mean_of_standard_deviation(S_scaled(trainset));
    for l1 = 1:length(S_scaled)
        S_scaled(l1).Data = S_scaled(l1).Data./(repmat(std_whole, size(S_scaled(l1).Data,1),1));
    end
elseif auto_mode == 2
    std_control = mean_of_standard_deviation(S_scaled(logical((Labels == 0).*trainset)));
    for l1 = 1:length(S_scaled)
        S_scaled(l1).Data = S_scaled(l1).Data./(repmat(std_control, size(S_scaled(l1).Data,1),1));
    end
elseif auto_mode == 3
    for l1 = 1:length(S_scaled)
        std_indiv = std(S_scaled(l1).Data);
        S_scaled(l1).Data = S_scaled(l1).Data./(repmat(std_indiv ,size(S_scaled(l1).Data,1),1));
    end
else
    warning('Data is not scaled!')
end

%%
pre_process_modes = [log_mode; center_mode; auto_mode];