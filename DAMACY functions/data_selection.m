function [S_selected, VariableNames, selection_mode] = data_selection(S, VariableNames, selection_mode)
%%
% Function to create a selection of the data.
% Input
% S                 raw data struct
% VariableNames     j x 1 cell with names of j variables
% selection_mode    3 x 1 or 4 x 1 cell containing the modes for selecting 
%                   the data with the first option being the variable 
%                   selection, the second the control label, the third the 
%                   challenged label and optional fourth one containing 
%                   the tube selection
%                   if left empty, questions appear for data selection
% 
% Output:
% S_selected        struct containing the selected data
% VariableNames     j* x 1 cell containing the selected j* variables
% selection_mode    containing the options chosen
%
% Written by G.H. Tinnevelt at Radboud University at 19-6-2015
%%

% select a Tube if any
if isfield(S(1), 'Tube')
    if isempty(selection_mode)
        [S, tubenr] = select_tube(S, []); 
        selection_mode{4} = tubenr;
    else
        [S] = select_tube(S, selection_mode{4}); 
    end
end
Labels = vertcat({S.Labels});
if isempty(selection_mode)
    check_var_used = 0;
    while check_var_used == 0 %used to create the warning
        var_string = 'What variables should be used (give in [] or as 1:5)?';
        for l1 = 1: length(VariableNames)
            var_string = [var_string '\n' num2str(l1) ': ' VariableNames{l1}];
        end
        var_string = [var_string '\n' 'Variables used: '];
        var_used = input(var_string);
        % check if the var_used selected are indeed possible// not too many
        % not too less var_used
        if length(var_used) <= length(VariableNames) && max(var_used) <= length(VariableNames) && min(var_used) > 0
            check_var_used = 1; %stop while loop
        else % if var_used is not correct, give warning and start over
            warning('Wrong input, please enter valid variables.')
        end
    end
    check_label = 0;
    while check_label == 0 || check_label == 1; %used to create the warning
        if check_label == 0 %ask for control label
            select_label(check_label+1) = {input('What is the label of the control in the data: ')};
        elseif check_label == 1 %ask for diseased/challenged label
            select_label(check_label+1) = {input('What is the label of the diseased in the data: ')};
        end
        if isa(Labels{1},'char') || isa(Labels{1}, 'cell')
            if sum(ismember(Labels, select_label{check_label+1})) > 0
                check_label = check_label + 1;
            else
                warning('Wrong input, please enter valid labels.')
            end
        elseif isa(Labels{1},'double') %if double instead of char or cell, first convert to cell using cellmat when comparing
            if sum(ismember(cell2mat(Labels), cell2mat(select_label(check_label+1)))) > 0
                check_label = check_label + 1;
            else
                warning('Wrong input, please enter valid labels.')
            end
        end
    end
    selection_mode{1} = var_used;
    selection_mode(2:3) = select_label;
else
    var_used = selection_mode{1};
    select_label = selection_mode(2:3);
end

%%
S_selected = S;
unselect = []; %create a matrix containing all samples that are not needed
for l1 = 1:length(S)
    S_selected(l1).Data = S(l1).Data(:, var_used); %remove all variables that are not needed
    if isa(Labels{1},'char')
        if strcmp(S(l1).Labels, select_label{1}) %controls
            S_selected(l1).Labels = 0; %give control labels a 0
        elseif strcmp(S(l1).Labels, select_label{2}) %case
            S_selected(l1).Labels = 1;%give challenged labels a 1
        else
            unselect = [unselect l1];
        end
    elseif isa(Labels{1},'double')
        if ismember(S(l1).Labels,  select_label{1}) %controls
            S_selected(l1).Labels = 0;%give control labels a 0
        elseif ismember(S(l1).Labels,  select_label{2}) %case
            S_selected(l1).Labels = 1;%give challenged labels a 1
        else
            unselect = [unselect l1];
        end
    end
end
S_selected(unselect) = []; %remove all samples that are not needed

VariableNames = VariableNames(var_used); %remove all variables that are not needed