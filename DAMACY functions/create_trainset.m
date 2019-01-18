function [trainset, testset] = create_trainset(S, cv)
%%
% Function creates two logical matrices where each column is an iteration
% and each row is a measurement.
%
% Input:
        % struct S with i being the number of measurements and with
            % the following fields:
                % Data      The preprocessed data with m x n with m being the number of cells
                            % measured in that measurement and n being the number of surface proteins
                            % measured
                % Labels    
                % ID 
                
        % 2x1 vector cv with the cv(1) being the number of fold and cv(2)
        % the number of iterations. If left empty the most optimal
        % crossvalidation is found. 
%
%
% Output: 
        % logical matrix train with i x j size with i being the number of
        % measurements and j being the number of iterations.
        %
        % logical matrix test with i x j size with i being the number of
        % measurements and j being the number of iterations.
%%

Labels = vertcat(S.Labels); % 0 is control, 1 is diseased

%% what type of cross validation should be used
if isempty(cv)
    cv = min([sum(Labels == 0), sum(Labels == 1)]);
    if cv > 50 % use 7 crossvalidation if the classes are well distributed
        cv = min([sum(Labels == 0), sum(Labels == 1), 7]);
        cv(2) = 20;
        x_message = ['Function performs ' num2str(cv) ' fold cross validation with 20 iterations.'];
        disp(x_message)
    else %if the number of subjects is lower than 50, always use leave one out cv
        x_message = 'Function performs leave atleast one of every class out.';
        disp(x_message)
        cv(2) = 1;
    end
else
    x_message = ['Function performs ' num2str(cv(1)) ' fold cross validation with ' num2str(cv(2)) ' iterations.'];
    disp(x_message)
end
%% create random case and control vectors
idx_case = find(Labels == 1);
idx_control = find(Labels == 0);
R_case = randperm(length(idx_case));
R_control = randperm(length(idx_control));
R_case = idx_case(R_case);
R_control = idx_control(R_control);


%% Divide data according to crossvalidation
nm_folds = cv(1); % number of folds

m = ones(nm_folds+1,1);
n = ones(nm_folds+1,1);
for l1 = 1:nm_folds
    m(l1) = round(length(R_case)/nm_folds*l1); %divide the case
    n(l1) = round(length(R_control)/nm_folds*l1); %divide the control
end
m(l1+1) = m(l1) + 1;
n(l1+1) = n(l1) + 1;


% perform crossvalidation
trainset = zeros(length(Labels),cv(2)*nm_folds);
for l2 = 1:cv(2) % number of iterations
    R_case = randperm(length(idx_case));
    R_control = randperm(length(idx_control));
    R_case = idx_case(R_case);
    R_control = idx_control(R_control);
    for l1 = 1:nm_folds % number of folds
        R_test = [R_case(m(l1):m(l1+1)-1); R_control(n(l1):n(l1+1)-1)];
        R_train = setdiff([R_control; R_case], R_test);
        trainset(R_train,(l2-1)*nm_folds+l1) = 1;
    end
end


trainset = logical(trainset);
testset = ~trainset;





