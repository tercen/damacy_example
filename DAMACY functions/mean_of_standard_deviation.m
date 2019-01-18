function [sigma] = mean_of_standard_deviation(S)
%%

%%

n_individuals = length(S);
sigma_square = 0;
for l1 = 1:n_individuals
    sigma_square = sigma_square + sum((S(l1).Data).^2)/(size(S(l1).Data,1)-1);
end
sigma = sqrt(sigma_square/n_individuals); 
end
