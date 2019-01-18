function [mu] = mean_of_means(S)
%%

%%
mu = 0; 
for l1 = 1:length(S)
    mu = mu + mean(S(l1).Data);
end
mu = mu/length(S);