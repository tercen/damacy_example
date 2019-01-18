function [mu] = mean_of_medians(S)
%%

%%
mu = 0; 
for l1 = 1:length(S)
    mu = mu + median(S(l1).Data);
end
mu = mu/length(S);