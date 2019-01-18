function [S, VariableNames] = importdata_simplestruct_big_csv_operator(filename)
%%

%%
delimiter = ',';
Data = dlmread(filename,delimiter, 1,0);
fid = fopen(filename);
headername = strsplit(fgetl(fid),delimiter); %load the first line of the data containing the header

Label_ID = Data(:, end-1:end);
unique_Label_ID = unique(Label_ID, 'rows'); 

for l1 = 1:size(unique_Label_ID,1)
    S(l1).Data = Data(ismember(Label_ID,unique_Label_ID(l1,:), 'rows'),1:end-2);
    S(l1).Labels = unique_Label_ID(l1,1); 
    S(l1).ID = unique_Label_ID(l1,2); 
end

VariableNames = headername(1:end-2);
end