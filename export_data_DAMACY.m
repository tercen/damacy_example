function export_data_DAMACY(output_data_filename, S, VariableNames, S_scaled, S_scores_DAMACY, DAMACY_output)


headername = [VariableNames VariableNames]; %create the VariableNames plus PC names as header

for l2 = 1:size(S_scaled(1).Data,2)
    headername = [headername {['PC' num2str(l2)]}];
end
headername = [headername 'Damacyweights'];
p = length(headername); %number of variables
headername{p+1} = 'Label';
headername{p+2} = 'ID';
x = zeros(DAMACY_output.binsize); 
x(DAMACY_output.indices) = DAMACY_output.weights; 

Data = cell(length(S),1); 
for l1 = 1:length(S)
    BinArray = [];
    tmp = S_scores_DAMACY(l1).Data; 
    for l2 = 1:size(DAMACY_output.edges,1)
        [~, BinArray(:, l2)] = histc(tmp(:, l2),[-inf DAMACY_output.edges(l2, 2:end-1) inf]); %position of cells in histogram
    end
    W = diag(x(BinArray(:,1), BinArray(:,2)));
    Data{l1} = [S(l1).Data S_scaled(l1).Data S_scaled(l1).Data*DAMACY_output.PCAloading' W];
    Data{l1}(:, p+1) = S(l1).Labels;
    Data{l1}(:, p+2) = S(l1).ID;
end
Data = cell2mat(Data); 

csvwrite(output_data_filename, Data);

fid = fopen(output_data_filename, 'w') ;
fprintf(fid, '%s,', headername{1,1:end-1}) ;
fprintf(fid, '%s\n', headername{1,end}) ;
fclose(fid) ;
dlmwrite(output_data_filename, Data, '-append') ;