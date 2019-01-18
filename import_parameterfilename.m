function P = import_parameterfilename(parameter_filename)

delimiter = ',';
fid = fopen(parameter_filename);

tmp = fgetl(fid); %load the first line of the parameters
while ischar(tmp)
    tmp = strsplit(tmp,delimiter);
    P.(tmp{1}) = cell2mat(cellfun(@str2num,(tmp(2:end)), 'un', 0));
    tmp = fgetl(fid); %load the first line of the parameters
end
fclose(fid);
