function export_CVperformance(output_CVperformance_filename, results_CV_DAMACY)

fid = fopen(output_CVperformance_filename, 'w')
fprintf(fid, '%s,', 'accuracy')
fprintf(fid, '%f', results_CV_DAMACY.facc)
fprintf(fid, '\n%s,', 'specificity')
fprintf(fid, '%f', results_CV_DAMACY.fspec);
fprintf(fid, '\n%s,', 'sensitivity')
fprintf(fid, '%f', results_CV_DAMACY.fsens);
fprintf(fid, '\n%s,', 'ModeLVn')
fprintf(fid, '%d', mode(results_CV_DAMACY.n_LV))
fclose(fid)