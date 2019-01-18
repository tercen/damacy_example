function DAMACY(parameter_filename, data_filename, output_data_filename, output_CVperformance_filename)

P = import_parameterfilename(parameter_filename); 
[S, VariableNames] = importdata_simplestruct_big_csv_operator(data_filename);
[trainset, testset] = create_trainset(S, []);
[results_CV_DAMACY, DAMACY_output, S_scores_DAMACY, Histograms_DAMACY] = MFC_crossvalidate_DAMACY_numcells(S,trainset,testset, P.pre_process_modes, P.pc_used,P.paired_data,P);
[S_scaled] = Pre_process_MFC_takingintoaccountthenumberofcells(S, P.paired_data, P.pre_process_modes);
export_data_DAMACY(output_data_filename, S, VariableNames, S_scaled, S_scores_DAMACY, DAMACY_output);
export_CVperformance(output_CVperformance_filename, results_CV_DAMACY); 