function o = import_31P_n_time_from_csv(csv_filepath)
% returns RNA conc and time from CSV exported by integrate_n_plot_31P

%% Params to run "locally"
%===========================
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    csv_filepath = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra/190415_IN115a_pA20_co-P50N525_100uM_298K_600/analysis/datasave/RNA_A20_P50N525.csv';
end

%% Parsing
%===========================
o = csvread(csv_filepath);

end