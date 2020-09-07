function [Kxa, ya] = dataset_copd_simul(pars)
% Simulation on COPD by generating multi-dimensional phenotypes

h2 = pars.h2; % heritability
sig = pars.sig; % noise variance from low dimension phenotype to high dimension
ydim = pars.ydim; % phenotype dimension
data_dir = pars.data_dir; % data folder

% set path
precision = sprintf('h2_%%.%dfS2_%%.%dfydim_%%d', digitNum(h2), digitNum(sig));
data_subdir = sprintf(precision, h2, sig, ydim);
data_path = fullfile(data_dir,data_subdir, 'pheno_data.mat');

% load GRM
GRMFile = fullfile(data_dir, 'GRM', 'CG10kNhwHg19Clean_v2_Mar2013_grm.grm');
GRMid = fullfile(data_dir, 'GRM', 'CG10kNhwHg19Clean_v2_Mar2013_grm.grm.id');
[Nsubj, subjId, Kxa] = ParseGRM(GRMFile, GRMid);
    
% load pheotype data or create a mat file for phenotype
if exist(data_path, 'file')
    load(data_path);
else
    files = dir(fullfile(data_dir, data_subdir, '*.txt'));
    ya = zeros(size(Kxa,1), ydim, length(files));
    % load phenotype data
    for i = 1:length(files)
        df = readtable(fullfile(data_dir, data_subdir, files(i).name));
        df = df(:,4:end);
        ya(:,:,i) = table2array(df);
    end
    save(data_path, 'ya', '-v7.3');
end
