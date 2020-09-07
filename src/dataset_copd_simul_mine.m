function [Kxa, ya] = dataset_copd_simul_mine(pars)
% Simulation on COPD by generating multi-dimensional phenotypes

h2 = pars.h2; % heritability
sigma = pars.sigma; % noise variance from low dimension phenotype to high dimension
sigma_ny = pars.sigma_ny; % noise variance from low dimension phenotype to high dimension
yzdim = pars.yzdim; % hidden phenotype dimension
ydim = pars.ydim; % phenotype dimension
data_dir = pars.data_dir; % data folder
N = pars.N; % sample size
T = pars.T; % number of replications

% set path
data_subdir = sprintf('h2_%.2f_sig%.2f_signy%.2f_dzy%d_dy%d_N%d.mat', h2,sigma,sigma_ny,yzdim,ydim,N);
data_path = fullfile(data_dir,data_subdir);
% data_path = fullfile(data_dir,'pheno_data.mat');

% load GRM
GRMFile = fullfile(data_dir, 'GRM', 'CG10kNhwHg19Clean_v2_Mar2013_grm.grm');
GRMid = fullfile(data_dir, 'GRM', 'CG10kNhwHg19Clean_v2_Mar2013_grm.grm.id');
[Nsubj, subjId, Kxa] = ParseGRM(GRMFile, GRMid);
    
% load pheotype data or create a mat file for phenotype
if exist(data_path, 'file')
    load(data_path);
else
    for t=1:T
        % generate Y
        y_h = zeros(N,yzdim);
        for i = 1:yzdim
            y_h(:,i) = mvnrnd(zeros(1,N), Kxa*sigma^2*h2 + sigma^2*(1-h2)*eye(N));
        end
        B = randn(ydim,yzdim);
        [V,D1] = eig(B*B');
        V = V(:,1:yzdim);
        ya(:,:,t) = y_h * V' + sigma_ny*randn(N,ydim);
    end
    save(data_path, 'ya', '-v7.3');
end
