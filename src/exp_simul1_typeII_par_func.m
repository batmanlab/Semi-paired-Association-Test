function exp_simul1_typeII_par_func(h2, sigma_n, dim, iter)
% Experiments on simulated data, type II error evaluation, for parallel computation on psc 
addpath(genpath('../../KCI-test/algorithms'));
addpath(genpath('../../KCI-test/gpml-matlab'));
out_path = '../outputs/simul1_typeII/';
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

%for sigma_n = [0.1 0.5]
    %for dim = [20 50 100]        
% dataset configuration
dataset.name = 'dataset_simul1';
dataset.h2 = h2; % heritability
dataset.sigma_nx = sqrt(sigma_n); % noise variance from low dimension genotype to high dimension
dataset.sigma_ny = sqrt(sigma_n); % noise variance from low dimension phenotype to high dimension
dataset.sigma = 1;
dataset.xzdim = 10; % genotype hidden dimension
dataset.yzdim = 10; % phenotype hidden dimension
dataset.xdim = dim; % genotype dimension
dataset.ydim = dim; % phenotype dimension
dataset.N = 1000; % sample size
dataset.T = 1000; % number of repetitions
dataset.LMM = 0;
% dataset.data_dir = '/home/gongm/sdc/Data/unpaired/LMM_Simulation/'; % data folder.
dataset.data_dir = '../../data/LMM_Simulation/'; % data folder.

% algorithm parameters
pars.It_max = 1000;
pars.ids = [200, 500, 1000];
pars.iter = iter;
pars.np = 100; % paired data
pars.rx = dataset.yzdim; % pheno
pars.ry = dataset.xzdim; % geno
pars.stId = 1;
pars.test = 1;
pars.uv = 'u';
pars.method = 'HSIC_Test';
pars.thresh = 0.05;

output_subdir = sprintf('h2_%.2f_sig%.2f_signx%.2f_signy%.2f_dzx%d_dzy%d_dx%d_dy%d_N%d_LMM%d_thres%.2f_%s', ...
            dataset.h2, dataset.sigma, dataset.sigma_nx, dataset.sigma_ny, dataset.xzdim, ...
            dataset.yzdim, dataset.xdim, dataset.ydim, dataset.N, dataset.LMM, pars.thresh, pars.uv);
pars.snapshot_subdir = output_subdir;        
pars.out_path = out_path;
      
% compute p-values
if ~exist(fullfile(out_path, [output_subdir, '_' num2str(iter) '.mat']),'file')
    [p_val0_1, p_val0_2, p_val_1, p_valSemi_1, Sta_1, StaSemi_1] = simul_typeII_par(dataset, pars);
end

resFiles = dir([out_path '*_dx' num2str(dim) '_dy' num2str(dim) '_N' num2str(dataset.N) '*.mat'])
if length(resFiles) ~= pars.It_max
    fprintf('file size %d\n', length(resFiles));
    return;               
end
%resFiles = dir(out_path);
%if length(resFiles) < pars.It_max+2
    %fprintf('file size %d\n', length(resFiles));
    %return;
%end

% summarize p-values, if all iterations completed
p_val0_1A = zeros(pars.It_max,length(pars.ids));
p_val0_2A = zeros(pars.It_max,length(pars.ids)+1);
p_valA = zeros(pars.It_max,length(pars.ids));
p_valSemiA = zeros(pars.It_max,length(pars.ids));
StaA = zeros(pars.It_max,length(pars.ids));
StaSemiA = zeros(pars.It_max,length(pars.ids));
for i = 1:length(resFiles)
    load(fullfile(out_path, [output_subdir, '_' num2str(i) '.mat']));
    p_val0_1A(i,:) = p_val0_1;
    p_val0_2A(i,:) = p_val0_2;
    p_valA(i,:) = p_val_1;
    p_valSemiA(i,:) = p_valSemi_1;
    StaA(i,:) = Sta_1;
    StaSemiA(i,:) = StaSemi_1;
end

% compute type II error (power)
p_val0_1 = p_val0_1A;
p_val0_2 = p_val0_2A;
p_val = p_valA;
p_valSemi = p_valSemiA;
Sta = StaA;
StaSemi = StaSemiA;

ids = pars.ids;
It_max = pars.It_max;
thresh = pars.thresh;
TypeII_Sup_All = zeros(1,length(ids)+1);
TypeII_Oracle_All = zeros(1,length(ids)+1);
TypeII_SemiNull_All = zeros(1,length(ids)+1);
TypeII_SemiDR_All = zeros(1,length(ids)+1);

for i = 1:length(ids)+1
    Error2_bs5(1) = sum(p_val0_1(1:It_max,1)>thresh)/It_max;
    Error2_bs5(2) = sum(p_val0_2(1:It_max,i)>thresh)/It_max;
    if i==1
        Error2_bs5(3) = sum(p_val0_1(1:It_max,1)>thresh)/It_max;
    else
        Error2_bs5(3) = sum(p_val(1:It_max,i-1)>thresh)/It_max;
    end
    if i==1
        Error2_bs5(4) = sum(p_val0_1(1:It_max,1)>thresh)/It_max;
    else
        Error2_bs5(4) = sum(p_valSemi(1:It_max,i-1)>thresh)/It_max;
    end
    TypeII_Sup_All(i) = Error2_bs5(1);
    TypeII_Oracle_All(i) = Error2_bs5(2);
    TypeII_SemiNull_All(i) = Error2_bs5(3);
    TypeII_SemiDR_All(i) = Error2_bs5(4);
end

h = figure;
idx = [pars.np,pars.ids];
h1 = plot(idx,1-TypeII_Sup_All,'k','LineWidth',3);
hold on, h2 = plot(idx,1-TypeII_Oracle_All,'r','LineWidth',3);
hold on, h3 = plot(idx,1-TypeII_SemiNull_All,'b','LineWidth',3);
hold on, h4 = plot(idx,1-TypeII_SemiDR_All,'m','LineWidth',3);
axis([min(idx) max(idx) 0 1]);
hh(1)=h1(1);hh(2)=h2(1);hh(3)=h3(1);hh(4)=h4(1);
legend(hh,{'Only paired data','Oracle method', 'Our method Null','Our method DR'},'FontSize',20);
xlabel(sprintf('Sample size N (pair data size n=%d)',pars.np),'FontSize',20);
ylabel('Test Power','FontSize',20);
h = tightfig(h);
savefig(h,fullfile(out_path,[output_subdir,'.fig']));
save(fullfile(out_path,[output_subdir,'_final.mat']),'TypeII_Sup_All','TypeII_Oracle_All','TypeII_SemiNull_All','TypeII_SemiDR_All','p_val0_1','p_val0_2','p_val','p_valSemi','Sta','StaSemi');
