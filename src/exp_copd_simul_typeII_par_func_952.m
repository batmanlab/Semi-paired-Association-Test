function exp_copd_simul_typeII_par_func_952(h2, sig, ydim, iter)
% COPD simulated phenotype, type II error evaluation
addpath('../../MEGHA-v1/');
addpath(genpath('../../KCI-test/algorithms'));
addpath(genpath('../../KCI-test/gpml-matlab'));
out_path = '../outputs/COPD_Simulation_typeII/';
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

%for h2=[0.1 0.3 0.5]
    %for sig=[0.5 1]
        %for ydim=[20,60,100]
% dataset configuration
dataset.name = 'dataset_copd_simul_mine';
dataset.h2 = h2; % heritability
dataset.sigma = 1; % noise std from low dimension phenotype to high dimension
dataset.sigma_ny = sqrt(sig);
dataset.ydim = ydim; % phenotype dimension
dataset.yzdim = 10;
dataset.N = 6670;
dataset.T = 1000;
dataset.data_dir = '../../data/COPD_Simulation/'; % data folder.

% precision = sprintf('h2_%%.%dfS2_%%.%dfydim_%%d', digitNum(dataset.h2), digitNum(dataset.sig));
% output_subdir = sprintf(precision, dataset.h2, dataset.sig, dataset.ydim);
%if exist(fullfile(out_path,[output_subdir,'.txt']),'file')
    %continue;
%else
    %fid = fopen(fullfile(out_path,[output_subdir,'.txt']), 'w');
    %fclose(fid);
%end

% algorithm parameters
pars.It_max = 951;
pars.ids = [4000, 5000, 6000];
pars.iter = iter;
pars.np = 3000; % paired data
pars.rx = 10; 
pars.ry = pars.np;         
pars.stId = 1;
pars.test = 1;
pars.uv = 'u';
pars.method = 'LMM_ScoreTest';
pars.thresh = 0.05;

% precision = sprintf('h2_%%.%dfS2_%%.%dfydim_%%d_%%.2f_%%s', digitNum(dataset.h2), digitNum(dataset.sig));
% output_subdir = sprintf(precision, dataset.h2, dataset.sig, dataset.ydim, pars.thresh, pars.uv);
output_subdir = sprintf('h2_%.2f_sig%.2f_signy%.2f_dzy%d_dy%d_thres%.2f_%s',dataset.h2,dataset.sigma,dataset.sigma_ny,dataset.yzdim,dataset.ydim,pars.thresh,pars.uv);
pars.snapshot_subdir = output_subdir; 
pars.out_path = out_path;

% compute p-values
tic;
if ~exist(fullfile(out_path, [output_subdir, '_' num2str(iter) '.mat']),'file')
    [p_val0_1, p_val0_2, p_val_1, p_valSemi_1, Sta_1, StaSemi_1] = simul_typeII_par(dataset, pars);
end
toc;
resFiles = dir([out_path '*_dzy' num2str(dataset.yzdim) '_dy' num2str(dataset.ydim) '*.mat'])
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
