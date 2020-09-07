function exp_copd_simul_typeI_par_func(h2, sig, ydim, iter)
% COPD simulated phenotype, type I error evaluation, using my own data
% generator
addpath('../../MEGHA-v1/');
addpath(genpath('../../KCI-test/algorithms'));
addpath(genpath('../../KCI-test/gpml-matlab'));
out_path = '../outputs/COPD_Simulation_typeI/';
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

%for sig=[0.5,1]
    %for ydim=[20,60,100]
% dataset configuration
dataset.name = 'dataset_copd_simul_mine';
dataset.h2 = 0; % heritability
dataset.sigma = 1; 
dataset.sigma_ny = sqrt(sig); % noise std from low dimension phenotype to high dimension
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
pars.It_max = 1000;
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

% precision = sprintf('h2_%%.%dfS2_%%.%dfydim_%%d_%%.2f_%%s', digitNum(dataset.h2), digitNum(dataset.sigma_ny));
% output_subdir = sprintf(precision, dataset.h2, dataset.sigma_ny, dataset.ydim, pars.thresh, pars.uv);
output_subdir = sprintf('h2_%.2f_sig%.2f_signy%.2f_dzy%d_dy%d_thres%.2f_%s',dataset.h2,dataset.sigma,dataset.sigma_ny,dataset.yzdim,dataset.ydim,pars.thresh,pars.uv);
pars.snapshot_subdir = output_subdir;        
pars.out_path = out_path;

% compute p-values
tic;
if ~exist(fullfile(out_path, [output_subdir, '_' num2str(iter) '.mat']),'file')
    [p_val0, p_val, p_valSemi, Sta, StaSemi] = simul_typeI_par(dataset, pars);
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
p_val0A = zeros(pars.It_max,length(pars.ids));
p_valA = zeros(pars.It_max,length(pars.ids));
p_valSemiA = zeros(pars.It_max,length(pars.ids));
StaA = zeros(pars.It_max,length(pars.ids));
StaSemiA = zeros(pars.It_max,length(pars.ids));
for i = 1:length(resFiles)
    load(fullfile(out_path, [output_subdir, '_' num2str(i) '.mat']));
    p_val0A(i,:) = p_val0;
    p_valA(i,:) = p_val;
    p_valSemiA(i,:) = p_valSemi;
    StaA(i,:) = Sta;
    StaSemiA(i,:) = StaSemi;
end

% compute type I error
p_val0 = p_val0A;
p_val = p_valA;
p_valSemi = p_valSemiA;
Sta = StaA;
StaSemi = StaSemiA;

TypeI_Sup_All = zeros(1,length(pars.ids));
TypeI_SemiNull_All = zeros(1,length(pars.ids));
TypeI_SemiDR_All = zeros(1,length(pars.ids));
for i = 1:length(pars.ids)
    Error1_bs5(1) = sum(p_val0(1:pars.It_max,i)<pars.thresh)/pars.It_max;
    Error1_bs5(2) = sum(p_val(1:pars.It_max,i)<pars.thresh)/pars.It_max;
    Error1_bs5(3) = sum(p_valSemi(1:pars.It_max,i)<pars.thresh)/pars.It_max;
    TypeI_Sup_All(i) = Error1_bs5(1);
    TypeI_SemiNull_All(i) = Error1_bs5(2);
    TypeI_SemiDR_All(i) = Error1_bs5(3);
end

h = figure;
idx = [pars.np,pars.ids];
idxLen = length(idx);
h1 = plot(idx,repmat(TypeI_Sup_All(1),1,idxLen),'k','LineWidth',3);
hold on, h2 = plot(idx,[TypeI_Sup_All(1),TypeI_SemiNull_All(1:end)],'b','LineWidth',3);
hold on, h3 = plot(idx,[TypeI_Sup_All(1),TypeI_SemiDR_All],'m','LineWidth',3);
hold on, h4 = plot(idx,pars.thresh*ones(1,idxLen),'r--','LineWidth',3);
axis([min(idx) max(idx) 0 0.1]);
hh(1)=h1(1);hh(2)=h2(1);hh(3)=h3(1);hh(4)=h4(1);
legend(hh,{'Only paired data','Our method', 'Our method DR','True Type I error'},'FontSize',20);
xlabel(sprintf('Sample size N (pair data size n=%d)',pars.np),'FontSize',20);
ylabel('Type I Error','FontSize',20);
h = tightfig(h);
save(fullfile(out_path,[output_subdir,'_final.mat']),'TypeI_Sup_All','TypeI_SemiNull_All','TypeI_SemiDR_All','p_val0','p_val','p_valSemi','Sta','StaSemi');
savefig(h,fullfile(out_path,[output_subdir,'.fig']));
