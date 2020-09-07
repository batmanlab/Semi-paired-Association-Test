function [Kxa, Kya, Kza, Kza_eig] = dataset_copd_lle(pars)
% COPD dataset

phe_no = pars.phe_no;
data_path = pars.data_path;
phe_path = fullfile(data_path, 'LLE.txt');
if ~exist('df','var')
    % load phenotype data
    df = readtable(phe_path);
    % load genotype data
    GRMFile = fullfile(data_path, 'CG10kNHWhg19clean_filt.grm');
    GRMid = fullfile(data_path, 'CG10kNHWhg19clean_filt.grm.id');
    [Nsubj, subjId, K] = ParseGRM(GRMFile, GRMid);
end

% extract paired data, x - phenotype, y - genotype
[pairId, pairId_Phe, pairId_Gen] = intersect(table2array(df(:,1)), subjId);
dfArray = table2array(df(:,phe_no));
xa = dfArray(pairId_Phe,:);

if size(phe_no,2)==1 && phe_no==15
    xa = log(xa);
end
if size(phe_no,2)==1
    xa = zscore(xa);
end

% removing confounding
if size(phe_no,2)==1
    za = table2array(df(pairId_Phe,[4,12,13]));
else
    za = table2array(df(pairId_Phe, [102,103,104]));
end
za_eig = table2array(df(pairId_Phe,17:22));

% calculate kernel matrix
Kxa = xa*xa';
if size(phe_no,2)~=1
    pairId_All = 1:size(K,1);
    pairId_All(pairId_Gen) = [];
    pairId_Gen = [pairId_Gen; pairId_All'];
end
Kya = K(pairId_Gen,pairId_Gen);
za = zscore(za);
Kza = za*za';
Kza_eig = za_eig*za_eig';
