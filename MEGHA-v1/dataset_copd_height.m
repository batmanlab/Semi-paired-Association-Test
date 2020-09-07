function [xa, Kya, za] = dataset_copd_height(pars)
% COPD dataset

phe_no = pars.phe_no;
data_path = pars.data_path;
phe_path = fullfile(data_path, 'cg10kNhwAllPhenoCovs.txt');
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
% xa = zscore(xa);

% removing confounding
za = table2array(df(pairId_Phe,[4,12,17:22]));
za_eig = table2array(df(pairId_Phe,17:22));

% calculate kernel matrix
Kxa = xa*xa';
if size(phe_no,2)~=1
    pairId_All = 1:size(K,1);
    pairId_All(pairId_Gen) = [];
    pairId_Gen = [pairId_Gen; pairId_All'];
end
Kya = K(pairId_Gen,pairId_Gen);
% za = zscore(za);