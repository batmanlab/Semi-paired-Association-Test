function [Kxa,Kya,Kza,Kza_eig,np,dimx,dimy] = dataset_copd_lle_geneExp(pars)
% Load copd LLE vs gene expression data
% phe_no = 2:101;
% data_path = '../../data/copd_geneExp';
% phe_path = fullfile(data_path, 'LLE.txt');

phe_no = pars.phe_no;
data_path = pars.data_path;
phe_path = fullfile(data_path, 'LLE.txt');
phe_cov_path = fullfile(data_path, 'cg10kNhwAllPhenoCovs.txt');

% load LLE phenotype
if length(phe_no) == 1
    df = readtable(phe_cov_path);
else
    df = readtable(phe_path);
end

% load covariates
phe_cov = df(:,102:104);

% load gene expression
gene_dat_path = fullfile(data_path, 'gene_express.csv');
gene_dat = readtable(gene_dat_path);
% size(gene_dat)

% match ids
[pairId, pairId_Phe, pairId_Gene] = intersect(table2array(df(:,1)), table2array(gene_dat(:,2)));

% convert table to array, for phenotype LLE
pairId_All = 1:size(df,1);
pairId_All(pairId_Phe) = [];
pairId_Phe = [pairId_Phe; pairId_All'];
dfArray = table2array(df(:,phe_no));
xa = dfArray(pairId_Phe,:);
za = table2array(phe_cov);
za = za(pairId_Phe,:);

% for gene expression
pairId_All1 = 1:size(gene_dat,1);
pairId_All1(pairId_Gene) = [];
pairId_Gene = [pairId_Gene; pairId_All1'];
gene_dat1 = gene_dat(:,3:8);
gene_dat2 = gene_dat(:,9:end);
gene_dat1_array = table2array(gene_dat1);
gene_dat2_array = table2array(gene_dat2);
gene_dat1_array = cell2mat(cellfun(@str2num,gene_dat1_array,'un',0));
ya = [gene_dat1_array,gene_dat2_array];
ya = ya(pairId_Gene,:);

% construct kernel matrix
Kxa = xa*xa';
Kya = ya*ya';
Kza = za*za';
Kza_eig = [];
np = length(pairId);
dimx = size(xa,2);
dimy = size(ya,2);