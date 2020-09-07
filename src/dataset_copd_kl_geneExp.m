function [Kxa,Kya,Kza,Kza_eig,np,dimx,dimy] = dataset_copd_kl_geneExp(pars)
% Load copd LLE vs gene expression data
% phe_no = 2:101;
% data_path = '../../data/copd_geneExp';
% phe_path = fullfile(data_path, 'LLE.txt');

data_path = pars.data_path;
phe_path = fullfile(data_path, 'div.mat');
pheId_path = fullfile(data_path, 'div_id.csv');
phe_cov_path = fullfile(data_path, 'cg10kNhwAllPhenoCovs.txt');

% load KL phenotype
load(phe_path);
Kxa = double(div);
df = readtable(pheId_path);

phe_cov = readtable(phe_cov_path);
if pars.fev1 == true
    phe_cov = phe_cov(:,{'FID','age_enroll','bmi','ats_packyears','fev1pp_utah','gender'});
else
    phe_cov = phe_cov(:,{'FID','age_enroll','bmi','ats_packyears','gender'});
end

% load gene expression
gene_dat_path = fullfile(data_path, 'gene_express.csv');
gene_dat = readtable(gene_dat_path);
% size(gene_dat)

% match ids
[pairId, pairId_Phe, pairId_Cov] = intersect(table2array(df(:,1)), table2array(phe_cov(:,1)));
df = df(pairId_Phe,:);
Kxa = Kxa(pairId_Phe, pairId_Phe);
phe_cov = phe_cov(pairId_Cov,:);

[pairId, pairId_Phe, pairId_Gene] = intersect(table2array(df(:,1)), table2array(gene_dat(:,2)));

% convert table to array, for phenotype LLE
pairId_All = 1:size(df,1);
pairId_All(pairId_Phe) = [];
pairId_Phe = [pairId_Phe; pairId_All'];
Kxa = Kxa(pairId_Phe, pairId_Phe);
za = table2array(phe_cov(:,2:end));
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
% Kxa = xa*xa';
Kya = ya*ya';
Kza = za*za';
Kza_eig = [];
np = length(pairId);
dimx = size(Kxa,2);
dimy = size(ya,2);