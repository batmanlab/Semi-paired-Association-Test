function [Kxa,Kya,Kza,Kza_eig,np,dimx,dimy] = dataset_copd_geneExpTseng(pars)
% gene expression and clinic indicators
data_path = pars.data_path;
phe_path = fullfile(data_path, pars.pheno_name);
phe_cov_path = fullfile(data_path, 'pheno_cov.csv');

% load phenotypes and covariates
df = readtable(phe_path);
df_cov = readtable(phe_cov_path);

% load gene
gene_dat_path = fullfile(data_path, pars.gene_name);
gene_dat = readtable(gene_dat_path);
% gene_dat(1:5,1:8)
% size(gene_dat)

% convert table to array, for phenotype and covariates
dfArray = table2array(df(:,2:end));
dfArray = cellfun(@tonan,dfArray,'un',0);
dfArray = cell2mat(cellfun(@str2num,dfArray,'un',0));
dfArrayRowsum = sum(dfArray,2);
xa = dfArray(~isnan(dfArrayRowsum),:);
% size(xa)

% needs further improvement, missing values
dfArray_cov = table2array(df_cov(:,3:end));
dfArray_cov = cellfun(@tonan,dfArray_cov,'un',0);
dfArray_cov = cell2mat(cellfun(@str2num,dfArray_cov,'un',0));
for i = 1:size(dfArray_cov,2)
    x_tmp = dfArray_cov(:,i);
    x_tmp(isnan(x_tmp))=mean(x_tmp(~isnan(x_tmp)));
    dfArray_cov(:,i) = x_tmp;
end
za = dfArray_cov(~isnan(dfArrayRowsum),:);
% size(za)

% for gene expression
geneArray = table2array(gene_dat(:,2:end));
ya = [geneArray(~isnan(dfArrayRowsum),:);geneArray(isnan(dfArrayRowsum),:)];
% size(ya)

% construct kernel matrix
Kxa = xa*xa';
Kya = ya*ya';
Kza = za*za';
% Kza = [];
Kza_eig = [];
np = size(xa,1);
if strcmp(pars.pheno_name, 'pheno_func.csv')
    np = 200;
end
if strcmp(pars.pheno_name, 'pheno_imaging.csv')
    np = 60;
end
if strcmp(pars.pheno_name, 'pheno_breathe.csv')
    np = 50;
end
% np = 60;% for imaging features
% np = 200; % for functional features
dimx = size(xa,2);
dimy = size(ya,2);