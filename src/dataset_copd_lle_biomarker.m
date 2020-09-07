function [Kxa,Kya,Kza,Kza_eig,np,dimx,dimy] = dataset_copd_lle_biomarker(pars)
% Load copd LLE (including FEV1) v.s. biomarker data

% phe_no = 2:101;
% data_path = '../../data/copd_biomarker';
% phe_path = fullfile(data_path, 'LLE.txt');

phe_no = pars.phe_no;
data_path = pars.data_path;
phe_path = fullfile(data_path, 'LLE.txt');
bio_cov_path = fullfile(data_path, 'cg10kNhwAllPhenoCovs.txt');
bio_dat_path = fullfile(data_path, 'biomarker_data.csv');

% load LLE phenotype and covariates
if length(phe_no) == 1
    df = readtable(bio_cov_path);
else
    df = readtable(phe_path);
end
bio_cov = readtable(bio_cov_path);
if pars.fev1 == true
    bio_cov = bio_cov(:,{'FID','age_enroll','bmi','ats_packyears','fev1pp_utah','gender'});
else
    bio_cov = bio_cov(:,{'FID','age_enroll','bmi','ats_packyears','gender'});
end

% load blood biomarker and covariates
bio_dat = readtable(bio_dat_path);

% match ids
[pairId, pairId_Phe, pairId_Cov] = intersect(table2array(df(:,1)), table2array(bio_cov(:,1)));
df = df(pairId_Phe,:);
bio_cov = bio_cov(pairId_Cov,:);
% df(1:5,1:8)
% bio_cov(1:5,:)

[pairId, pairId_Phe, pairId_Bio] = intersect(table2array(df(:,1)), table2array(bio_dat(:,2)));

% convert table to array, for phenotype LLE
pairId_All = 1:size(df,1);
pairId_All(pairId_Phe) = [];
pairId_Phe = [pairId_Phe; pairId_All'];
dfArray = table2array(df(:,phe_no));
xa = dfArray(pairId_Phe,:);
za = table2array(bio_cov(:,2:end));
za = za(pairId_Phe,:);

% for blood biomarker, fill in nan values
pairId_All1 = 1:size(bio_dat,1);
pairId_All1(pairId_Bio) = [];
pairId_Bio = [pairId_Bio; pairId_All1'];
ya = table2array(bio_dat(:,3:end));
ya = ya(pairId_Bio,:);
ya_zero = ya;
ya_zero(isnan(ya)) = 0;
ya_mean = repmat(mean(ya_zero),size(ya,1),1);
ya_mean = ya_mean.*double(isnan(ya));
ya = ya_zero+ya_mean;

% construct kernel matrix
Kxa = xa*xa';
Kya = ya*ya';
Kza = za*za';
% Kza = [];
Kza_eig = [];
np = length(pairId);
dimx = size(xa,2);
dimy = size(ya,2);