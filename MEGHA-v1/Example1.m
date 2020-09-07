%% Examples

%% MEGHA.m
data_path = '/home/gongm/sdc/gitProjects/unpaired/data/copd_geneExp/';
PhenoFile = 'phenotype.txt';
header = 1;
CovFile = 'gender.txt';
% CovFile = '';
delimiter = ',';
GRMFile = fullfile(data_path, 'CG10kNHWhg19clean_filt.grm');
GRMid = fullfile(data_path, 'CG10kNHWhg19clean_filt.grm.id');
Nperm = 0;
WriteStat = 1;
OutDir = './';
%
[Pval, h2, SE, PermPval, PermFWEcPval, Nsubj, Npheno, Ncov] = MEGHA(PhenoFile, header, CovFile, delimiter, GRMFile, GRMid, Nperm, WriteStat, OutDir);
%% MEGHAmat.m
% load Pheno; load Cov; load K;
% Nperm = 1000;
% %
% [Pval, h2, SE, PermPval, PermFWEcPval, Nsubj, Npheno, Ncov] = MEGHAmat(Pheno, Cov, K, Nperm);
%% MEGHAsurf.m
% SurfDir = '/path/SUBJECT_DIR/';
% ImgSubj = 'path/ImgID.txt';
% ImgFileLh = 'lh.fsaverage.thickness.fwhm20.mgh';
% ImgFileRh = 'rh.fsaverage.thickness.fwhm20.mgh';
% FSDir = '/path/freesurfer/';
% CovFile = '/path/qcovar.txt';
% delimiter = '\t';
% GRMFile = '/path/GRM.grm';
% GRMid = '/path/GRM.grm.id';
% WriteImg = 1;
% OutDir = '/OutPutDirectory/';
% Nperm = 1000;
% Pthre = 0.01;
% %
% [PvalLh, PvalRh, h2Lh, h2Rh, SE, ClusPLh, ClusPRh, PeakLh, ClusLh, ClusidLh, PeakRh, ClusRh, ClusidRh, Nsubj, NvetLh, NvetRh, Ncov] = ...
%     MEGHASurf(SurfDir, ImgSubj, ImgFileLh, ImgFileRh, FSDir, CovFile, delimiter, GRMFile, GRMid, WriteImg, OutDir, Nperm, Pthre);
%% MEGHAsurfmat.m
% load PhenoSurf; load Cov; load K;
% FSDir = '/path/freesurfer/';
% Nperm = 1000;
% Pthre = 0.01;
% %
% [PvalLh, PvalRh, h2Lh, h2Rh, SE,  ClusPLh, ClusPRh, PeakLh, ClusLh, ClusidLh, PeakRh, ClusRh, ClusidRh, Nsubj, NvetLh, NvetRh, Ncov] = ...
%     MEGHASurfmat(FSDir, SurfLh, SurfRh, Cov, K, Nperm, Pthre);
%%