function [Kxa, Kya, Kza, Kza_eig, np, dimx, dimy] = dataset_uganda(pars)
% Uganda dataset

phe_no = pars.phe_no;
data_dir = pars.data_path;
% load phenotype data
if ~exist('df','var')
    df = readtable(fullfile(data_dir,'gwas_phenotypes_28Oct14.txt'));
    % load genotype data
%     load('../../data/QC-ed/Uganda_id_table.mat');
%     load('../../data/QC-ed/Uganda_SNP.mat');
%     grm = double(grm);
    GRMFile = fullfile(data_dir, 'uganda_grm.grm');
    GRMid = fullfile(data_dir, 'uganda_grm.grm.id');
    [Nsubj, subjId, grm] = ParseGRM(GRMFile, GRMid);
end
% subjIdTab = struct2table(subjId.ids{1});
% for i = 2:length(subjId.ids)
%     subjIdTab = [subjIdTab; struct2table(subjId.ids{i})];
% end
% subjId = subjIdTab;

% extract paired data, x - phenotype, y - genotype
% [pairId, pairId_Phe, pairId_Gen] = intersect(table2array(df(:,2)), table2array(subjId(:,1)));
[pairId, pairId_Phe, pairId_Gen] = intersect(table2array(df(:,2)), subjId(:,1));

emptyIndex = zeros(length(pairId_Phe),1);
for i=1:length(phe_no)
    dfArray_i = table2array(df(:,phe_no(i)));
    if iscell(dfArray_i) && phe_no(i) < 25
        dfArray_i = cellfun(@str2num,dfArray_i,'un',0);
        dfArray_i = dfArray_i(pairId_Phe);
        dfArray_i1 = cell2mat(dfArray_i');
        emptyIndex = emptyIndex | cellfun(@isempty, dfArray_i);
        dfArray_i(emptyIndex) = {mean(dfArray_i1)};
        dfArray_i = cell2mat(dfArray_i')';
        dfArray(:,i) = dfArray_i;
    elseif ~iscell(dfArray_i) && phe_no(i) < 25
        dfArray_i = dfArray_i(pairId_Phe);
        dfArray(:,i) = dfArray_i;
    else
        dfArray_i = cellfun(@str2num,dfArray_i,'un',0);
        dfArray_i = dfArray_i(pairId_Phe);
        emptyIndex = emptyIndex | cellfun(@isempty, dfArray_i);
        dfArray0{i} = dfArray_i;
    end
end
for i=1:length(phe_no)
    dfArray_i = table2array(df(:,phe_no(i)));
    if iscell(dfArray_i) && phe_no(i) < 25
        continue;
    elseif ~iscell(dfArray_i) && phe_no(i) < 25
        continue;
    else
        dfArray(:,i) = cell2mat(dfArray0{i}(~emptyIndex)')';
    end
end
xa = zscore(dfArray);

% removing confounding
% z = table2array(df(pairId_Phe,[4,12,13,17:22]));
% z = table2array(df(pairId_Phe,[12,13,17:22]));
if max(phe_no) < 25
    Id_Phe = 1:size(dfArray,1);
else
    Id_Phe = [find(emptyIndex==0);find(emptyIndex==1)];
end

dfArray = table2array(df(:,3));
dfArray = cellfun(@sex2num,dfArray,'un',0);
dfArray = cell2mat(dfArray')';
dfArray = dfArray(pairId_Phe);
za = dfArray(Id_Phe);

dfArray = table2array(df(:,4));
dfArray = dfArray(pairId_Phe);
za(:,2) = dfArray(Id_Phe);

% calculate kernel matrix
Kxa = xa*xa';

Kya = grm;
Kya = Kya(pairId_Gen,pairId_Gen);
Kya = Kya(Id_Phe,Id_Phe);

za = zscore(za);
Kza = za*za';
np = size(Kxa, 1);

Kza_eig = [];
dimx = 0;
dimy = 0;
