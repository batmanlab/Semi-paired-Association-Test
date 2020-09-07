Nperm = 0;
dataset.phe_no = 11;
dataset.data_path = '~/sdc/gitProjects/unpaired/data/QC-ed/';
[xa, Kya, za, Kza_eig, np, dimx, dimy] = dataset_uganda_height(dataset);
[u,d]=getEigen(Kya);
d(1:6,1:6)=0;
KyaP = u*d*u';

[Pval_1, h2_1, SE, PermPval_1, PermFWEcPval_1, Nsubj, Npheno, Ncov] = MEGHAmat(xa, ones(size(za,1),1), Kya, Nperm);
[Pval_2, h2_2, SE, PermPval_2, PermFWEcPval_2, Nsubj, Npheno, Ncov] = MEGHAmat(xa, [za(:,1),ones(size(za,1),1)], Kya, Nperm);
[Pval_3, h2_3, SE, PermPval_3, PermFWEcPval_3, Nsubj, Npheno, Ncov] = MEGHAmat(xa, [za(:,2),ones(size(za,1),1)], Kya, Nperm);
[Pval_4, h2_4, SE, PermPval_4, PermFWEcPval_4, Nsubj, Npheno, Ncov] = MEGHAmat(xa, [za(:,1:2),ones(size(za,1),1)], Kya, Nperm);
[Pval_5, h2_5, SE, PermPval_5, PermFWEcPval_5, Nsubj, Npheno, Ncov] = MEGHAmat(xa, [ones(size(za,1),1)], KyaP, Nperm);
[Pval_6, h2_6, SE, PermPval_6, PermFWEcPval_6, Nsubj, Npheno, Ncov] = MEGHAmat(xa, [za,ones(size(za,1),1)], KyaP, Nperm);

fprintf('no confounder h2=%f, pval=%f\n', h2_1, Pval_1);
fprintf('gender confounder h2=%f, pval=%f\n', h2_2, Pval_2);
fprintf('age confounder h2=%f, pval=%f\n', h2_3, Pval_3);
fprintf('gender+age confounder h2=%f, pval=%f\n', h2_4, Pval_4);
fprintf('pc1-6 confounder h2=%f, pval=%f\n', h2_5, Pval_5);
fprintf('gender+age+pc1-6 confounder h2=%f, pval=%f\n', h2_6, Pval_6);
