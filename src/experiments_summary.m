%% Summarization of exploration experiments
%% Biomarker
% 1. LLE 
% 2. KL_div
% exp_copd_lle_biomarker_func(false);
% exp_copd_kl_biomarker_func(false);

%% Gene Expression
% 1. LLE
% 2. KL_div
% exp_copd_lle_geneExp_func(false)
exp_copd_kl_geneExp_func(false);

%% Tseng Gene Expression
% 1. Imaging
% 2. Breathe
% 3. Functional
% 4. Query
% 5. Blood test
% 6. Covariates
% for i=1
%     for j=1:4
%         exp_copd_geneExpTseng_func(i,j);
%     end
% end
%% Uganda Real data
% for i = 25:44
%     exp_uganda_real_func(i);
% end
%% Uganda mimic data
% for i = 6:24
%     exp_uganda_mimic_func(i);
% end
