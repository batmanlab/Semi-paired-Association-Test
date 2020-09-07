function exp_copd_geneExpTseng_func(pheno_id, gene_id)
% Experiments on copd LLE phenotype & blood biomarker data.
% clear;
out_path = '../outputs/copd_geneExp_Tseng';
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

% set parameters
pheno_names = {'pheno_imaging.csv','pheno_breathe.csv','pheno_func.csv','pheno_query.csv','pheno_blood.csv'};
gene_names = {'gene_kegg1.csv','gene_kegg2.csv','gene_kegg3.csv','gene_kegg4.csv'};

dataset.data_path = '../../data/gene_tseng';
% dataset.pheno_name = 'pheno_imaging.csv';
% dataset.pheno_name = 'pheno_blood.csv'; % not correlated
% dataset.pheno_name = 'pheno_breathe.csv';
% dataset.pheno_name = 'pheno_func.csv';
% dataset.pheno_name = 'pheno_query.csv'; % not correlated
% dataset.gene_name = 'gene_kegg1.csv'; 
% dataset.gene_name = 'gene_kegg2.csv'; 
% dataset.gene_name = 'gene_kegg3.csv'; 
% dataset.gene_name = 'gene_kegg4.csv'; 

dataset.pheno_name = pheno_names{pheno_id};
dataset.gene_name = gene_names{gene_id};
dataset.loader = @dataset_copd_geneExpTseng;

switch dataset.pheno_name
    case 'pheno_imaging.csv'
        pars.np = 177; % paired data
        pars.rx = 103;
        pars.ry = 264;
    case 'pheno_breathe.csv'
        pars.np = 114; % paired data
        pars.rx = 15;
        pars.ry = 264;
    case 'pheno_func.csv'
        pars.np = 248; % paired data
        pars.rx = 7;
        pars.ry = 264;
end

pars.stId = 1;
pars.test = 1;
pars.uv = 'u';
pars.ratio = 1;

% get p-values
[p_val0,p_val,p_valSemi,Sta,StaSemi,nlSel,nlSelUp] = exp_exploration_real(dataset,pars);

% plots
h = figure;
set(gca,'fontsize',20);
id = nlSelUp+nlSel;
h1 = errorbar(id,mean(repmat(p_val0(:,1),1,length(id))),std(repmat(p_val0(:,1),1,length(id))),'k','LineWidth',3);
hold on, h2 = errorbar(id,mean(p_val),std(p_val),'g','LineWidth',3);
hold on, h3 = errorbar(id,mean(p_valSemi),std(p_valSemi),'b','LineWidth',3);hold off;
% axis([nlSel id(end) 0 1]);
hh(1)=h1(1);hh(2)=h2(1);hh(3)=h3(1);
legend(hh,{'Only paired data','Our method (Only improve Null Dstr)','Our method (improve Null and test stat)'},'FontSize',20);
xlabel(sprintf('Sample size N (paired data size n=%d)', nlSel),'FontSize',20);
ylabel('p-value','FontSize',20);
% h = tightfig(h);

% save data and plots
save(fullfile(out_path,['pheno_' num2str(dataset.pheno_name) '_' num2str(dataset.gene_name)  '_' num2str(pars.ratio) '_real.mat']),'p_val0','p_val','p_valSemi','Sta','StaSemi','nlSel','nlSelUp');
savefig(h,fullfile(out_path,['pheno_' num2str(dataset.pheno_name) '_' num2str(dataset.gene_name) '_' num2str(pars.ratio) '_real.fig']));
