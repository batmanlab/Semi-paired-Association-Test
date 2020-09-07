function exp_copd_lle_geneExp_func(fev1_cov)
out_path = '../outputs/copd_lle_geneExp';
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

% set parameter
dataset.data_path = '../../data/copd_geneExp/';
dataset.phe_no = 2:101;
% dataset.phe_no = 15; % FEV1
dataset.fev1 = fev1_cov;
dataset.loader = @dataset_copd_lle_geneExp;

pars.stId = 1;
pars.test = 1;
pars.uv = 'u';
% pars.ratio = 0.99;  
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
axis([nlSel id(end) 0 1]);
hh(1)=h1(1);hh(2)=h2(1);hh(3)=h3(1);
legend(hh,{'Only paired data','Our method (Only improve Null Dstr)','Our method (improve Null and test stat)'},'FontSize',20);
xlabel(sprintf('Sample size N (paired data size n=%d)', nlSel),'FontSize',20);
ylabel('p-value','FontSize',20);
h = tightfig(h);

% save data and plots
save(fullfile(out_path,['pheno_' num2str(dataset.phe_no(1)) '_' num2str(dataset.phe_no(end)) '_' num2str(pars.ratio) '_real.mat']),'p_val0','p_val','p_valSemi','Sta','StaSemi','nlSel','nlSelUp');
savefig(h,fullfile(out_path,['pheno_' num2str(dataset.phe_no(1)) '_' num2str(dataset.phe_no(end)) '_' num2str(pars.ratio) '_real.fig']));
