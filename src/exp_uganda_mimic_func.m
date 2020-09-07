function exp_uganda_mimic_func(pno)
% Experiments on uganda data, mimic the missingness
    out_path = '../outputs/uganda_mimic';
    if ~exist(out_path, 'dir')
        mkdir(out_path);
    end

    % set parameters
    dataset.data_path = '../../data/QC-ed/';
    dataset.phe_no = pno;
    dataset.loader = @dataset_uganda;
    
    pars.np = 1778;
    pars.rx = 1;
    pars.ry = 1778;
    pars.stId = 7;
    pars.test = 1;
    pars.uv = 'u';
    pars.ratio = 1;
    pars.T_BS = 100000;
    
    % get p-values
    [p_val0All,p_val0,p_val,p_valSemi,Sta,StaSemi,nlSel,nlSelUp] = exp_exploration_mimic(dataset,pars);
    
    % plots
    h = figure;
    set(gca,'fontsize', 20);
    id = nlSelUp + nlSel;
    h1 = errorbar(id,mean(repmat(p_val0(:,1),1,length(id))),std(repmat(p_val0(:,1),1,length(id))),'k','LineWidth',3);
    hold on, h2 = errorbar(id,mean(p_val),std(p_val),'g','LineWidth',3);
    hold on, h3 = errorbar(id,mean(p_val0All),std(p_val0All),'r','LineWidth',3);
    hold on, h4 = errorbar(id,mean(p_valSemi),std(p_valSemi),'b','LineWidth',3);hold off;
%     axis([nlSel id(end) 0 1]);
    hh(1)=h1(1);hh(2)=h2(1);hh(3)=h3(1);hh(4)=h4(1);
    % ax = ancestor(h3,'axes');
    % ax.YAxis.FontSize = 10;
    % ax.XAxis.FontSize = 10;
    legend(hh,{'Only paired data','Our method (Only improve Null Dstr)','Oracle','Our method (improve Null and test stat)'},'FontSize',20);
    xlabel(sprintf('Sample size N (paired data size n=%d)', nlSel),'FontSize',20);
    ylabel('p-value','FontSize',20);
%     h = tightfig(h);
    save(fullfile(out_path,['pheno_' num2str(dataset.phe_no) '_' num2str(pars.ratio) '_simul.mat']),'p_val0','p_val','p_val0All','p_valSemi','Sta','StaSemi','nlSel','nlSelUp');
    savefig(h,fullfile(out_path,['pheno_' num2str(dataset.phe_no) '_' num2str(pars.ratio) '_simul.fig']));
