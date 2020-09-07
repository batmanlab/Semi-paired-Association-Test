function [TypeI_Sup_All, TypeI_SemiNull_All, TypeI_SemiDR_All, p_val0, p_val, p_valSemi, Sta, StaSemi] = simul_typeI(dataset, pars)
% Calculate type I error from simulated data

% load data
[xa, ya] = feval(dataset.name, dataset);

ids = pars.ids;
It_max = min(pars.It_max,size(ya,3));
p_val0 = zeros(It_max,length(ids));
p_val = zeros(It_max,length(ids));
p_valSemi = zeros(It_max,length(ids));
Sta = zeros(It_max,length(ids));
StaSemi = zeros(It_max,length(ids));
out_path = pars.out_path;
snapshot_subdir = pars.snapshot_subdir;
snapshot_path = fullfile(out_path, [snapshot_subdir,'.mat']);
if exist(snapshot_path, 'file')
    load(snapshot_path);
end
if exist('iter', 'var')
    iter0 = iter+1;
else
    iter0 = 1;
end

% test under different unpaired sample size
for iter = iter0:It_max
    if ~mod(iter,10)
        fprintf('%d replications in finding the type I errors in Case I...\n', iter),
    end
    if size(xa,3)>1
        xa_i = xa(:,:,iter);
    else
        xa_i = xa;
    end
    ya_i = ya(:,:,iter);
    Kya_i = ya_i*ya_i';

    if strcmp(dataset.name, 'dataset_simul1')
        Kxa_i = xa_i*xa_i';
    else
        Kxa_i = xa_i;
    end 

    for i = 1:length(ids)
        N = ids(i);
        % testing...
        Kzx = [];
        Kzy = [];
        if strcmp(pars.method, 'LMM_ScoreTest')
            Ky = Kxa_i(1:pars.np,1:pars.np);
        else
            Ky = Kxa_i(1:N,1:N);
        end
        Kx = Kya_i(1:N,1:N);
        [p_val0(iter,i), p_val(iter,i), p_valSemi(iter,i), Sta(iter,i), StaSemi(iter,i)] = feval(pars.method, Kx, Ky, Kzx, Kzy, pars);
    end
    save(snapshot_path, 'iter', 'p_val0','p_val','p_valSemi','Sta','StaSemi');
end

thresh = pars.thresh;
TypeI_Sup_All = zeros(1,length(ids));
TypeI_SemiNull_All = zeros(1,length(ids));
TypeI_SemiDR_All = zeros(1,length(ids));
for i = 1:length(ids)
    Error1_bs5(1) = sum(p_val0(1:It_max,i)<thresh)/It_max;
    Error1_bs5(2) = sum(p_val(1:It_max,i)<thresh)/It_max;
    Error1_bs5(3) = sum(p_valSemi(1:It_max,i)<thresh)/It_max;
    TypeI_Sup_All(i) = Error1_bs5(1);
    TypeI_SemiNull_All(i) = Error1_bs5(2);
    TypeI_SemiDR_All(i) = Error1_bs5(3);
end
