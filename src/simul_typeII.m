function [TypeII_Sup_All, TypeII_Oracle_All, TypeII_SemiNull_All, TypeII_SemiDR_All, p_val0_1, p_val0_2, p_val_1, p_valSemi_1, Sta_1, StaSemi_1] = simul_typeII(dataset, pars)
% Calculate type I error from simulated data

% load data
[xa, ya] = feval(dataset.name, dataset);

ids = pars.ids;
It_max = min(pars.It_max,size(ya,3));
p_val0_1 = zeros(It_max,length(ids));
p_val_1 = zeros(It_max,length(ids));
p_valSemi_1 = zeros(It_max,length(ids));
Sta_1 = zeros(It_max,length(ids));
StaSemi_1 = zeros(It_max,length(ids));
p_val0_2 = zeros(It_max,length(ids)+1);
p_val_2 = zeros(It_max,length(ids)+1);
p_valSemi_2 = zeros(It_max,length(ids)+1);
Sta_2 = zeros(It_max,length(ids)+1);
StaSemi_2 = zeros(It_max,length(ids)+1);
np = pars.np;
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

    pars.np = np;
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
        [p_val0_1(iter,i), p_val_1(iter,i), p_valSemi_1(iter,i), Sta_1(iter,i), StaSemi_1(iter,i)] = feval(pars.method, Kx, Ky, Kzx, Kzy, pars);
    end
    
    % oracle method
    idsAg = [np, ids];
    for i = 1:length(ids)+1
%         N = ids(end);
        pars.np = idsAg(i);
        N = pars.np;
        % testing...
        Kzx = [];
        Kzy = [];
        if strcmp(pars.method, 'LMM_ScoreTest')
            Ky = Kxa_i(1:pars.np,1:pars.np);
        else
            Ky = Kxa_i(1:N,1:N);
        end
        Kx = Kya_i(1:N,1:N);
        [p_val0_2(iter,i), p_val_2(iter,i), p_valSemi_2(iter,i), Sta_2(iter,i), StaSemi_2(iter,i)] = feval(pars.method, Kx, Ky, Kzx, Kzy, pars);
    end
    save(snapshot_path, 'iter', 'p_val0_1','p_val_1','p_valSemi_1','Sta_1','StaSemi_1',...
        'p_val0_2','p_val_2','p_valSemi_2','Sta_2','StaSemi_2');
end

thresh = pars.thresh;
TypeII_Sup_All = zeros(1,length(ids)+1);
TypeII_Oracle_All = zeros(1,length(ids)+1);
TypeII_SemiNull_All = zeros(1,length(ids)+1);
TypeII_SemiDR_All = zeros(1,length(ids)+1);

for i = 1:length(ids)+1
    Error2_bs5(1) = sum(p_val0_1(1:It_max,1)>thresh)/It_max;
    Error2_bs5(2) = sum(p_val0_2(1:It_max,i)>thresh)/It_max;
    if i==1
        Error2_bs5(3) = sum(p_val0_1(1:It_max,1)>thresh)/It_max;
    else
        Error2_bs5(3) = sum(p_val_1(1:It_max,i-1)>thresh)/It_max;
    end
    if i==1
        Error2_bs5(4) = sum(p_val0_1(1:It_max,1)>thresh)/It_max;
    else
        Error2_bs5(4) = sum(p_valSemi_1(1:It_max,i-1)>thresh)/It_max;
    end
    TypeII_Sup_All(i) = Error2_bs5(1);
    TypeII_Oracle_All(i) = Error2_bs5(2);
    TypeII_SemiNull_All(i) = Error2_bs5(3);
    TypeII_SemiDR_All(i) = Error2_bs5(4);
end
