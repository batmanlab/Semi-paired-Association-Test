function [p_val0_1, p_val0_2, p_val_1, p_valSemi_1, Sta_1, StaSemi_1] = simul_typeII_par(dataset, pars)
% Calculate type II error from simulated data, for parallel computation on psc

% load data
[xa, ya] = feval(dataset.name, dataset);
fprintf('Dataset loaded successfully\n');

iter = pars.iter;
ids = pars.ids;
It_max = 1;
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
snapshot_path = fullfile(out_path, [snapshot_subdir, '_' num2str(iter) '.mat']);

if exist(snapshot_path, 'file')
    load(snapshot_path);
end

% test under different unpaired sample size
fprintf('iter%d\n', iter);
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
    [p_val0_1(i), p_val_1(i), p_valSemi_1(i), Sta_1(i), StaSemi_1(i)] = feval(pars.method, Kx, Ky, Kzx, Kzy, pars);
end

% oracle method
idsAg = [np, ids];
for i = 1:length(ids)+1
    % N = ids(end);
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
    [p_val0_2(i), p_val_2(i), p_valSemi_2(i), Sta_2(i), StaSemi_2(i)] = feval(pars.method, Kx, Ky, Kzx, Kzy, pars);
end
save(snapshot_path, 'iter', 'p_val0_1','p_val_1','p_valSemi_1','Sta_1','StaSemi_1',...
    'p_val0_2','p_val_2','p_valSemi_2','Sta_2','StaSemi_2');
