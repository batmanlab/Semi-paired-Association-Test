function [p_val0, p_val, p_valSemi, Sta, StaSemi] = simul_typeI_par(dataset, pars)
% Calculate type I error from simulated data

% load data
[xa, ya] = feval(dataset.name, dataset);
fprintf('Dataset loaded successfully\n');

iter = pars.iter;
ids = pars.ids;
It_max = 1;
p_val0 = zeros(It_max,length(ids));
p_val = zeros(It_max,length(ids));
p_valSemi = zeros(It_max,length(ids));
Sta = zeros(It_max,length(ids));
StaSemi = zeros(It_max,length(ids));
out_path = pars.out_path;
snapshot_subdir = pars.snapshot_subdir;
snapshot_path = fullfile(out_path, [snapshot_subdir, '_' num2str(iter) '.mat']);

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

for i = 1:length(ids)
    N = ids(i);
    Kzx = [];
    Kzy = [];
    if strcmp(pars.method, 'LMM_ScoreTest')
        Ky = Kxa_i(1:pars.np,1:pars.np);
    else
        Ky = Kxa_i(1:N,1:N);
    end
    Kx = Kya_i(1:N,1:N);
    [p_val0(i), p_val(i), p_valSemi(i), Sta(i), StaSemi(i)] = feval(pars.method, Kx, Ky, Kzx, Kzy, pars);
end
save(snapshot_path, 'iter', 'p_val0','p_val','p_valSemi','Sta','StaSemi');
