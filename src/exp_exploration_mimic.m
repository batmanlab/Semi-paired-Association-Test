function [p_val0All,p_val0,p_val,p_valSemi,Sta,StaSemi,nl,nlSelUp] = exp_exploration_mimic(dataset,pars)
% Mimic missinging on real dataset, using all the paired data as a oracle
% method. Retunrs p-values for all compared methods.
addpath('../../MEGHA-v1/');
addpath(genpath('../../KCI-test/algorithms'));
addpath(genpath('../../KCI-test/gpml-matlab'));

[Kxa,Kya,Kza] = dataset.loader(dataset);
nax = size(Kxa,1);
nay = size(Kya,1);
naz = size(Kza,1);
nl = pars.np;
if nax~=nay
    error('Data should be paired');
end
nlSelUp = 0:floor((nax-nl)/3):(nax-nl);
nlSelP = nl + nlSelUp;

% select examples
for iter0 = 1:5
    fprintf('iter%d...\n',iter0);
    rng(iter0);
     idSelP = randperm(nax);
     for iter = 1:length(nlSelP)
         nlSel = nlSelP(iter);
         idSelAll = idSelP(1:nlSel);
         
         KxaSel = Kxa(idSelAll,idSelAll);
         KyaSel = Kya(idSelAll,idSelAll);
         
         if isempty(Kza)
             KzaSel_x = [];
             KzaSel_y = [];
         else
             KzaSel_x = Kza(idSelAll,idSelAll);
             KzaSel_y = Kza(idSelAll,idSelAll);
         end
         % run HISC tests
         pars.np = nlSel; % paired data
         p_val0All(iter0,iter) = HSIC_Test(KxaSel, KyaSel, KzaSel_x, KzaSel_y, pars);
     end
    
    % select examples
    rng(iter0);
    idSelUp = randperm(nax);
    for iter = 1:length(nlSelUp)
        nlSel = nl;
        nlSelUp_xy = nlSelUp(iter);
        idSelAll = [idSelUp(1:nlSel),idSelUp(1+nlSel:nlSelUp_xy+nlSel)];
        
        KxaSel = Kxa(idSelAll,idSelAll);
        KyaSel = Kya(idSelAll,idSelAll);
        if isempty(Kza)
            KzaSel_x = [];
            KzaSel_y = [];
        else
            KzaSel_x = Kza(idSelAll,idSelAll);
            KzaSel_y = Kza(idSelAll,idSelAll);
        end
        % run HISC tests
        pars.np = nlSel; % paired data
%         [ux,dx] = getEigen(KxaSel);
%         [uy,dy] = getEigen(KyaSel);
%         energy_x = cumsum(diag(dx));
%         energy_y = cumsum(diag(dy));
%         pars.rx = find(energy_x>energy_x(end)*pars.ratio);
%         pars.ry = find(energy_y>energy_y(end)*pars.ratio);
%         pars.rx = pars.rx(1);
%         pars.ry = pars.ry(1);

        pars.rx = min(rank(KxaSel),nlSel*pars.ratio);
        pars.ry = min(rank(KyaSel),nlSel*pars.ratio);
        [p_val0(iter0,iter),p_val(iter0,iter),p_valSemi(iter0,iter),Sta(iter0,iter),StaSemi(iter0,iter)] = HSIC_Test(KxaSel,KyaSel,KzaSel_x,KzaSel_y,pars);
    end
end
