function [p_val0,p_val,p_valSemi,Sta,StaSemi,nl,nlSelUp] = exp_exploration_real(dataset,pars)
% Experiments on real dataset with missing data. Retunrs p-values for all compared methods.
addpath('../../MEGHA-v1/');
addpath(genpath('../../KCI-test/algorithms'));
addpath(genpath('../../KCI-test/gpml-matlab'));

[Kxa, Kya, Kza,~,nl,dimx,dimy] = dataset.loader(dataset);
nax = size(Kxa,1);
nay = size(Kya,1);
naz = size(Kza,1);
if nax-nl==0
    nlSelUp_x = zeros(1,4);
else
    nlSelUp_x = 0:floor((nax-nl)/3):nax-nl;
end
if nay-nl==0
    nlSelUp_y = zeros(1,4);
else
%     nlSelUp_y = 0:floor((nay-nl)/3):nay-nl;
    nlSelUp_y = zeros(1,4);
end
minLen = min(length(nlSelUp_x),length(nlSelUp_y));
nlSelUp_x = nlSelUp_x(1:minLen);
nlSelUp_y = nlSelUp_y(1:minLen);
if nlSelUp_x(end)>nlSelUp_y(end)
    nlSelUp = nlSelUp_x;
else
    nlSelUp = nlSelUp_y;
end

% random select examples
for iter0 = 1:5
    fprintf('iter%d...\n',iter0);
    rng(iter0);
    idSelP_x = randperm(nax-nl);
    idSelP_y = randperm(nay-nl);
    
    for iter = 1:length(nlSelUp_x)
        nlSel = nl;
        idSelAll_x = [1:nlSel,idSelP_x(1:nlSelUp_x(iter))+nlSel];
        idSelAll_y = [1:nlSel,idSelP_y(1:nlSelUp_y(iter))+nlSel];
        
        KxaSel = Kxa(idSelAll_x,idSelAll_x);
        KyaSel = Kya(idSelAll_y,idSelAll_y);
        if isempty(Kza)
            KzaSel_x = [];
            KzaSel_y = [];
        elseif naz==(nax+nay-nl)
            idSelAll_zy = [1:nlSel,idSelP_y(1:nlSelUp_y(iter))+nax];
            KzaSel_x = Kza(idSelAll_x,idSelAll_x);
            KzaSel_y = Kza(idSelAll_zy,idSelAll_zy);
        elseif naz==nax
            KzaSel_x = Kza(idSelAll_x,idSelAll_x);
            KzaSel_y = [];
        else
            KzaSel_y = Kza(idSelAll_y,idSelAll_y);
            KzaSel_x = [];
        end

        % run HSIC tests
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
%         [p_val0(iter0,iter), p_val(iter0,iter), p_valSemi(iter0,iter), Sta(iter0,iter), StaSemi(iter0,iter)] = HSIC_Test(KxaSel, KyaSel, KzaSel_x, KzaSel_y, pars);
        [p_val0(iter0,iter), p_val(iter0,iter), p_valSemi(iter0,iter), Sta(iter0,iter), StaSemi(iter0,iter)] = HSIC_Test(KxaSel, KyaSel, [], [], pars);
    end
end
