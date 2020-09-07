function [p_val0, p_val, p_valSemi, Sta, StaSemi] = LMM_ScoreTest(Kx, Ky, Kzx, Kzy, pars)
% Variance Component Score Test (VCST) using semi-paired data. X-phenotype, Y-genotype
% Inputs: 
%       Kx - kernel matrix on x (NxN, the first np x np block contains paired data)
%       Ky - kernel matrix on y (MxM, the first np x np block contains paired data)
%       Kzx - kernel matrix on covariate z (for x, empty if not exist)
%       Kzy - kernel matrix on covariate z (for y, empty if not exist)
%       pars - hyperparameters
%           .np - number of paired data
%           .rx - reduced dimension of x, default pars.np
%           .ry - reduced dimension of y, default pars.np
%           .stId - remove the top (stId-1) eigenvalues (for genotype data), defaut 1
%           .T_BS - bootstrap sample size, default 10000
%           .uv - 'u' or 'v' statistics, default 'u'
%           .test - compute p-values if true, default 1
% Outputs:
%       p_val0 - p value of the original VCST using only paired data
%       p_val - p value of our Semi-paired test (SAT), only improve null distribution
%       p_valSemi - p value of our SAT, improve both test statistics and null distribution
%       Sta - test statistic of the original VCST using only paired data
%       StaSemi - test statistic of our SAT

Tx = length(Kx); % the sample size
Ty = length(Ky);
rx = pars.rx;
ry = pars.np;
Tn = pars.np;
stId = pars.stId;

T_BS = 10000;
Thresh = 1E-8;

H = eye(Tn)-ones(Tn,Tn)/Tn; % for centering of the data in feature space
Hzx = eye(Tx)-ones(Tx)/Tx;
if ~isempty(Kzx)
    Kzx = Hzx*Kzx*Hzx;
    epsilon = 1e-5;
    Rx = epsilon*inv((Kzx+epsilon*eye(Tx)));
    Kx = Rx*Kx*Rx;
end
Kyl = Ky(1:Tn,1:Tn);
if ~isempty(Kzy)
    Kzyl = Kzy(1:Tn,1:Tn);
    Kzyl = H*Kzyl*H;
    epsilon = 1e-5;
    Ry = epsilon*inv((Kzyl+epsilon*eye(Tn)));
    Kyl = Ry*Kyl*Ry;
end

Kxl = Kx(1:Tn,1:Tn);
Kxlc = H*Kxl*H;
% Kyl = Ky(1:Tn,1:Tn);
% Kylc = H*Kyl*H;

% if ~isempty(Kz_eig)
%     Kzl_eig = Kz_eig(1:Tn,1:Tn);
%     Kzl_eig = H*Kzl_eig*H;
%     epsilon = 1e-5;
%     R_eig = epsilon*inv((Kzl_eig+epsilon*eye(Tn)));
% %     Kxl = R_eig*Kxl*R_eig;
%     Kyl = R_eig*Kyl*R_eig;
% end

% remove top eigs of Y (SNPs)
KylC = H*Kyl*H;
[uy1, dy1] = getEigen(KylC);

dy1T = dy1;
dy1T(1:stId-1,1:stId-1) = 0;
KylPC = uy1*dy1T*uy1';

if strcmp(pars.uv,'v')
    Sta = trace(Kxlc*KylPC);
elseif strcmp(pars.uv,'u')
    Sta = trace(Kxlc*KylPC) - 1/Tn*trace(Kxlc)*trace(KylPC);
end

if pars.test == 1
    p_val0 = cal_pval(Sta, Kxlc, KylPC, Tn, T_BS, Tn, Tn, Thresh, 1, pars.uv);
else
    p_val0 = nan;
end
% Sta0 = trace(Kxl * Kyl)
% p_val00 = cal_pval(Sta0, Kxl, Kyl, Bootstrap, Approximate, Tn, T_BS, Tn, Tn, Thresh)

if Tn == Tx 
    p_val = p_val0;
    p_valSemi = p_val0;
    StaSemi = Sta;
    return;
end

Hx = eye(Tx)-ones(Tx)/Tx;
KxC = Hx*Kx*Hx;
[ux1, dx1] = getEigen(KxC);

if strcmp(pars.uv,'v')
    Sta = trace(Kxlc*KylPC);
elseif strcmp(pars.uv,'u')
    Sta = trace(Kxlc*KylPC) - 1/Tx*trace(KxC)*trace(KylPC);
end

% original HSIC + modified null distr
if pars.test == 1
    p_val = cal_pval(Sta, KxC, KylPC, Tn, T_BS, Tx, Tn, Thresh, stId, pars.uv);
else
    p_val = nan;
end

if rx == Tx 
    p_valSemi = p_val;
    StaSemi = Sta;
    return;
end

% modified HSIC + modified null distr
reg = 1e-5;
KxS = Kx(1:Tn,:)*Hx*ux1(:,1:rx)/(dx1(1:rx,1:rx)+reg*eye(rx))*ux1(:,1:rx)'*Hx*Kx(:,1:Tn);
KxSC = H*KxS*H;
if strcmp(pars.uv,'v')
    StaSemi = trace(KxSC*KylPC);
elseif strcmp(pars.uv,'u')
    StaSemi = trace(KxSC*KylPC) - 1/Tx*trace(dx1(1:rx,1:rx))*trace(KylPC);
end

if pars.test == 1
    p_valSemi = cal_pval(StaSemi, KxC, KylPC, Tn, T_BS, rx, Tn, Thresh, stId, pars.uv);
else 
    p_valSemi = nan;
end
end


function p_val = cal_pval(Sta, Kx, Ky, Tn, T_BS, Num_eigX, Num_eigY, Thresh, stId, uv)
% calculate the eigenvalues that will be used later
% Due to numerical issues, Kx and Ky may not be symmetric:
Tx = length(Kx);
Ty = length(Ky);
[eig_Kx, eivx] = eigdec((Kx+Kx')/2,Num_eigX);
[eig_Ky, eivy] = eigdec((Ky+Ky')/2,Num_eigY);
Num_eigY = Num_eigY-stId+1;
% calculate Cri...
% first calculate the product of the eigenvalues
eig_Kx = eig_Kx(1:Num_eigX)/Tx;
eig_Ky = eig_Ky(stId:end);
eig_prod = stack( (eig_Kx * ones(1,Num_eigY)) .* (ones(Num_eigX,1) * eig_Ky'));
II = find(eig_prod > max(eig_prod) * Thresh);
eig_prod = eig_prod(II); %%% new method

% use mixture of F distributions to generate the Null dstr
if strcmp(uv, 'u')
    f_rand1 = chi2rnd(1,length(eig_prod),T_BS) - 1;
elseif strcmp(uv, 'v')
    f_rand1 = chi2rnd(1,length(eig_prod),T_BS);
end
Null_dstr = eig_prod' * f_rand1; %%%%Problem
p_val = sum(Null_dstr>Sta)/T_BS;
end
