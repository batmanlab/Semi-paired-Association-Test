function [p_val0, p_val, p_valSemi, Sta, StaSemi] = HSIC_Test(Kx, Ky, Kzx, Kzy, pars)
% Hilbert Schmidt Independence Criterion (HSIC) test using semi-paired data
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
%           .test - compute p-values if true, default true
% Outputs:
%       p_val0 - p value of the original HSIC using only paired data
%       p_val - p value of our Semi-paired test (SAT), only improve null distribution
%       p_valSemi - p value of our SAT, improve both test statistics and null distribution
%       Sta - test statistic of the original HSIC using only paired data
%       StaSemi - test statistic of our SAT

Tx = length(Kx); % the sample size
Ty = length(Ky);
rx = pars.rx;
ry = pars.ry;
Tn = pars.np;
stId = pars.stId;

% boostrap parameters
if isfield(pars, 'T_BS')
    T_BS = pars.T_BS;
else
    T_BS = 10000;
end
Thresh = 1E-8;

% original HSIC
H = eye(Tn) - ones(Tn,Tn)/Tn; % for centering of the data in feature space

Hzx = eye(Tx)-ones(Tx)/Tx;
if ~isempty(Kzx)
    Kzx = Hzx*Kzx*Hzx;
    epsilon = 1e-5;
    Rx = epsilon*inv((Kzx+epsilon*eye(Tx)));
    Kx = Rx*Kx*Rx;
end
Hzy = eye(Ty)-ones(Ty)/Ty;
if ~isempty(Kzy)
    Kzy = Hzy*Kzy*Hzy;
    epsilon = 1e-5;
    Ry = epsilon*inv((Kzy+epsilon*eye(Ty)));
    Ky = Ry*Ky*Ry;
end

Kxl = Kx(1:Tn,1:Tn);
Kxlc = H*Kxl*H;
% Kyl = Ky(1:Tn,1:Tn);
% Kylc = H*Kyl*H;

% remove top eigs of Y (SNPs)
Hy = eye(Ty)-ones(Ty)/Ty;
KyC = Hy*Ky*Hy;
[uy1, dy1] = getEigen(KyC);

reg = 1e-5;
% remove top eigenvalues
KylP = Ky(1:Tn,:)*Hy*uy1(:,stId:rank(KyC))*inv(dy1(stId:rank(KyC),stId:rank(KyC))+reg*eye(1+rank(KyC)-stId))*uy1(:,stId:rank(KyC))'*Hy*Ky(:,1:Tn);

if strcmp(pars.uv,'v')
    KylPC = H*KylP*H;
    Sta = trace(Kxlc * KylPC)/Tn^2;
elseif strcmp(pars.uv,'u')
    onev = ones(Tn,1);
    Kxl_zero = Kxl;
    Kxl_zero(logical(eye(Tn))) = 0;
    Kyl_zero = KylP;
    KylPC = H*KylP*H;
    Kyl_zero(logical(eye(Tn))) = 0;
    Sta1 = trace(Kxl_zero*Kyl_zero)/Tn/(Tn-3);
    Sta2 = onev'*Kxl_zero*onev*onev'*Kyl_zero*onev/Tn/(Tn-1)/(Tn-2)/(Tn-3);
    Sta3 = onev'*Kxl_zero*Kyl_zero*onev/Tn/(Tn-2)/(Tn-3);
    Sta = Sta1 + Sta2 - 2*Sta3;
end
if pars.test == 1
    p_val0 = cal_pval(Sta, Kxlc, KylPC, Tn, T_BS, Tn, Tn, Thresh, stId, pars.uv);
else
    p_val0 = nan;
end
% Sta0 = trace(Kxl * Kyl)
% p_val00 = cal_pval(Sta0, Kxl, Kyl, Bootstrap, Approximate, Tn, T_BS, Tn, Tn, Thresh)

if Tn == Tx && Tn == Ty
    p_val = p_val0;
    p_valSemi = p_val0;
    StaSemi = Sta;
    return;
end
Hx = eye(Tx)-ones(Tx)/Tx;
KxC = Hx*Kx*Hx;
[ux1, dx1] = getEigen(KxC);

% original HSIC + modified null distr
if pars.test == 1
    p_val = cal_pval(Sta, KxC, KyC, Tn, T_BS, Tx, Ty, Thresh, stId, pars.uv);
else
    p_val = nan;
end

if rx == Tx && ry == Ty
    p_valSemi = p_val;
    StaSemi = Sta;
    return;
end

% modified HSIC + modified null distr
reg = 1e-5;
KxS = Kx(1:Tn,:)*Hx*ux1(:,1:rx)/(dx1(1:rx,1:rx)+reg*eye(rx))*ux1(:,1:rx)'*Hx*Kx(:,1:Tn);
KyS = Ky(1:Tn,:)*Hy*uy1(:,stId:ry)/(dy1(stId:ry,stId:ry)+reg*eye(ry-stId+1))*uy1(:,stId:ry)'*Hy*Ky(:,1:Tn);

if strcmp(pars.uv,'v')
    KxSC = H*KxS*H;
    KySC = H*KyS*H;
    StaSemi = trace(KxSC * KySC)/Tn^2;
elseif strcmp(pars.uv,'u')
    onev = ones(Tn,1);
    KxS_zero = KxS;
    KxS_zero(logical(eye(Tn))) = 0;
    KyS_zero = KyS;
    KyS_zero(logical(eye(Tn))) = 0;
    Sta1 = trace(KxS_zero *KyS_zero)/Tn/(Tn-3);
    Sta2 = onev'*KxS_zero*onev*onev'*KyS_zero*onev/Tn/(Tn-1)/(Tn-2)/(Tn-3);
    Sta3 = onev'*KxS_zero*KyS_zero*onev/Tn/(Tn-2)/(Tn-3);
    StaSemi = Sta1 + Sta2 - 2*Sta3;
end

if pars.test == 1
    p_valSemi = cal_pval(StaSemi, KxC, KyC, Tn, T_BS, rx, ry, Thresh, stId, pars.uv);
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
eig_Ky = eig_Ky(stId:end);
eig_prod = stack( (eig_Kx * ones(1,Num_eigY)) .* (ones(Num_eigX,1) * eig_Ky'));
II = find(eig_prod > max(eig_prod) * Thresh);
eig_prod = eig_prod(II); %%% new method
% sum((eig_Kx(1)/Tx).^2)
% sum((eig_Ky(1:21)/Ty).^2)

% use mixture of F distributions to generate the Null dstr
if strcmp(uv, 'u')
    f_rand1 = chi2rnd(1,length(eig_prod),T_BS) - 1;
elseif strcmp(uv, 'v')
    f_rand1 = chi2rnd(1,length(eig_prod),T_BS);
end
Null_dstr = eig_prod'*(1/Tn/Tx/Ty) * f_rand1; %%%%Problem
p_val = sum(Null_dstr>Sta)/T_BS;
end
