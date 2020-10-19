import numpy as np
from numpy.linalg import inv, matrix_rank, multi_dot
from scipy.sparse.linalg import eigs
from scipy.linalg import eigh as eig

def HSIC_Test(Kx, Ky, Kzx, Kzy, pars):
    # Hilbert Schmidt Independence Criterion (HSIC) test using semi-paired data
    # Inputs: 
    #       Kx - kernel matrix on x (NxN, the first np x np block contains paired data)
    #       Ky - kernel matrix on y (MxM, the first np x np block contains paired data)
    #       Kzx - kernel matrix on covariate z (for x, empty if not exist)
    #       Kzy - kernel matrix on covariate z (for y, empty if not exist)
    #       pars - hyperparameters
    #           .np - number of paired data
    #           .rx - reduced dimension of x, default pars.np
    #           .ry - reduced dimension of y, default pars.np
    #           .stId - remove the top (stId) eigenvalues (for genotype data), defaut 0
    #           .T_BS - bootstrap sample size, default 10000
    #           .uv - 'u' or 'v' statistics, default 'u'
    #           .test - compute p-values if true, default true
    # Outputs:
    #       p_val0 - p value of the original HSIC using only paired data
    #       p_val - p value of our Semi-paired test (SAT), only improve null distribution
    #       p_valSemi - p value of our SAT, improve both test statistics and null distribution
    #       Sta - test statistic of the original HSIC using only paired data
    #       StaSemi - test statistic of our SAT

    Tx = len(Kx) # The sample size
    Ty = len(Ky)
    rx = pars['rx']
    ry = pars['ry']
    Tn = pars['np']
    stId = pars['stId']

    # Boostrap parameters
    T_BS = pars.get('T_BS', 10000)
    Thresh = 1e-8

    # Original HSIC
    H = np.eye(Tn) - np.ones((Tn, Tn)) / Tn
    Hzx = np.eye(Tx) - np.ones((Tx, Tx)) / Tx

    if Kzx.size != 0:
        Kzx = multi_dot([Hzx, Kzx, Hzx])
        epsilon = 1e-5
        Rx = epsilon*inv((Kzx+epsilon*np.eye(Tx)))
        Kx = multi_dot([Rx, Kx, Rx])

    if Kzy.size != 0:
        Kzy = multi_dot([Hzy, Kzy, Hzy])
        epsilon = 1e-5
        Ry = epsilon*inv((Kzy+epsilon*np.eye(Ty)))
        Ky = multi_dot([Ry, Ky, Ry])

    Kxl = Kx[0:Tn,0:Tn].copy()
    Kxlc = multi_dot([H, Kxl, H])

    # remove top eigs of Y (SNPs)
    Hy = np.eye(Ty)-np.ones((Ty, Ty))/Ty
    KyC = multi_dot([Hy, Ky, Hy])
    uy1, dy1 = getEigen(KyC)

    reg = 1e-5

    # remove top eigenvalues
    KylP = multi_dot([\
             Ky[:Tn,:],\
             Hy,\
             uy1[:,stId:matrix_rank(KyC)],\
             inv(dy1[stId:matrix_rank(KyC),stId:matrix_rank(KyC)]+\
               reg*np.eye(matrix_rank(KyC)-stId)),\
             np.transpose(uy1[:,stId:matrix_rank(KyC)]),\
             Hy,\
             Ky[:,:Tn]])

    if pars['uv'] == 'v':
        KylPC = multi_dot([H, KylP, H])
        Sta = np.trace(np.dot(Kxlc, KylPC))/Tn^2
    elif pars['uv'] == 'u':
        onev = np.ones((Tn,1))
        Kxl_zero = Kxl
        Kxl_zero[np.eye(Tn).astype(np.bool)] = 0
        Kyl_zero = KylP
        KylPC = multi_dot([H, KylP, H])
        Kyl_zero[np.eye(Tn).astype(np.bool)] = 0
        Sta1 = np.trace(np.dot(Kxl_zero, Kyl_zero))/Tn/(Tn-3)
        Sta2 = multi_dot([np.transpose(onev), Kxl_zero, onev, np.transpose(onev), Kyl_zero, onev])/Tn/(Tn-1)/(Tn-2)/(Tn-3)
        Sta3 = multi_dot([np.transpose(onev), Kxl_zero, Kyl_zero, onev])/Tn/(Tn-2)/(Tn-3)
        Sta = Sta1 + Sta2 - 2*Sta3

    if pars['test'] == 1:
        p_val0 = cal_pval(Sta, Kxlc, KylPC, Tn, T_BS, Tn, Tn, Thresh, stId, pars['uv'])
    else:
        p_val0 = None

    if Tn == Tx and Tn == Ty:
        p_val = p_val0
        p_valSemi = p_val0
        StaSemi = Sta
        return p_val0, p_val, p_valSemi, Sta, StaSemi

    Hx = np.eye(Tx)-np.ones((Tx, Tx)) / Tx
    KxC = multi_dot([Hx, Kx, Hx])
    ux1, dx1 = getEigen(KxC)

    # original HSIC + modified null distr
    if pars['test'] == 1:
        p_val = cal_pval(Sta, KxC, KyC, Tn, T_BS, Tx, Ty, Thresh, stId, pars['uv'])
    else:
        p_val = None

    if rx == Tx and ry == Ty:
        p_valSemi = p_val
        StaSemi = Sta
        return p_val0, p_val, p_valSemi, Sta, StaSemi

    # modified HSIC + modified null distr
    reg = 1e-5
    KxS = multi_dot([Kx[:Tn,:],
                     Hx,\
                     right_divide(ux1[:,:rx], (dx1[:rx,:rx]+reg*np.eye(rx))),\
                     np.transpose(ux1[:,:rx]),\
                     Hx,\
                     Kx[:,:Tn]])
    KyS = multi_dot([Ky[:Tn,:],
                     Hy,
                     right_divide(uy1[:,stId:ry], (dy1[stId:ry,stId:ry]+reg*np.eye(ry-stId))),
                     uy1[:,stId:ry].transpose(),
                     Hy,
                     Ky[:,:Tn]])
    #print(KxS)
    #print(KyS)

    if pars['uv'] == 'v':
        KxSC = multi_dot([H, KxS, H])
        KySC = multi_dot([H, KyS, H])
        StaSemi = np.trace(np.dot(KxSC, KySC)) / (Tn^2)
    elif pars['uv'] == 'u':
        onev = np.ones((Tn,1))
        KxS_zero = KxS
        KxS_zero[np.eye(Tn).astype(np.bool)] = 0
        KyS_zero = KyS
        KyS_zero[np.eye(Tn).astype(np.bool)] = 0
        Sta1 = np.trace(np.dot(KxS_zero, KyS_zero)) / Tn / (Tn-3)
        Sta2 = multi_dot([np.transpose(onev),
                          KxS_zero,
                          onev,
                          np.transpose(onev),
                          KyS_zero,
                          onev]) /Tn/(Tn-1)/(Tn-2)/(Tn-3)
        Sta3 = multi_dot([np.transpose(onev),
                          KxS_zero,
                          KyS_zero,
                          onev]) /Tn/(Tn-2)/(Tn-3)
        #print(Sta1)
        #print(Sta2)
        #print(Sta3)
        StaSemi = Sta1 + Sta2 - 2*Sta3

    if pars['test'] == 1:
        p_valSemi = cal_pval(StaSemi, KxC, KyC, Tn, T_BS, rx, ry, Thresh, stId, pars['uv'])
    else:
        p_valSemi = None
    return p_val0, p_val, p_valSemi, Sta, StaSemi

def right_divide(A, B):
    # Get x in xB = A
    return np.linalg.solve(B.conj().T, A.conj().T).conj().T

def getEigen(M):
    # compute and sort eigenvalues
    M = 1 / 2 * (M + np.transpose(M))
    eigenValues, eigenVectors = eig(M)
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    return eigenVectors, np.diag(eigenValues)

def eigdec(x, N):
    # EIGDEC Sorted eigendecomposition
    if N / x.shape[1] > 0.04:
        temp_evals, temp_evec = eig(x)
    else:
        temp_evals, temp_evec = eigs(x, k=N, which='LM')
    idx = temp_evals.argsort()[::-1]
    temp_evals = temp_evals[idx]
    temp_evec = temp_evec[:,idx]
    return temp_evals[:N], temp_evec[:,:N]

def cal_pval(Sta, Kx, Ky, Tn, T_BS, Num_eigX, Num_eigY, Thresh, stId, uv):
    # calculate the eigenvalues that will be used later
    # Due to numerical issues, Kx and Ky may not be symmetric:
    Tx = len(Kx)
    Ty = len(Ky)
    eig_Kx, eivx = eigdec((Kx+np.transpose(Kx))/2, Num_eigX)
    eig_Ky, eivy = eigdec((Ky+np.transpose(Ky))/2, Num_eigY)
    Num_eigY = Num_eigY-stId
    # calculate Cri...
    # first calculate the product of the eigenvalues
    eig_Ky = eig_Ky[stId:]
    eig_prod = np.multiply( np.dot(eig_Kx.reshape((-1,1)), np.ones((1, Num_eigY))),\
                            np.dot(np.ones((Num_eigX,1)), eig_Ky.reshape(1,-1))).reshape((-1,))
    #print(eig_Kx)
    #print(eig_prod.shape)
    II = (eig_prod > np.max(eig_prod) * Thresh).astype(np.bool)
    eig_prod = eig_prod[II] # new method
    #print(eig_prod.shape)
    # sum((eig_Kx(1)/Tx).^2)
    # sum((eig_Ky(1:21)/Ty).^2)

    # use mixture of F distributions to generate the Null dstr
    if uv == 'u':
        f_rand1 = np.random.chisquare(1, (len(eig_prod), T_BS)) - 1
    elif uv == 'v':
        f_rand1 = np.random.chisquare(1, (len(eig_prod), T_BS))

    Null_dstr = np.dot(np.transpose(eig_prod)*(1/Tn/Tx/Ty), f_rand1) # Problem
    p_val = (Null_dstr>Sta).astype(np.int).sum()/T_BS
    return p_val
