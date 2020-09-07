% standard deviation: sqrt(sn_est)
%h2 = 0.5, hf2 = 0.1; sn_est = 0.36, 0.45, 0.63, 0.8944 for
%P=300,200,100,50
%h2 = 0.5, hf2 = 0.2; sn_est = 0.22, 0.27, 0.39, 0.5477 for
%P=300,200,100,50

% variance: (sn_est), use this
%h2 = 0.5, hf2 = 0.1; sn_est = 0.133, 0.2, 0.4, 0.8, 2 for
%P=300,200,100,50,20
%h2 = 0.5, hf2 = 0.2; sn_est = 0.05, 0.075, 0.15, 0.3, 0.75 for
%P=300,200,100,50,20

h2 = 0.5;
hf2 = 0.2;
p = 10;
P = 50;
A = randn(P,p);
[U,D] = eig(A*A');
U = U(:,1:p);
Sigma_g2 = h2*eye(p);
sn_est = (h2*p/hf2-p)/P
sqrt(sn_est)
Sigma_e2 = U*(1-h2)*eye(p)*U' + sn_est*eye(P);
hf = trace(Sigma_g2)/(trace(Sigma_g2)+trace(Sigma_e2));
hf1 = h2*p/(p+sn_est*P);