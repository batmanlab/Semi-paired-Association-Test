function [xa, ya] = dataset_simul1(pars)
% Simulate low-dimesional data from random effect model, to test type I error

h2 = pars.h2; % heritability
sigma_nx = pars.sigma_nx; % noise variance from low dimension genotype to high dimension
sigma_ny = pars.sigma_ny; % noise variance from low dimension phenotype to high dimension
sigma = pars.sigma; 
dzx = pars.xzdim; % genotype hidden dimension
dzy = pars.yzdim; % phenotype hidden dimension
dx = pars.xdim; % genotype dimension
dy = pars.ydim; % phenotype dimension
N = pars.N; % sample size
T = pars.T; % number of repetitions
LMM = pars.LMM;
data_dir = pars.data_dir; % data folder

if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

% set path
data_subdir = sprintf('h2_%.2f_sig%.2f_signx%.2f_signy%.2f_dzx%d_dzy%d_dx%d_dy%d_N%d_LMM%d.mat', h2, sigma, sigma_nx, sigma_ny, dzx, dzy, dx, dy, N, LMM);
data_path = fullfile(data_dir,data_subdir);

if exist(data_path, 'file')
	load(data_path);
else
	ya = zeros(N,dy,T);
	rng(0);
    if LMM==0
        xa = zeros(N,dx,T);
		for t=1:T
			% generate X
			zx = randn(N,dzx);
			A = randn(dx,dzx);
			[U,D] = eig(A*A');
			U = U(:,1:dzx);
			% U = A./repmat(sum(A.^2,2),1,dzx);
			x = zx*U' + sigma_nx*randn(N,dx);
			
			% generate Y
			y_h = zeros(N,dzy);
			for i = 1:dzy
			    y_h(:,i) = mvnrnd(zeros(1,N), zx*U'*U*zx'*sigma^2*h2/dx + sigma^2*(1-h2)*eye(N));
			end
			B = randn(dy,dzy);
			[V,D1] = eig(B*B');
			V = V(:,1:dzy);
			y = y_h * V' + sigma_ny*randn(N,dy);
			
			% calculate gram matrix
% 			Kxa(:,:,t) = x*x';
			ya(:,:,t) = y;
            xa(:,:,t) = x;
		end
    else
		% generate X
		zx = randn(N,dzx);
		A = randn(dx,dzx);
		[U,D] = eig(A*A');
		U = U(:,1:dzx);
		% U = A./repmat(sum(A.^2,2),1,dzx);
		x = zx*U' + sigma_nx*randn(N,dx);	
        xa = x;
		for t=1:T
            % generate Y
			y_h = zeros(N,dzy);
			for i = 1:dzy
			    y_h(:,i) = mvnrnd(zeros(1,N), x*x'*sigma^2*h2/dx + sigma^2*(1-h2)*eye(N));
			end
			B = randn(dy,dzy);
			[V,D1] = eig(B*B');
			V = V(:,1:dzy);
			ya(:,:,t) = y_h * V' + sigma_ny*randn(N,dy);
		end
    end
    save(data_path, 'ya', 'xa', '-v7.3');
end
