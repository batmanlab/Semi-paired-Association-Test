function K = gaussrbf( X, Xt, varargin )
%rbf_kernel - the Gaussian RBF kernel, <psi, psi>

if isempty(varargin)
    sig = 1;  % the default value of sigma
else
    sig = varargin{1};
end

D = squaredist(X, Xt);
K = exp(-0.5.*D./(sig^2));
