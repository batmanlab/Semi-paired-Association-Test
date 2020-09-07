function [d,u] = sortEigen(d,u)
% sort eigenvalues in non-increasing order
if ~issorted(diag(d),'descend')
    [d,I] = sort(diag(d),'descend');
    d = diag(d);
    u = u(:,I);
end