function [u, d] = getEigen(M)
% compute and sort eigenvalues
% try 
%    [u,d] = svd(M);
% catch ErrorInfo 
%     disp(ErrorInfo);
%     disp(ErrorInfo.identifier);
%     disp(ErrorInfo.message);
%     disp(ErrorInfo.stack);
%     disp(ErrorInfo.cause);
% end
M = 1/2*(M+M');
[u,d] = eig(M);
[d,I] = sort(diag(d),'descend');
d = diag(d);
u = u(:,I);
end
