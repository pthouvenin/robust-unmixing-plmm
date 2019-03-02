function [L,S] = pcp(M,delta)
% Principal component pursuit solved by ADMM, corresponding to the
% following optimization problem
%
% min_{L,S} |L|_{*} + lambda|S|_{1} s.t. L + S = M.
%-------------------------------------------------------------------------%
% Input:
% < M      data matrix (model : M = L + S)
% < delta  tolerance of the convergence criterion
%
% Output
% > L      dense compoent of M
% > S      sparse compoent of M
%-------------------------------------------------------------------------%
% Reference
% [1] Candès, Li, Ma, Wright (2011) Robust Principal Component Analysis ?
%-------------------------------------------------------------------------%
%%
[n1,n2] = size(M);
lambda = 1/sqrt(max([n1,n2]));
norm_fro = @(x) sqrt(sum(x(:).^2)); % Frobenius norm
shrink = @(x,l) sign(x).*max(bsxfun(@minus,abs(x),l),0);

% Initialization
S = zeros(n1,n2);
Y = zeros(n1,n2); % Lagange mutliplier
mu = n1*n2/(4*sum(abs(M(:))));
inv_mu = 1/mu;
err = 1;
normF_M = norm_fro(M);

% Algorithm
while err    
    % Update L
    L = M - S + Y*inv_mu;
    [U,D,V] = svd(L);
    L = U*shrink(D, inv_mu)*(V');     
    % Update S
    S = shrink(M - L + inv_mu*Y, lambda*inv_mu);
    % Update Y
    Y = Y + mu*(M - L - S); 
    % Test convergence criterion
    err = (norm_fro(M - L - S) > delta*normF_M);
end

end

