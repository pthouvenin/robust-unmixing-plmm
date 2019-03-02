function M = sample_M(Y,M,A,dM,X,Sig,xi)
%Sample the endmember matrix M s.t. M+dM >=0 (see (47)-(50)).
%
%-------------------------------------------------------------------------%
%%
% Author : Pierre-Antoine Thouvenin, 2016.
% [Code verification: 06/06/2016]
%-------------------------------------------------------------------------%
%%
% Inputs:
% > Y       pixels (hyperspectral data cube reshaped in matrix form) [L,N,T]
% > M       endmember matrix [L,R]
% > A       abundance matrix [R,N,T]
% > dM      variability matrices [L,R,T]
% > X       outlier terms [L,N,T] 
% > Sig     noise variances \sig_t^2 [T,1]
% > xi      hyperparameter: prior variance of the endmember matrix [1]
%
% Outputs:
% < M       endmember matrix [L,R]
%-------------------------------------------------------------------------%
%%

[L,R,T] = size(dM);
mu0 = sum(bsxfun(@rdivide,mtimesx(Y - mtimesx(dM,A) - X,A,'t'),reshape(Sig,[1,1,T])),3);
S = 1./sqrt( sum(sum(bsxfun(@rdivide,A.^2,reshape(Sig,[1,1,T])),2),3) + 1/xi); % [R,1]

temp = zeros(R-1,R);
id = true(R,1);
for r = 1:R
    id(r) = false;
    temp(:,r) = sum(sum(bsxfun(@times,A(id,:,:),bsxfun(@rdivide,A(r,:,:),reshape(Sig,[1,1,T]))),2),3);
    id(r) = true;
end

% Constraint
a = min(dM,[],3);
a = bsxfun(@max,0,-a);

% b = max(dM,[],3);
% b = bsxfun(@min,1,1-b);

for l = 1:L
    for r = 1:R
        id(r) = false;          
        mu_lr = (S(r)^2)*( mu0(l,r) - M(l,id)*temp(:,r));
        % M(l,r) = rtnorm(a(l,r),b(l,r),mu_lr,S(r));
        M(l,r) = rtnorm(a(l,r),Inf,mu_lr,S(r));
        id(r) = true;
    end
end

end
