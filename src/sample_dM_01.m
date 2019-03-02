function dM = sample_dM_01(Y,M,A,dM,X,Sig,Psi,nu_dM)
%Sample the variability matrix dM s.t. 0 <= M+dM <= 1 (see (52)-(54)).
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
% > Psi     variability prior variances \psi_{\ell,r}^2 [L,R] 
%           (time smoothness prior)
% > nu_dM   hyperparameter: prior variance of the endmember matrix 
%           dM(:,:,1) [1]
%
% Outputs:
% < M       endmember matrix [L,R]
%-------------------------------------------------------------------------%
%%
% Constraint M+dM \in [0, 1]^{L, R}

[L,R,T] = size(dM);
mu0 = bsxfun(@rdivide,mtimesx(Y - mtimesx(M,A) - X,A,'t'),reshape(Sig,[1,1,T])); % [L,R,T]

S1 = bsxfun(@rdivide,sum(A.^2,2),reshape(Sig,[1,1,T])); % [R,1,T]
S = zeros(R,L,T);

S(:,:,2:T-1) = bsxfun(@plus,S1(:,1,2:T-1),2./Psi');
S(:,:,[1,T]) = bsxfun(@plus,S1(:,1,[1,T]),1./Psi');
S(:,:,1) = bsxfun(@plus,S(:,:,1),1/nu_dM);
S = 1./sqrt(S);

temp = zeros(R-1,R,T);
id = true(R,1);
for r = 1:R
    id(r) = false;
    temp(:,r,:) = sum(bsxfun(@times,A(id,:,:),bsxfun(@rdivide,A(r,:,:),reshape(Sig,[1,1,T]))),2); % [R-1,1,T]
    id(r) = true;
end

% t = 1
for l = 1:L
    for r = 1:R
        id(r) = false;          
        mu_lr = (S(r,l,1)^2)*( mu0(l,r,1) - dM(l,id,1)*temp(:,r,1) + dM(l,r,2)/Psi(l,r));
        dM(l,r,1) = rtnorm(-M(l,r),1-M(l,r),mu_lr,S(r,l,1));
        id(r) = true;
    end
end

% t > 1 & t < T
for t = 2:T-1
    for l = 1:L
        for r = 1:R
            id(r) = false;          
            mu_lr = (S(r,l,t)^2)*( mu0(l,r,t) - dM(l,id,t)*temp(:,r,t) + (dM(l,r,t+1) + dM(l,r,t-1))/Psi(l,r));
            dM(l,r,t) = rtnorm(-M(l,r),1-M(l,r),mu_lr,S(r,l,t));
            id(r) = true;
        end
    end
end

% t = T
for l = 1:L
    for r = 1:R
        id(r) = false;          
        mu_lr = (S(r,l,T)^2)*( mu0(l,r,T) - dM(l,id,T)*temp(:,r,T) + dM(l,r,T-1)/Psi(l,r));
        dM(l,r,T) = rtnorm(-M(l,r),1-M(l,r),mu_lr,S(r,l,T));
        id(r) = true;
    end
end

end
