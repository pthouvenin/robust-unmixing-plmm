function M = sample_M_subspace(Y,M,A,dM,X,Sig,xi,U,P,Up,Um,idp,idm,y_bar)
%Sample the endmember matrix M in the dimensionality-reduced space (PCA) 
%under the constraint M+dM >=0 (see (67) and (68)). 
%
%-------------------------------------------------------------------------%
%%
% Author : Pierre-Antoine Thouvenin, 2016.
% [Code verification: 10/06/2016]
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
% > U       back projector onto the orignal space [L,K]
% > P       projector onto the PCA-subspace [K,L]
% > Up      binary matrix: positive elements in U [L,K]
% > Um      binary matrix: negative elements in U [L,K] 
% > idp     any(Up,1);
% > idm     any(Um,1);
% > y_bar   empirical average of the data over space and time [L,1]
%
% Outputs:
% < M       endmember matrix [L,R]
%-------------------------------------------------------------------------%
%%

% Auxiliary variables
[R,~,T] = size(A);
K = size(U,2); % K = R-1 (dimension of the subspace)
idk = true(1,K);
idr = true(1,R);
mu_tmp = bsxfun(@rdivide,Y - X - mtimesx(bsxfun(@plus,dM,y_bar),A),reshape(Sig,[1,1,T])); % [L,N,T]
lambda_tmp = sum(bsxfun(@rdivide,sum(A.^2,2),reshape(Sig,[1,1,T])),3); % [R,1]
B = bsxfun(@min,0,min(dM,[],3)); % [07/09/2016] (bounds)

% Projection onto the PCA space (M -> Tm)
Tm = P*bsxfun(@minus,M,y_bar);

% Sampling
for r = 1:R
    
    idr(r) = false;    
    
    % Mean and standard deviation
    lambda_r = (lambda_tmp(r) + 1/xi)*eye(K); % simplification U'*U = I
    mu_r = lambda_r\( (U')*(sum(mtimesx(mu_tmp,A(r,:,:),'t'),3) - y_bar/xi) - sum(mtimesx(bsxfun(@rdivide,mtimesx(Tm(:,idr),A(idr,:,:)),reshape(Sig,[1,1,T])),A(r,:,:),'t'),3) );
    lambda_r = inv(lambda_r);

    % Interval bounds
    for k = 1:K        
        %  !! The bounds evolve as t_kr is updated !!
        %  (truncated gaussian to be included inside the loop!)
        tm = -Inf;
        tp = Inf;
        idk(k) = false;       
        v = - bsxfun(@rdivide, y_bar + U(:,idk)*Tm(idk,r) + B(:,r), U(:,k)); % [correction : 07/09/2016] 
        idk(k) = true;
        
        % t-
        if idp(k)
           tm = max(v(Up(:,k)));
        end
        % t+
        if idm(k)
           tp = min(v(Um(:,k)));
        end
     
        % Sample T(k,r) (truncated Gaussian distribution)
        Tm(k,r) = gibbs_tgr(mu_r,lambda_r,Tm(:,r),tm,tp,k);
    end  
        
    idr(r) = true;
end

% Back-projection onto the original space (Tm -> M)
M = bsxfun(@plus,U*Tm,y_bar);

end

