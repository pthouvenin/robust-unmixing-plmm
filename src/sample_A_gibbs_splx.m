function A = sample_A_gibbs_splx(Y,M,A,dM,Z,X,Eps,Sig)
%Sample the abudance coefficients a_{n,t} using (42)-(46) with a Gibbs 
%sampler.
%
%-------------------------------------------------------------------------%
%%
% Author : Pierre-Antoine Thouvenin, 2016.
% [Code verification: 02/05/2016]
%-------------------------------------------------------------------------%
%%
% Inputs:
% > Y       pixels (hyperspectral data cube reshaped in matrix form) [L,N,T]
% > M       endmember matrix [L,R]
% > A       abundance matrix [R,N,T]
% > dM      variability matrices [L,R,T]
% > Z       current labels [N,T]
% > X       outlier terms [L,N,T] 
% > Eps     abundance prior variance [1,N] (see (8)) (temporal smoothing)
% > Sig     noise variances \sig_t^2 [T,1]
%
% Outputs:
% < A       abundance matrix [R,N,T]
%-------------------------------------------------------------------------%
%%
% Version with the sum-to-one constraint

% Preliminary computations
[R,N,T] = size(A);
MdM = bsxfun(@plus,M,dM);

% Auxiliary variables
% Lambda0 = mtimesx(MdM,'t',bsxfun(@rdivide,MdM,reshape(Sig,[1,1,T])));
% mu0 = mtimesx(MdM,'t',bsxfun(@rdivide,Y - X,reshape(Sig,[1,1,T])));
id = true(1,R);
r = ceil(R*rand(1));
id(r) = false;

% Abundance sampling (scalar)
for t = 1:T
    for n = 1:N
        if (Z(n,t) == 0)
            % Auxiliary variables (initialization)
            alpha_nt = 0;
            mu_tmp = zeros(R-1,1);
            
            % Reset id + randomly determine the position of element which is not
            % sampled
            id(r) = true;
            r = ceil(R*rand(1));
            id(r) = false;
            
            % Find \delta_n1 and \delta_n2
            if t > 1
                delta = find(~Z(n,1:t-1),1,'last');
                if ~isempty(delta)
                    alpha_nt = alpha_nt + 1; 
                    mu_tmp = mu_tmp + bsxfun(@plus,A(id,n,delta),1-A(r,n,delta))/Eps(n);
                end
            end 
            if t < T
                delta = find(~Z(n,t+1:end),1,'first');
                if ~isempty(delta)
                    delta = t + delta;
                    alpha_nt = alpha_nt + 1; 
                    mu_tmp = mu_tmp + bsxfun(@plus,A(id,n,delta),1-A(r,n,delta))/Eps(n);
                end
            end
            
            % Compute mean and covariance matrix
            temp = bsxfun(@minus,MdM(:,id,t),MdM(:,r,t));
            lambda_tmp = temp'*temp/Sig(t) + alpha_nt*(eye(R-1) + ones(R-1))/Eps(n);
            mu_tmp = lambda_tmp\(temp'*(Y(:,n,t) - MdM(:,r,t))/Sig(t) + mu_tmp);
            lambda_tmp = inv(lambda_tmp);
            
            % Sample A (Gaussian truncated on the (R-1)-simplex)
            A(id,n,t) = gibbs_tg_splx(mu_tmp,lambda_tmp,A(id,n,t)); 
            A(r,n,t) = 1 - sum(A(id,n,t)); % sum-to-one constraint
        else
            % Sample A (Gaussian truncated on the intersection of the l1 ball (0,1) and of the positive orthant)
            % Compute mean and covariance matrix
            temp = bsxfun(@minus,MdM(:,id,t),MdM(:,r,t));
            lambda_tmp = temp'*temp/Sig(t);
            mu_tmp = lambda_tmp\(temp'*(Y(:,n,t) - MdM(:,r,t) - X(:,n,t))/Sig(t));
            lambda_tmp = inv(lambda_tmp);
            
            % Sample A (Gaussian truncated on the (R-1)-simplex)
            A(id,n,t) = gibbs_tg_splx(mu_tmp,lambda_tmp,A(id,n,t)); 
            A(r,n,t) = 1 - sum(A(id,n,t)); % sum-to-one constraint
        end
    end    
end

end
