function X = sample_X(mu,Lambda,Z)
%Sample the outlier terms x_{n,t} using (55).
%
%-------------------------------------------------------------------------%
%%
% Author : Pierre-Antoine Thouvenin, 2016.
% [Code verification: 02/05/2016]
%-------------------------------------------------------------------------%
%%
% Inputs:
% > mu      mean of the conditional distribution of the outliers [L,N,T]
% > Lambda  variance for the conditional disctribution of the outliers [T,1]
% > Z       binary outlier labels [N,T]
%
% Outputs:
% < X       outlier terms [L,N,T] 
%-------------------------------------------------------------------------%
%%
% Version with a single variance per spectral band

L = size(mu,1);
[N,T] = size(Z);
X = zeros(L,N,T);

for t = 1:T
    if any(Z(:,t))
        ind = find(Z(:,t));
        for l = 1:L
            for n = 1:numel(ind)
                X(l,ind(n),t) = rtnorm(0,Inf,mu(l,ind(n),t),sqrt(Lambda(t))); 
            end
        end
    end
end

end
