function Sig = sample_Sig(Y,M,A,dM,X,alpha_n,beta_n)
%Sample the noise variances \sigma_{t}^2 using (60).
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
% > X       outlier term [L,N,T] 
% > alpha_n hyperparameter of the prior on the noise variances (29)
% > beta_n  hyperparameter of the prior on the noise variances (29)
%
% Outputs:
% < Sig     noise variances [1,T]
%-------------------------------------------------------------------------%
%%

[L,N,T] = size(Y);
temp = reshape(sum(sum((Y - mtimesx(bsxfun(@plus,M,dM),A) - X).^2,2),1),[1,T]); % size [1,T]
Sig = 1./gamrnd((N*L)/2 + alpha_n, 1./(temp/2 + beta_n)); % warning : MATLAB definition of gamma distribution 

end
