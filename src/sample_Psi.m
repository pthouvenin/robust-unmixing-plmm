function Psi = sample_Psi(dM,a_psi,b_psi,T)
%Sample the variability variances \psi_{\ell,r}^2 using (61).
%
%-------------------------------------------------------------------------%
%%
% Author : Pierre-Antoine Thouvenin, 2016.
% [Code verification: 22/06/2016]
%-------------------------------------------------------------------------%
%%
% Inputs:
% > dM      variability matrices [L,R,T]
% > a_psi   hyperparameter of the prior on the outlier variances (31)
% > b_psi   hyperparameter of the prior on the outlier variances (31)
% > T       number of time instants
%
% Outputs:
% < Psi     variability variances [L,R]
%-------------------------------------------------------------------------%
%%

Psi = 1./gamrnd(a_psi + (T-1)/2, 1./bsxfun(@plus,b_psi,0.5*sum(diff(dM,1,3).^2,3))); % [L,R]
% warning : beware MATLAB definition of the gamma distribution 

end
