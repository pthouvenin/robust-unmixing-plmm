function s = sample_s(a_s,b_s,X,Z)
%Sample the outlier variances \s_{t}^2 using (59).
%
%-------------------------------------------------------------------------%
%%
% Author : Pierre-Antoine Thouvenin, 2016.
% [Code verification: 03/05/2016]
%-------------------------------------------------------------------------%
%%
% Inputs:
% > a_s     hyperparameter of the prior on the outlier variances (31)
% > b_s     hyperparameter of the prior on the outlier variances (31)
% > X       outlier term [L,N,T] 
% > Z       binary outlier labels [N,T]
%
% Outputs:
% < s       outlier variances [1,T]
%-------------------------------------------------------------------------%
%%

[L,~,T] = size(X);
s = 1./gamrnd(L*sum(Z,1)/2 + a_s, 1./reshape((b_s + 0.5*sum(sum(X.^2,1),2)),[1,T])); % [1,T]
% warning : beware MATLAB definition of the gamma distribution

end
