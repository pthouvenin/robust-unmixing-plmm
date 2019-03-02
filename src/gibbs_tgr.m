function xr = gibbs_tgr(mu,sigma,x,ar,br,r)
%Sample an element from the conditional law of the rth vector component,
% following a truncated Gaussian distribution[(ar,br) depending on x]
%
%-------------------------------------------------------------------------%
%%
% Author : Pierre-Antoine Thouvenin, 2016.
% [Code verification: 02/05/2016]
%-------------------------------------------------------------------------%
%%
% Inputs:
% > mu      avergage of the truncated Gaussian distribution [R,1]
% > sigma   covariance matrix [R,R]
% > x       initial vector [R, 1]
% > ar, br  bounds of the truncated interval [R,1]
% > r       component w.r.t. the element is sampled
%
%
% Outputs:
% < x       sampled vector [R, 1]
%-------------------------------------------------------------------------%
%%

R = size(x,1);
id = true(R,1);

id(r) = false;
sr = sigma(id,r); %sigma(r,id);
mu_r = mu(r) + sr'*(sigma(id,id)\(x(id) - mu(id))); % mu_r = mu(r) + (sr'/sigma(id,id))*(x(id) - mu(id));
sigma_r = sigma(r,r) - sr'*(sigma(id,id)\sr);       % sigma(r,r) - (sr'/sigma(id,id))*(sr); %sigma(r,r) - (sr*inv_sigma)*(sr'); %

% Generate a scalar according to a truncated Gaussian distribution (code by Vincent Mazet)
xr = rtnorm(ar,br,mu_r,sqrt(sigma_r));