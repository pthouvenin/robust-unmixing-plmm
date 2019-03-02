function x = gibbs_tg_splx(mu,sigma,x)
%Sample a vector x according to a multi-dimensional Gaussian distribution
%truncated to the unit simplex (using a Gibbs sampler).
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
%
%
% Outputs:
% < x       sampled vector [R, 1]
%-------------------------------------------------------------------------%
%%

R = size(x,1);
id = true(R,1);

for r = 1:R
    id(r) = false;
    sr = sigma(id,r); %sigma(r,id);
    mu_r = mu(r) + sr'*(sigma(id,id)\(x(id) - mu(id))); % mu_r = mu(r) + (sr'/sigma(id,id))*(x(id) - mu(id));
    sigma_r = sigma(r,r) - sr'*(sigma(id,id)\sr);       % sigma(r,r) - (sr'/sigma(id,id))*(sr); %sigma(r,r) - (sr*inv_sigma)*(sr'); %
    
    % Generate a scalar according to a truncated Gaussian distribution (code by Vincent Mazet)
    x(r) = rtnorm(0,1-sum(x(id)),mu_r,sqrt(sigma_r));
    id(r) = true;
end

end

