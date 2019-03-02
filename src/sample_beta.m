function [Beta,id_accept] = sample_beta(Beta,Cn,Z,H,W,s_beta)
%Sample the granularity parameters \beta_t using a Metropolis-Hastings 
%procedure (closed-form expression for the normalizing constant of the 
%Ising field in the case of square images), using Section IV.I and (62).
%
%-------------------------------------------------------------------------%
%%
% Author : Pierre-Antoine Thouvenin, 2016.
% [Code verification: 21/06/2016]
%-------------------------------------------------------------------------%
%%
% Inputs:
% > Beta    granularity parameter in the Potts model [T,1]
% > Cn      auxiliary parameter used in the sampling of \beta_t [N,1]
%           (see Section IV.I)
% > Z       current labels [N,T]
% > H       heigth of the image
% > W       width of the image; (N = W*H)
% > s_beta  variances for the proposals used in the MH step [T,1]
%
% Outputs:
% < Z       outlier labels [N,T]
% < id_accept  boolean vector indicating the elements updated in Beta 
%              [T, 1] (used to adjust s_beta)
%-------------------------------------------------------------------------%
%%
% Sources: 
% [1] J. F. Giovannelli, "Ising field parameter estimation from incomplete 
%     and noisy data," 2011 18th IEEE International Conference on Image 
%     Processing, Brussels, 2011, pp. 1853-1856.
% [2] J. F. Giovannelli, "Estimation of the Ising field parameter thanks to
%     the exact partition function," 2010 IEEE International Conference on 
%     Image Processing, Hong Kong, 2010, pp. 1441-1444.
%-------------------------------------------------------------------------%
%% 
% REMARKS:  
% - Exact formula for square lattices, used as an approximation in the 
% general case (approximation valid for sufficiently/reasonably large 
% images);
% - Beware definition of the neighbourhood adopted in [1]
% - Prior: beta ~ U[0,2]; 
% - Proposal: Gaussian random walk beta_star = beta(t) + e_t, e_t ~ N(0,s_{\beta}^2),
% s_{\beta}^2 adjusted to obtain an acceptance rate between 40% and 60%.
%-------------------------------------------------------------------------%
%%

N = numel(Cn); % Cn : [N,1]
T = length(Beta);
beta_p = Beta + sqrt(s_beta).*randn(1,T);
id = bsxfun(@and,beta_p >= 0, beta_p <= 2);

% Vertical/horizontal neighbourhood
Z_im = reshape(Z,[H,W,T]);
Dv = (diff(Z_im,1,1) == 0);
Dh = (diff(Z_im,1,2) == 0);
% 4 neighbourhood (Left, Right, Up, Down)
D = sum(reshape(vertcat(zeros(1,W,T),Dv) + vertcat(Dv,zeros(1,W,T)) + horzcat(zeros(H,1,T),Dh) + horzcat(Dh,zeros(H,1,T)),[N,T]),1); % [1,T]

% Log acceptance ratio
log_rho = (beta_p - Beta).*(D/2 - N) + 0.5*( N*log(sinh(Beta)./sinh(beta_p)) ...
          + sum( acosh(bsxfun(@minus,((cosh(Beta)).^2)./sinh(Beta),Cn)) - ...
          acosh(bsxfun(@minus,((cosh(beta_p)).^2)./sinh(beta_p),Cn)), 1)); % MRF contribution
          % warning : D/2, related to the definition adopted for the neighbourhood
% log_rho(~id) = NaN;

% Accept/reject procedure
id_accept = (log(rand(1,T)) < log_rho);
id_accept = bsxfun(@and,id_accept,id);
Beta(id_accept) = beta_p(id_accept);

end
