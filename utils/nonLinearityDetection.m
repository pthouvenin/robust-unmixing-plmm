function [map,sigma2] = nonLinearityDetection(Y,M,pfa)
%Non-lineariy detector (adapted from [1]).
% H0: y = Ma
% H1: y = Ma + r
%-------------------------------------------------------------------------%
%%
% Input
% > Y      [L,N]    hyperspectral image
% > pfa             error of the first kind (false alarm probability)
% > R               endmember number
% > sigma2          noise variance estimate
%
% Ouput
% < map    [1,N]    non-linearity detection map (1: H1 model, 0: H0)
%-------------------------------------------------------------------------%
%%
% Reference:
%
% [1]  Y. Altmann, N. Dobigeon, J.-Y. Tourneret and J. C. M. Bermudez, 
% "A robust test for nonlinear mixture detection in hyperspectral images," 
% in Proc. IEEE Int. Conf. Acoust., Speech, and Signal Processing (ICASSP),
% Vancouver, Canada, June 2013, pp. 2149-2153. 
%-------------------------------------------------------------------------%
%% Code: Pierre-Antoine Thouvenin, [24/06/2016]

[L,N] = size(Y);
R = size(M,2);
K = L - R + 1;
tc = chi2inv(1-pfa,K);

% sigma2 estimation
Yc = bsxfun(@minus,Y,mean(Y,2));
S = Yc*(Yc')/N; % empirical covariance
lambda = svd(S);
sigma2 = mean(lambda(R:end)); % average of the L - R + 1 smallest eigenvalues

Q = bsxfun(@minus,M(:,1:end-1),M(:,end));
E = (eye(L) - Q*((Q'*Q)\(Q')))* bsxfun(@minus,Y,M(:,end));
% E = (eye(L) - Q*inv(Q'*Q)*(Q'))* bsxfun(@minus,Y,M(:,1));
T = sum(E.^2,1)/sigma2;

map = (T > tc); % 1: non-linearity detected

end
