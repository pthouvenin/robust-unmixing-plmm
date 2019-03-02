function z = sample_labels(w)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% Debug: 18/02/16

% w : [N,T,K] probability of belonging to a given class k, K number of
% classes.

% Draw zn in {1,...,K} with probability {w1,...,wK}
u = unifrnd(0,1,[size(w,1),size(w,2)]); % [ok]
W(:,:,1) = 1./(1 + exp(w(:,:,2) - w(:,:,1)));
W(:,:,2) = 1./(1 + exp(w(:,:,1) - w(:,:,2)));
t = cumsum(W,3);
z = sum(bsxfun(@gt,u,t(:,:,end-1)),3); % (u > t) [label values from 0 to K-1]

end
