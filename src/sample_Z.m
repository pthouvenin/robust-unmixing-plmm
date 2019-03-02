function [Z,Lambda,mu] = sample_Z(Y,M,A,dM,Z,Sig,Beta,s,H,W)
%Sample the binary outlier labels z_{n,t} using (56) and (57).
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
% > Sig     noise variances \sig_t^2 [T,1]
% > Beta    granularity parameter in the Potts model [T,1]
% > s       outlier variances \s_t^2 [T,1]
% > H       heigth of the image
% > W       width of the image; (N = W*H)
%
% Outputs:
% < Z       outlier labels [N,T]
% < Lambda  variance for the conditional disctribution of the outliers [T,1]
% < mu      mean of the conditional distribution of the outliers [L,N,T]
%-------------------------------------------------------------------------%
%%

[L,N,T] = size(Y);

% Auxiliary variables (weights)
tmp = Y - mtimesx(bsxfun(@plus,M,dM),A);
Lambda = bsxfun(@times,Sig,s)./bsxfun(@plus,Sig,s); % if NaN values inside Lambda (Inf/Inf) -> increase value of as, bs
mu = bsxfun(@times,tmp,reshape(bsxfun(@rdivide,s,Sig + s),[1,1,T])); % [L,N,T]
Weight = 0.5*bsxfun(@plus,bsxfun(@times,squeeze(sum(tmp.^2,1)),s./((Sig + s).*Sig)),L*(log(Sig./(Sig + s)))); % [N,T]
id = (s == Inf);
Weight(:,id) = -Inf;

%% Checkerboard sampling scheme
I = reshape(1:N,H,W);

% Indices of the "black" elements
temp1 = I(1:2:end,1:2:end);
temp2 = I(2:2:end,2:2:end);
id_b = [temp1(:); temp2(:)];
n_b = numel(id_b); 

% Indices of the "white" elements
temp1 = I(2:2:end,1:2:end);
temp2 = I(1:2:end,2:2:end);
id_w = [temp1(:); temp2(:)];
n_w = numel(id_w);
clear temp1 temp2;

%% Update "black" elements
% Sampling using a checkerboard procedure
Z_im = reshape(Z,[H,W,T]);
id = (Z == 0); % labels in the set {0,1} only -> extension needed for the more general case of K > 2 labels.

% Vertical/horizontal neighbourhood
Dv = (diff(Z_im,1,1) == 0);
Dh = (diff(Z_im,1,2) == 0);
D = reshape(vertcat(zeros(1,W,T),Dv) + vertcat(Dv,zeros(1,W,T)) + horzcat(zeros(H,1,T),Dh) + horzcat(Dh,zeros(H,1,T)),[N,T]); % 4 neighbourhood (Left, Right, Up, Down)

% Potts
Potts_b(:,:,1) = reshape(id(id_b,:).*D(id_b,:),n_b,T);   % labels == 0 in the neighbourhood
Potts_b(:,:,2) = reshape(~id(id_b,:).*D(id_b,:),n_b,T);  % labels == 1 in the neighbourhood

% Exponential term [N,T]
P1_b = 1./(1 + exp(bsxfun(@times,Beta,Potts_b(:,:,1)) - Weight(id_b,:) - bsxfun(@times,Beta,Potts_b(:,:,2)))); % probability that the label is equal to 1
Z(id_b,:) = (rand(n_b,T) < P1_b); % if P0 : Z(id_b,:,:) = (rand(n_b,T) > P0_b)

%% Update "white" elements
% Sampling using a checkerboard procedure
Z_im = reshape(Z,[H,W,T]);
id = (Z == 0); % labels in the set {0,1} only -> extension needed for the more general case of K > 2 labels.

% Vertical/horizontal neighbourhood
Dv = (diff(Z_im,1,1) == 0);
Dh = (diff(Z_im,1,2) == 0);
D = reshape(vertcat(zeros(1,W,T),Dv) + vertcat(Dv,zeros(1,W,T)) + horzcat(zeros(H,1,T),Dh) + horzcat(Dh,zeros(H,1,T)),[N,T]); % 4 neighbourhood (Left, Right, Up, Down)

% Potts
Potts_w(:,:,1) = reshape(id(id_w,:).*D(id_w,:),n_w,T);   % labels == 0 in the neighbourhood
Potts_w(:,:,2) = reshape(~id(id_w,:).*D(id_w,:),n_w,T);  % labels == 1 in the neighbourhood

% Exponential term [N, T]
P1_w = 1./(1 + exp(bsxfun(@times,Beta,Potts_w(:,:,1)) - Weight(id_w,:) - bsxfun(@times,Beta,Potts_w(:,:,2)))); % probability that the label is equal to 1
Z(id_w,:) = (rand(n_w,T) < P1_w);

end
