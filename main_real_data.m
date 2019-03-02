%%
%=========================================================================%
%       TEMPORAL VARIABILITY IN HYPERSPECTRAL UNMIXING  (rPLMM-MCMC)      %
%=========================================================================%
%% FILE DESCRIPTION
% File : main_real_data.m
% Author : P.A. Thouvenin [18/02/2016]
% Last modified : [02/03/2019]
clc, clear all, close all, format compact;
addpath data;
addpath utils; addpath src; addpath src_real_data;
addpath lib
addpath mtimesx/src;
%-------------------------------------------------------------------------%
%%
% Reference:
% P.-A. Thouvenin, N. Dobigeon, J.-Y. Tournret, "A hierarchical Bayesian 
% model accounting for endmember variability and abrupt spectral changes 
% to unmix multitemporal hyperspectral images," IEEE Trans. Comput. Imag.,
% vol. 4, no. 1, pp. 32-45, Mar. 2018.
%-------------------------------------------------------------------------%
%% REMARK:
% Reshape hyperspectral data cube "data" into a matrix Y:
% - MATLAB ordering (column-wise): [H,W,L] -> [L,H*W] Y = reshape(data,H*W,L)';
%                                  [L,H*W] -> [H,W,L] data = reshape(Y',H,W,L);
%=========================================================================%
%%
% dbstop in file if expression (naninf)
% dbclear all
load('Series_mcmc_150','Y','H','W');
proj_l1ball_array = @(y) max(bsxfun(@minus,abs(y),max(max(bsxfun(@rdivide,cumsum(sort(abs(y),1,'descend'),1)-(1-1e-5),(1:size(y,1))'),[],1),0)),0).*sign(y);
%--------------------------------------------------------------
% General parameters
%--------------------------------------------------------------
T = 6;          % number of time instants
R = 3;          % number of endmembers
L = size(Y,1);  % number of spectral bands
N = H*W;        % total number of pixels
Parameter = struct('nMCMC',500,'nBi',450,'update_freq',20);
Hyperparameter = struct('beta',1.7,'xi',1,'as',1e-3,'bs',1e-3,'Psi',1e-2*ones(L,R),'Eps',1e-2*ones(1,N),'nu_dM',1e-5,'s_beta',1e-4);

%--------------------------------------------------------------
% Initialization (uncomment to initialize the method, comment once initialized)
%--------------------------------------------------------------
% M = vca(reshape(Y,[L,T*N]),'Endmembers',R,'verbose','off');
% M = sisal(reshape(Y,[L,T*N]),R,'spherize','yes','VERBOSE',0);%,'TOLF',1e-5,'MM_ITERS',500);
[U,Ut,P,Pt,y_bar,Y_proj,d] = pca(reshape(Y,[L,T*N]),R-1);
M = vca(Y(:,:,3),'Endmembers',R,'verbose','off'); % initialize M from the third image of the sequence (endmembers relatively well represented)
M = abs(M);

A = zeros(R,N,T);
for t = 1:T
    % FCLS by ADMM
    A(:,:,t) = sunsal(M,Y(:,:,t),'POSITIVITY','yes'); % 'ADDONE','yes'
    %A(:,:,t) = max(bsxfun(@minus,A(:,:,t),max(bsxfun(@rdivide,cumsum(sort(A(:,:,t),1,'descend'),1)-1,(1:R)'),[],1)),0);
    %A(:,:,t) = proj_l1ball_array(A(:,:,t));
end
% 
A = bsxfun(@rdivide,abs(A),sum(A,1) + 1e-3);
dM = zeros(L,R,T);
X = zeros(L,N,T);
Z = false(N,T);
Sig = 1e-4*ones(1,T);
s = 5e-3*ones(1,T);

save('Init_rd_v150','M','A','dM','X','Z','Sig','s','U','P','y_bar')
return

%--------------------------------------------------------------
% Unmixing
%--------------------------------------------------------------
load('Init_rd_v150.mat');
Hyperparameter.asig = 1e-3;
Hyperparameter.bsig = 1e-3;

tic
[M_mse,A_mse,dM_mse,X_mse,Z_map,Sig_mse,s_mse,Psi_mse,beta_mse,accept_rate,beta_c,Sigc,s_beta_c] = rPLMM_real_data(Y,M,A,dM,X,Z,Sig,s,P,U,y_bar,H,W,Parameter,Hyperparameter);
run_time = toc;

%--------------------------------------------------------------
% Reconstruction errors
%--------------------------------------------------------------
RE = 0.5*sum(sum(sum((Y - mtimesx(bsxfun(@plus,M_mse,dM_mse),A_mse) - X_mse).^2,1),2),3)/numel(Y);

% Save estimates
save('rd_v150','Hyperparameter','A_mse','M_mse','dM_mse','X_mse','Z_map','Sig_mse','Sigc','s_mse','accept_rate','Parameter','run_time','beta_mse','beta_c','s_beta_c','Psi_mse','RE', '-v7.3');
 
figure;
for t = 1:T
    hold on;
    plot(M_mse);
    hold on; plot(M_mse + dM_mse(:,:,t),'-.')
end

for r = 1:R
    figure;
    for t = 1:T
        hold on;
        plot(M_mse(:,r));
        hold on; plot(M_mse(:,r) + dM_mse(:,r,t),'-.')
    end
end
