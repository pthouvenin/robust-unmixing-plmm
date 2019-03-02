function [M_mse,A_mse,dM_mse,X_mse,Z_map,Sig_mse,s_mse,Psi_mse,beta_mse,accept_rate,beta_c,Sigc,s_beta_c] = rPLMM_real_data(Y,M,A,dM,X,Z,Sig,s,P,U,y_bar,H,W,Parameter,Hyperparameter)
%Robust unmixing based on the Perturbed Linear Mixing Model (PLMM).
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
% > X       outlier terms [L,N,T] 
% > Z       current labels [N,T]
% > Sig     noise variances \sig_t^2 [1,T]
% > s       outlier variances [1,T]
% > P       projector onto the PCA-subspace [K,L]
% > U       back projector onto the orignal space [L,K]
% > y_bar   empirical average of the data over space and time [L,1]
% > H       heigth of the image
% > W       width of the image; (N = W*H)
% > Parameter  structure containing the following algorithmic parameters
% >> nBi          number of brun-in iterations
% >> nMCMC        total number of MCMC iterations
% >> update_freq  frequency to update the variance of the proposals 
%                 in the \beta update step
% > Hyperparameter  structure containing the following hyperparameters
% >> nu_dM   variability prior variance for dM(:,:,1) [1]
% >> Psi     variability prior variances (temporal smoothness) [L,R]
% >> Eps     abundance prior variance (see (8)) (temporal smoothness) [1,N]
% >> xi      hyperparameter: prior variance of the endmember matrix [1]
% >> as,bs   hyperparameters for the IG prior of the outlier variances [1]
% >> asig,bsig  hyperparameters for the IG prior of the noise variances [1]
% >> beta    initial value granularity parameter in the Potts model [1]
% >> s_beta  variances of the Gaussian proposals in the update of Beta [1]
%
% Outputs:
% > M_mse        MMSE estimate of the endmembers [L,R]
% > A_mse        MMSE estimate of the abundances [R,N,T]
% > dM_mse       MMSE estimate of the variability [L,R,T]
% > X_mse        MMSE estimate of the outlier terms [L,N,T]
% > Z_map        mMAP estimate of the outlier labels [N,T]
% > Sig_mse      MMSE estimate of the noise variances [1,T]
% > s_mse        MMSE estimate of the outlier variances [1]
% > Psi_mse      MMSE estimate of the variability variances [L,R]
% > beta_mse     MMSE estimate of the granularity parameters \beta_t [1,T]
% > accept_rate  final acceptance rate for the samlped vlaues of beta
% > beta_c       full chain of the sampled beta values
% > Sigc         full chain of the sampled Sig parameter
% > s_beta_c     full chain of the s_beta parameter
%
%-------------------------------------------------------------------------%
%%
% Reference:
% P.-A. Thouvenin, N. Dobigeon, J.-Y. Tournret, "A hierarchical Bayesian 
% model accounting for endmember variability and abrupt spectral changes 
% to unmix multitemporal hyperspectral images," IEEE Trans. Comput. Imag.,
% vol. 4, no. 1, pp. 32-45, Mar. 2018.
%-------------------------------------------------------------------------%
%% Hybrid Gibbs sampler

[L,N,T] = size(Y);
R = size(M,2);

% Initialization / memory allocation
M_mse = zeros(L,R);
A_mse = zeros(R,N,T);
dM_mse = zeros(L,R,T);
X_mse = zeros(L,N,T);
Z_mse = zeros(N,T);
Sig_mse = zeros(1,T);
beta_mse = zeros(1,T);
Psi_mse = zeros(L,R);
s_mse = 0;
Sigc = zeros([Parameter.nMCMC,T]);
beta_c = zeros([Parameter.nMCMC,T]);

% Parameters to update the granularity parameter
lgth = floor(Parameter.nBi/Parameter.update_freq);
k = 1;
counter = 0;
s_beta = Hyperparameter.s_beta*ones(1,T);
s_beta_c = zeros(lgth,T);
accept_rate = zeros(lgth,T);

% Parameters to sample the granularity parameter beta
Cn = cos((2*(1:N)' - 1)*pi/(2*N)); % [N,1]
Beta = Hyperparameter.beta*ones(1,T);

% Hyperparameters to sample the variability dM
Psi = Hyperparameter.Psi;

% For projection onto the subspace (after PCA)
Up = (U > 0);
Um = (U < 0);
idp = any(Up,1);
idm = any(Um,1);

for q = 1:Parameter.nMCMC
    
    % Sample M (endmembers)
    % M = Gibbs_MdMp(Y,M,A,dM,X,Sig,Hyperparameter.xi);
    M = sample_M_subspace(Y,M,A,dM,X,Sig,Hyperparameter.xi,U,P,Up,Um,idp,idm,y_bar);

    % Sample dM (variability)
    % dM = sample_dM_01(Y,M,A,dM,X,Sig,Psi,nu_dM); % Hyperparameter.nu_dM % constraint [0, 1]
    dM = sample_dM(Y,M,A,dM,X,Sig,Psi,Hyperparameter.nu_dM); % non-negativity constraint only, to be used when the projection onto subspace

    % Sample A (abundance coefficients)
    A = sample_A(Y,M,A,dM,Z,X,Hyperparameter.Eps,Sig); % [ok: 30/03/16]

    % Sample Z (label map)
    [Z,Lambda,mu] = sample_Z(Y,M,A,dM,Z,Sig,Beta,s,H,W); % [ok: 25/04/16]

    % Sample X (outlier contribution)
    X = sample_X(mu,Lambda,Z);

    % Sample s_t^2 (outlier variance)
    s = sample_s(Hyperparameter.as,Hyperparameter.bs,X,Z);

    % Sample \sigma_t^2 (noise variance)
    Sig = sample_Sig(Y,M,A,dM,X,Hyperparameter.asig,Hyperparameter.bsig);

    % Sample \psi_{\ell,r}^2 (variability variance)
    Psi = sample_Psi(dM,Hyperparameter.as,Hyperparameter.bs,T); % a_psi = b_psi = 1e-3 = as = bs

    % Sample nu_dM (only when the constraint M + dM >= 0 is not active, difficult otherwise)
    % nu_dM = sample_nu(dM(:,:,1),Hyperparameter.as,Hyperparameter.bs,L,R);

    % Sample beta
    [Beta,id_accept] = sample_beta(Beta,Cn,Z,H,W,s_beta); % ok [22/06/16]

    % Update variance for the Beta proposal
    if (q <= Parameter.nBi)
        counter = counter + 1;
        if k <= lgth
            accept_rate(k,:) = accept_rate(k,:) + id_accept;
        end
        if (counter == Parameter.update_freq)
            % Acceptance ratio
            accept_rate(k,:) = accept_rate(k,:)/Parameter.update_freq;
            % Parameter adjustment
            id = (accept_rate(k,:) < 0.4);
            s_beta(id) = 0.8*s_beta(id); 
            id = (accept_rate(k,:) > 0.6);
            s_beta(id) = 1.25*s_beta(id);
            s_beta_c(k,:) = s_beta;
            counter = 0;
            k = k+1;
        end
    else % compute MMSE estimates (for q > nBi, aggregate information)
        M_mse = M_mse + M;
        A_mse = A_mse + A;
        dM_mse = dM_mse + dM;
        X_mse = X_mse + X;
        Z_mse = Z_mse + Z;
        Sig_mse = Sig_mse + Sig;
        s_mse = s_mse + s;
        beta_mse = beta_mse + Beta;
        Psi_mse = Psi_mse + Psi;
    end

    % Convergence monitoring
    Sigc(q,:) = Sig;
    beta_c(q,:) = Beta;

end

% Normalize MMSE estimates
M_mse = M_mse/(Parameter.nMCMC - Parameter.nBi);
A_mse = A_mse/(Parameter.nMCMC - Parameter.nBi);
dM_mse = dM_mse/(Parameter.nMCMC - Parameter.nBi);
X_mse = X_mse/(Parameter.nMCMC - Parameter.nBi);
Z_map = (Z_mse > (Parameter.nMCMC - Parameter.nBi)/2);
Sig_mse = Sig_mse/(Parameter.nMCMC - Parameter.nBi);
s_mse = s_mse/(Parameter.nMCMC - Parameter.nBi);
beta_mse = beta_mse/(Parameter.nMCMC - Parameter.nBi);
Psi_mse = Psi_mse/(Parameter.nMCMC - Parameter.nBi);

end
