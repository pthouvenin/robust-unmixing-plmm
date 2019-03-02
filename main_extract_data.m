%% Create .mat files (matrix cube) for the image time series
% Extract the hyperspectral data cube from the raw data files contained in
% data/raw_data.
%-------------------------------------------------------------------------%
%%
% Author : P.A. Thouvenin [18/02/2016]
%-------------------------------------------------------------------------%
clc; clear all; close all
format compact;
addpath utils
addpath data/raw_data
%-------------------------------------------------------------------------%
%%

% Extract time series parameters (from the first image)
FileName = strcat('Im1','.hdr');
PathName = strcat(cd,'/data/raw_data/');
[W, H, L_true, interleave, offset, byte_order, data_type, wavelength, wavelength_unit] = extract_parameters(PathName,FileName); 
data_type = strcat(data_type,'=>','double');
N = H*W;

% Set extra parameters for data extraction
T = 6; % number of datasets (time instants)
mask = [3:102,115:150,185:213,218:221]; % spectral mask
L = numel(mask);

Y = zeros(L,N,T);

for t = 1:T
	FileName = strcat('Im',num2str(t));
	data = multibandread(fullfile(PathName,FileName),[H,W,L_true],data_type,offset,interleave,byte_order,{'Column',1:W},{'Row',1:H},{'Band',1:L_true});
	Y(:,:,t) = reshape(data(:,:,mask),N,L)'/10000; % column-major ordering
%     figure;
%     plot(Y(:,:,t))
	
	%--
    % Debug tools to determine apporpriate number of endmembers and see if
	% poor bands with poor SNR are properly removed: uncomment if needed
	%--
	% R = 3; % endmember numbers
	% [lambda] = scree_plot(Y);
	% 
	% [w,Rw] = estNoise(Y,'additive','off');
	% [K,Es] = hysime(Y,w,Rw,'off');
	% [M, V, U, Y_bar, endm_proj, Y_proj] = find_endm(Y,R,'vca');
	% A = sunsal(M,Y,'POSITIVITY','yes','ADDONE','yes'); % fully constrained least squares by ADMM
	% A = max(bsxfun(@minus,A,max(bsxfun(@rdivide,cumsum(sort(A,1,'descend'),1)-1,(1:R)'),[],1)),0);
	%--
end

save('data/Series_mcmc_150','Y','H','W','L','L_true','wavelength','wavelength_unit','mask','-v7.3');
disp('... DONE');
disp('---------------------------------------------------------------------------');
