function nu_dM = sample_nu(dM1,a_nu,b_nu,L,R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nu_dM = 1./gamrnd(a_nu + L*R/2, 1./bsxfun(@plus,b_nu,0.5*sum(dM1(:).^2))); % [L|R]

end

