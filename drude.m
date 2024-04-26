function [metEps] = drude(gamma, metEps0, lamP, lam0Arr) 
%loads permittivity data for dielectric and anisotropic metal layers
c = 299792458;

wp = 2*pi*c/lamP/1e-6;

wArr = 2*pi*c./lam0Arr/1e-6;
metEps = conj(metEps0*(1-wp^2./(wArr.^2+1i*wArr.*gamma)));


end