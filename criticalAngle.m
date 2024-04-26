function [cAngle] = criticalAngle(study, lam0Arr)


    % [dieEps, metEpsRPhi, metEpsZZ] = loadEps(study.dataName, lam0Arr, study.lamPlasma);
metEps = drude(study.gamma, study.metEps0, study.lamPlasma, lam0Arr);
dieEps = study.dieEps;

unitHt = study.dieHt + study.metHt;
epsRPhi = (dieEps*study.dieHt+metEps*study.metHt)/(unitHt);
epsZZ = unitHt./(study.dieHt./dieEps+study.metHt./metEps);

cAngle = real(atan(sqrt(-(real(epsRPhi)./real(epsZZ)))));

end