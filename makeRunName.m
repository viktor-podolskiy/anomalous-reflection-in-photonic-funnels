function fileName = makeRunName(study)
    fileName = '';

    if(study.curve == 0)
        fileName = [fileName, 'Linear'];
    elseif(study.curve == 1)
        fileName = [fileName, 'Bezier'];
        fileName = [fileName, '.Br=', num2str(study.bezR)];
        fileName = [fileName, '.Bz=', num2str(study.bezZ)];
        fileName = [fileName, '.Bw=', num2str(study.bezW)];
    end

    fileName = [fileName, '.RBot=', num2str(study.RBot)];
    fileName = [fileName, '.RTop=', num2str(study.RTop)];

    fileName = [fileName, '.bilayers=', num2str(study.numLayers)];

    if(study.funnelType == 2)
        fileName = [fileName, '.Dielectric'];
    else
        if(study.funnelType == 0)
            fileName = [fileName, '.Hyperbolic'];
        elseif(study.funnelType == 1)
            fileName = [fileName, '.Hybrid'];
            fileName = [fileName, '.dielectricRadius=', num2str(study.dielectricRadius)];
        end
        fileName = [fileName, '.Drude'];
        fileName = [fileName, '.lamPl=', num2str(study.lamPlasma)];
        fileName = [fileName, '.metEps0=', num2str(study.metEps0)];
        fileName = [fileName, '.dieEps=', num2str(study.dieEps)];
    end

    if(study.lossRatio ~= 1)
        fileName = [fileName, '.lossRatio=', num2str(study.lossRatio)];
    end

    fileName = [fileName, '.gold=', num2str(study.gold)];
    if(study.gold ~= 0)
        fileName = [fileName, '.goldMinR=', num2str(study.goldMinR)];
    end


    if(length(study.lam0Arr) > 1)
        fileName = [fileName, '.lam0Arr=',num2str(study.lam0Arr(1)),'to',num2str(study.lam0Arr(end)),'by',num2str(study.lam0Arr(2)-study.lam0Arr(1))];
    else
        fileName = [fileName, '.lam0=',num2str(study.lam0Arr)];
    end


end