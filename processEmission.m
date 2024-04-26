function processEmission(sweepName, options)

    arguments
        sweepName % path to directory of study
        options.lam0 double = 0 % wavelength used for distribution plot
        options.funnelHeight = 100 % funnel height used for distribution plot
    end

    fileNameArr = load([sweepName,'/','fileNameArr.mat'],'fileNameArr').fileNameArr;
    ref = load([sweepName,'/','purcelReference.mat'],'lam0Arr','emittedArr','emittedBottomArr');
    study = load([sweepName,'/','purcelReference.mat'],'study').study;
    mkdir([sweepName,'/','figs']);

    lam0Arr = ref.lam0Arr;

    critAngleX = lam0Arr;
    critAngleY = criticalAngle(study, critAngleX);
    critAngleY = study.RBot./tan(critAngleY);
    critAngleZ = 0.*critAngleX + 1e27;

    funnelHtArr = zeros(1,length(fileNameArr));
    purcMat = zeros(length(fileNameArr),length(lam0Arr));
    transMat = zeros(length(fileNameArr),length(lam0Arr));
    for fileNameArrIndex = 1:1:length(fileNameArr)
        fileName = fileNameArr(fileNameArrIndex);
        data = load(sweepName + "/" + fileName,'lam0Arr','emittedArr','emittedBottomArr','funnelHt');
        funnelHtArr(fileNameArrIndex) = data.funnelHt;
        purc = data.emittedArr./ref.emittedArr;
        trans = -data.emittedBottomArr./data.emittedArr;
        purcMat(fileNameArrIndex,:) = purc;
        transMat(fileNameArrIndex,:) = trans;

    end

    [~, lam0Index] = min(abs(lam0Arr - options.lam0));
    [~, fileIndex] = min(abs(funnelHtArr - options.funnelHeight));
    %%
    htArr = study.numLayersArr*(study.metHt + study.dieHt)./1000;
    [lamGrid, htGrid] = meshgrid(lam0Arr, htArr);

    figure(1)
    clf
    surf(lamGrid, htGrid, purcMat,'EdgeColor','none')
    view(2)
    set(gca,'ColorScale','log')
    colorbar
    hold on
    plot3(critAngleX,critAngleY,critAngleZ,'Color',[0.1,1,0.2],'LineWidth',3,'LineStyle','--');
    ylim([min(htArr),max(htArr)])
    grid off;
    box on;
    set(gca,'Layer','top')
    xlabel('λ_0 (µm)');
    ylabel('h (µm)');
    set(gca, 'FontSize', 18);
    xlim(lam0Arr([1,end]))
    ylimits = ylim;
    clim([1e2,1e3])
    xticks(8:2:20)
    exportgraphics(gcf, [sweepName, '/figs/purcell_test.png'],'Resolution',600);

    figure(2)
    clf
    surf(lamGrid, htGrid, transMat,'EdgeColor','none')
    view(2)
    set(gca,'ColorScale','log')
    colorbar
    hold on
    plot3(critAngleX,critAngleY,critAngleZ,'Color',[0.1,1,0.2],'LineWidth',3,'LineStyle','--');
    ylim([min(htArr),max(htArr)])
    grid off;
    box on;
    set(gca,'Layer','top')
    xlabel('λ_0 (µm)');
    ylabel('h (µm)');
    set(gca, 'FontSize', 18);
    xlim(lam0Arr([1,end]))
    xticks(8:2:20)
    exportgraphics(gcf, [sweepName, '/figs/extraction_red_test.png'],'Resolution',600);

    bestData = load(sweepName + "/" + fileNameArr(fileIndex),'normEC','ErC', 'EphiC', 'EzC', 'BrC', 'BphiC', 'BzC', 'rGridFields','zGridFields');
    mu0 = 4*pi*1e-7;
    
    normE = bestData.normEC{lam0Index};
    %%
    um = 1e-6;
    raFields = 0*um; % start of r
    rbFields = max(max(bestData.rGridFields)); % end of r
    nRFields = 100;
    drFields = (rbFields-raFields)/(nRFields-1); % increment of r

    zaFields = 0*um; % start of z
    zbFields = 6*um; % end of z
    nZFields = 200;
    dzFields = (zbFields - zaFields)/(nZFields-1); % increment of z

    rArrFields =(raFields:drFields:rbFields);
    zArrFields = (zaFields:dzFields:zbFields);

    [rGridFields,zGridFields] = meshgrid(rArrFields,zArrFields);
    normE2Plot = abs(interp2(bestData.rGridFields,bestData.zGridFields,normE,rGridFields, zGridFields,'spline')).^2;

    figure(3)
    surf(rGridFields/(1e-6), zGridFields/(1e-6), normE2Plot./max(max(normE2Plot)),'EdgeColor','none');
    view(2)
    set(gca,'ColorScale','log')
    colorbar
    xlim('tight')
    daspect([1,1,1])
    box on
    grid off
    set(gca,'Layer','top')

    xlabel('r (µm)');
    ylabel('z (µm)');
    set(gca,'FontSize', 18);
    exportgraphics(gcf, [sweepName, '/figs/Distribution', num2str(lam0Arr(lam0Index)) ,'um.', 'layers.png'], "Resolution",600);


end