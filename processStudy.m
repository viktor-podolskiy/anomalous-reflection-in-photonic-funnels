function processStudy(sweepName, deltaZ, options)
    arguments
        sweepName % path to directory of study
        deltaZ double = 50 % distance (nm) above funnel tip to analyze intensity
        options.lam0 double = 0 % wavelength used for distribution plot
        options.funnelHeight = 100 % funnel height used for distribution plot
        options.LineStyle string = "-";
        options.Marker string = "diamond";

    end

    %%
    fileList = load([sweepName, '/fileNameArr.mat']);
    fileNameArr = fileList.fileNameArr;
    mkdir([sweepName,'/','figs']);

    disp(sweepName);

    runs = length(fileNameArr);
    numLayers = zeros(1,runs);
    funnelHtArr = zeros(1,runs);

    displayNames = string.empty;
    data = load(sweepName + "/" + fileNameArr(1));
    lam0Arr = data.lam0Arr;

    subRef = data.subRef;
    geomWd = data.geomWd;
    E0 = sqrt(2*2.99792e8*4*pi*1e-7/(subRef*pi*(geomWd.*1e-6)^2));

    noFunnelI1 = (E0*2*subRef./(1 + subRef)).^2;

    [~, lam0Index] = min(abs(lam0Arr - options.lam0));
    criticalRadiiI1Mat = zeros(runs,length(lam0Arr));
    I1Mat = zeros(runs,length(lam0Arr));

    transmissionMat = zeros(runs,length(lam0Arr));

    dzMeasurement = deltaZ/1000;

    rGrid = data.rGridFields;
    zGrid = data.zGridFields;
    rArr = (rGrid(1,2):(rGrid(1,2)-rGrid(1,1)):rGrid(1,end));
    for runIndex = 1:runs
        disp([num2str(runIndex), ' of ', num2str(runs)]);
        runHue = (runIndex - 1)*(1.0/runs)*0.9;
        runColor = hsv2rgb([runHue,1,0.7]);
        data = load(sweepName + "/" + fileNameArr(runIndex));
        lam0Arr = data.lam0Arr;
        numLayers(runIndex) = data.numLayers;
        funnelHtArr(runIndex) = data.funnelHt;
        transmissionMat(runIndex,:) = data.powTransArr;

        z1 = (data.funnelHt + dzMeasurement).*1e-6 + 0*rArr;
        I = real(cell2mat(data.normEC(lam0Index)));
        I = I.^2;
        I = interp2(rGrid,zGrid,I,rArr, z1);

        rEffI1 = rArr(find(I <= I(1)/exp(2),1));
        if(numel(rEffI1) == 0)
            rEffI1 = Inf;
        end

        criticalRadiiI1Mat(runIndex) = rEffI1;
        I1Mat(runIndex) = I(1);

        figure(1);
        hold on;
        h = semilogy(rArr/(1e-6), I,'LineWidth',2,'Color',runColor,'LineStyle',options.LineStyle);
        h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Layers', repmat({data.numLayers},size(h.XData)) );

        criticalRadiiI1 = zeros(1,length(lam0Arr));
        I1 = zeros(1,length(lam0Arr));
        maxI = zeros(1,length(lam0Arr));
        for lam0ArrIndex=1:length(lam0Arr)
            I = real(cell2mat(data.normEC(lam0ArrIndex)));
            I = I.^2;
            I = interp2(rGrid,zGrid,I,rArr, z1);
            rEffI1 = rArr(find(I <= I(1)/exp(2),1));
            if(numel(rEffI1) == 0)
                rEffI1 = Inf;
            end

            criticalRadiiI1(lam0ArrIndex) = rEffI1;
            criticalRadiiI1Mat(runIndex, lam0ArrIndex) = rEffI1;

            I1(lam0ArrIndex) = I(1);
            maxI(lam0ArrIndex) = max(I);

            I1Mat(runIndex, lam0ArrIndex) = I(1);

        end
        figure(3);
        hold on;
        h = plot(lam0Arr, criticalRadiiI1/(1e-6),'LineWidth',2,'linestyle',options.LineStyle, 'Color',runColor);
        h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Layers', repmat({data.numLayers},size(h.XData)) );

        figure(4);
        hold on;
        h = semilogy(lam0Arr, I1./noFunnelI1,'LineWidth',2,"LineStyle",options.LineStyle, 'Color',runColor);
        h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Layers', repmat({data.numLayers},size(h.XData)) );
        h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Layers', repmat({data.numLayers},size(h.XData)) );


        figure(5);
        hold on
        h = plot(lam0Arr, data.powTransArr./max(data.powTransArr),'LineWidth',2, 'LineStyle',options.LineStyle, 'Color',runColor);
        h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Layers', repmat({data.numLayers},size(h.XData)) );

    end

    figure(1);
    title(['Intensity at ', num2str(lam0Arr(lam0Index)),' µm']);
    xlabel('r (µm)');
    set(gca, 'FontSize', 18);

    figure(3);
    title 'Effective Focal Radius';
    xlabel('λ (µm)');
    ylabel('Effective Radius (µm)');
    ylim([0,2]);
    xlim([lam0Arr(1),lam0Arr(end)]);
    xticks((lam0Arr(1):2:lam0Arr(end)))
    grid on
    set(gca, 'FontSize', 18);
    exportgraphics(gcf, [sweepName, '/figs/EffectiveFocalRadiusByHeight', '.png'])

    figure(4);
    title 'Enhancement';
    xlabel('λ (µm)');
    ylabel('$I/I_0$','Interpreter','latex');
    ylim('padded');
    xlim([lam0Arr(1),lam0Arr(end)]);
    xticks((lam0Arr(1):2:lam0Arr(end)))
    grid on
    set(gca, 'FontSize', 18);
    exportgraphics(gcf, [sweepName, '/figs/MaxIntensityByHeight', '.png'])


    figure(5);
    xlabel('λ (µm)');
    ylabel('Transmission (AU)');
    ylim('padded');
    ylims = ylim();
    ylim([0,ylims(2)]);
    xlim([lam0Arr(1),lam0Arr(end)]);
    xticks((lam0Arr(1):2:lam0Arr(end)))
    set(gca, 'FontSize', 16);
    exportgraphics(gcf, [sweepName, '/figs/TransmissionByHeight', '.png'])

    figure(2);
    hold on;
    displayNames = string.empty;
    for runIndex = 1:runs
        runHue = (runIndex - 1)*(1.0/runs)*0.9;
        runColor = hsv2rgb([runHue,1,0.7]);
        h = semilogy(criticalRadiiI1Mat(runIndex,lam0Index)/(1e-6), I1Mat(runIndex,lam0Index)./noFunnelI1, 'Marker',options.Marker,'Color',runColor, 'MarkerFaceColor',runColor);
        h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Layers', numLayers(runIndex));
        if(mod(runIndex,4) == 1)
            displayNames(end+1) = num2str(numLayers(runIndex));
        else
            displayNames(end+1) = '';
        end
        hold on;

    end

    figure(2);
    title(['At ', num2str(lam0Arr(lam0Index)), ' µm']);
    xlabel('Effective Radius (µm)');
    ylabel('I_{max}/I_0');
    xlim([0,1.1*max(xlim)]);
    ylim('padded');
    lgd = legend(displayNames);
    lgd.Title.String = 'numLayers';
    lgd.Location = 'northwest';
    set(gca, 'FontSize', 18);
    exportgraphics(gcf, [sweepName, '/figs/maxIntensityVFocalRadius', num2str(lam0Arr(lam0Index)), '.png']);

    figure(6);
    hold off;
    surf(lam0Arr,funnelHtArr,I1Mat./noFunnelI1,'EdgeColor','none');
    xlim('tight');
    ylim('tight');
    X = lam0Arr;
    Y = criticalAngle(data.study, X);
    Y = data.funnelRBot./tan(Y)./(1-data.funnelRTop/data.funnelRBot);
    Z = 0.*X + 1e27;
    xlimits = xlim;
    ylimits = ylim;
    hold on
    plot3(X,Y,Z,'Color',[0.1,1,0.2],'LineWidth',3,'LineStyle','--'); % Critical angle overlay
    view(2);
    xlim(xlimits)
    ylim(ylimits)
    box on;
    set(gca,'Layer','top')
    grid off
    set(gca, 'ColorScale', 'log');
    title 'Enhancement';
    xlabel('λ_0 (µm)');
    ylabel('h (µm)');
    colorbar;
    set(gca, 'FontSize', 18);

    figure(7);
    hold off;
    surf(lam0Arr,funnelHtArr,criticalRadiiI1Mat*1e6,'EdgeColor','none');
    hold on
    plot3(X,Y,Z,'Color',[0.1,1,0.2],'LineWidth',3,'LineStyle','--'); % Critical angle overlay
    view(2);
    xlim(xlimits)
    ylim(ylimits)
    box on;
    set(gca,'Layer','top')
    grid off
    % set(gca, 'ColorScale', 'log');
    title 'Confinement';
    xlabel('λ_0 (µm)');
    ylabel('h (µm)');
    colorbar;
    set(gca, 'FontSize', 18);

    figure(8);
    
    [~, fileIndex] = min(abs(funnelHtArr - options.funnelHeight));
    data = load(sweepName + "/" + fileNameArr(fileIndex));
    um = 1e-6;
    raFields = 0*um; % start of r
    rbFields = 4*um; % end of r
    drFields = 0.01*um; % increment of r

    zaFields = 0*um; % start of z
    zbFields = 5*um; % end of z
    dzFields = 0.01*um; % increment of z

    rArrFields =(raFields:drFields:rbFields);
    zArrFields = (zaFields:dzFields:zbFields);

    [rGridFields,zGridFields] = meshgrid(rArrFields,zArrFields);
    normE2Plot = real(interp2(rGrid,zGrid,cell2mat(data.normEC(lam0Index)),rGridFields, zGridFields,'spline')).^2;
    hold off
    surf(rGridFields/(1e-6), zGridFields/(1e-6), normE2Plot./noFunnelI1,'EdgeColor','none');
    
    hold on
    
    view(2);
    colorbar;
    
    ylim('tight');
    xlabel('r (µm)');
    ylabel('z (µm)');
    set(gca,'FontSize', 18);
    set(gca,'Layer', 'top');
    box on;
    grid off;
    daspect([1,1,1]);
    set(gca, 'FontSize', 18);



    save([sweepName,'/processedData.mat'],'I1Mat','criticalRadiiI1Mat', 'transmissionMat', 'runs','numLayers','funnelHtArr','lam0Arr','rGrid','zGrid', 'subRef', 'geomWd');

end