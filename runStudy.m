function runStudy(study, model)
    %% setup sweep
    % input arrays
    %lam0Arr = (4:0.5:10);
    lam0Arr = study.lam0Arr;
    % output arrays
    normEC = cell(1,length(lam0Arr));
    ErC = normEC;
    EphiC = normEC;
    EzC = normEC;
    BrC = normEC;
    BphiC = normEC;
    BzC = normEC;
    epsRRC = normEC;
    epsPhiPhiC = normEC;
    epsZZC = normEC;
    powTransArr = 0*lam0Arr;
    powReflArr = 0*lam0Arr;
    powTransFarArr = 0*lam0Arr;
    powReflFarArr = 0*lam0Arr;

    emittedArr = 0*lam0Arr;
    emittedRightArr = 0*lam0Arr;
    emittedTopArr = 0*lam0Arr;
    emittedBottomArr = 0*lam0Arr;

    %% sweep
    funnelCurve = study.curve;   % controls funnel curve, 0 is straight line, 1 is bezier curve
    
    lossRatio = study.lossRatio; % fraction of loss used

    funnelType = study.funnelType;%1; % 0 for hyperbolic, 1 for hybrid, 2 for dielectric
    % dielectricLam0 = study.dielectricLam0;%0; % wavelength for which dielectric fills up to critical radius

    gold = study.gold;
    goldMinR = study.goldMinR;

    funnelRBot = study.RBot; % [um]
    funnelRTop = study.RTop; % [um]

    numLayersArr = study.numLayersArr;
    numRuns = length(numLayersArr);
    fileNameArr = strings(length(numLayersArr),1);
    sweepName = study.sweepName;

    mkdir(study.sweepName);
    
    if(study.sim == 0)
        % port input
        model.physics('ewfd').feature('port1').active(1);
        model.physics('ewfd').feature('ecd1').active(0);
        Ex = num2str(study.port.Ex);
        Ey = num2str(study.port.Ey);
        Ez = num2str(study.port.Ez);
        model.physics('ewfd').feature('port1').set('E0', {Ex, Ey, Ez});
        model.param().set('emitterEps', '1');
    elseif(study.sim == 1)
        % emitter input
        if(study.emitter.polarization == 'X')
            study.m = -1;
            study.emitter.Jx = 1;
            study.emitter.Jy = 1i;
            study.emitter.Jz = 0;
        elseif(study.emitter.polarization == 'Z')
            study.m = 0;
            study.emitter.Jx = 0;
            study.emitter.Jy = 0;
            study.emitter.Jz = 1;
        else
            error("Emitter polarization must be X or Z: study.emitter.polarization");
        end
        model.physics('ewfd').feature('port1').active(0);
        model.physics('ewfd').feature('ecd1').active(1);
        Jx = num2str(study.emitter.Jx);
        Jy = num2str(study.emitter.Jy);
        Jz = num2str(study.emitter.Jz);
        model.physics('ewfd').feature('ecd1').set('Je', {Jx, Jy, Jz});
        model.param().set('emitterEps', num2str(study.dieEps));
        numRuns = numRuns+1; % Do an extra run with all air except for emitter
    end
    model.physics('ewfd').prop('outofplanewavenumber').set('mFloquet', study.m);

    for numLayersArrIndex = 1:numRuns
        disp(['Run: ', num2str(numLayersArrIndex), ' of ', num2str(numRuns)]);
        tStart = tic;
        for lam0Index = 1:length(lam0Arr)
            disp(['  lam0Index: ', num2str(lam0Index), ' of ', num2str(length(lam0Arr))]);
            %% model setup
            % study properties
            lam0 = lam0Arr(lam0Index); % vacuum wavelength [um]
            subHt = 10; % substrate thickness [um]
            supHt = 10; %45; % superstrate thickness [um]
            geomWd = 15; %20; % width of simulation space [um]
            pmlSz = 10; % thickness of PMLs [um]

            % material properties
            % load permitivity data from file
            
            metEps = drude(study.gamma, study.metEps0, study.lamPlasma, lam0);
            dieEps = study.dieEps;
            metEps = real(metEps) + lossRatio*1i*imag(metEps); % just for low loss testing

            subRef = study.subRef;%sqrt(epsSubstrate(lam0)); % refractive index of substrate
            supRef = 1; %refractive index of superstrate

            % mesh properties
            %meshMaxMinRatio = 5; % maximum grid element size / minimum grid element size

            meshAirRatio = 10; % vacuum wavelength / maximum grid element size of Air
            meshAirMax = lam0 / meshAirRatio; % maximum grid element size of Air
            %meshAirMin = meshAirMax / meshMaxMinRatio; % minimum grid element size of Air

            meshSubstrateRatio = 10*subRef; %45; %20; % vacuum wavelength / maximum grid element size of Size
            meshSubstrateMax = lam0 / meshSubstrateRatio; % maximum grid element size of Size
            %meshSubstrateMin = meshSizeMax / meshMaxMinRatio; % minimum grid element size of Size

            meshFunnelRatio = 50; % vacuum wavelength / maximum grid element size of Funnel
            meshFunnelMax = lam0 / meshFunnelRatio; % maximum grid element size of Funnel
            %meshSmallMin = meshSmallMax / meshMaxMinRatio; % minimum grid element size of Funnel

            meshSmallRatio = 1600; % vacuum wavelength / maximum grid element size of Small
            meshSmallMax = lam0 / meshSmallRatio; % maximum grid element size of Small
            %meshSmallMin = meshSmallMax / meshMaxMinRatio; % minimum grid element size of Small

            % funnel geometry
            %funnelRBot = 3; % radius at bottom of funnel [um]
            %funnelRTop = 0.125; % radius at top of funnel [um]

            % smallMeshWd = 2*funnelRTop; % width of small mesh region % NOW SET WITH COMSOL
            % smallMeshHt = smallMeshWd; % height of small mesh region % NOW SET WITH COMSOL

            funLayers = numLayersArr(min(numLayersArrIndex,length(numLayersArr)));
            numLayers = numLayersArr(min(numLayersArrIndex,length(numLayersArr))); % number of bilayers in funnel
            study.numLayers = numLayers;

            dieHt = study.dieHt; % thickness of dielectric portion of bilayer [nm]
            metHt = study.metHt; % thickness of metal portion of bilayer [nm]
            funnelHt = (dieHt + metHt)*(numLayers)/1000; % total height of funnel [um]

            goldOffset = 1*dieHt/1000; % distance from top of gold to top of funnel [um]
            goldHt = funnelHt - goldOffset; % height of gold cladding [um]
            %gold = true; % whether or not gold cladding is present

            bezierR = study.bezR; % radius of bezier control point [um]
            % study.bezZ = study.bezZratio*funnelHt;
            bezierZ = study.bezZ; % height of bezier control point [um]
            bezierW = study.bezW; % weight of bezier control point


            if(funnelType == 0)
                dielectricRadius = 3*funnelRBot;
            elseif(funnelType == 1)
                % dielectricRadius = alpha * dielectricLam0 / (2 * pi * study.subRef);
                dielectricRadius = study.dielectricRadius;
            elseif(funnelType == 2)
                dielectricRadius = 0;
            end
            funnelDielectricHeight = supHt;

            if(numLayersArrIndex == length(numLayersArr) + 1)
                dieEps = 1;
                metEps = 1;
                subRef = 1;
                supRef = 1;
                gold = 0;
            end
            % set parameters and build model
            model.param().set('lam0', [num2str(lam0), '[um]']);
            model.param().set('supHt', [num2str(supHt), '[um]']);
            model.param().set('subHt', [num2str(subHt), '[um]']);
            model.param().set('geomWd', [num2str(geomWd), '[um]']);
            model.param().set('pmlSz', [num2str(pmlSz), '[um]']);

            %
            model.param().set('dieEps', num2str(dieEps));
            model.param().set('metEps', num2str(metEps));

            model.param().set('subRef', num2str(subRef));
            model.param().set('supRef', num2str(supRef));


            %
            %model.param().set('meshAirMin', [num2str(meshAirMin), '[um]']);
            model.param().set('meshAirMax', [num2str(meshAirMax), '[um]']);

            %model.param().set('meshSubstrateMin', [num2str(meshSubstrateMin), '[um]']);
            model.param().set('meshSubstrateMax', [num2str(meshSubstrateMax), '[um]']);

            %model.param().set('meshFunnelMin', [num2str(meshFunnelMin), '[um]']);
            model.param().set('meshFunnelMax', [num2str(meshFunnelMax), '[um]']);

            %model.param().set('meshSmallMin', [num2str(meshSmallMin), '[um]']);
            model.param().set('meshSmallMax', [num2str(meshSmallMax), '[um]']);

            %model.mesh('mesh1').feature('size').set('hgrad', '1.1');

            %
            model.param().set('funnelRBot', [num2str(funnelRBot),'[um]']);
            model.param().set('funnelRTop', [num2str(funnelRTop),'[um]']);

            model.param().set('emitterActiveHt', [num2str(study.emitter.activeHt), '[nm]']);
            model.param().set('emitterPassiveHt', [num2str(study.emitter.passiveHt), '[nm]']);

            model.param().set('numLayers', num2str(numLayers));

            model.param().set('dieHt', [num2str(dieHt), '[nm]']);
            model.param().set('metHt', [num2str(metHt), '[nm]']);

            model.param().set('goldHt', [num2str(goldHt),'[um]']);
            model.param().set('goldMinR', [num2str(goldMinR),'[um]']);
            model.physics('ewfd').feature('pec2').active(gold);

            model.param().set('funnelCurve', num2str(funnelCurve));

            model.param().set('bezierR', [num2str(bezierR),'[um]']);
            model.param().set('bezierZ', [num2str(bezierZ),'[um]']);
            model.param().set('bezierW', num2str(bezierW));

            model.param().set('dielectricRadius', [num2str(dielectricRadius), '[um]']);
            model.param().set('dielectricHeight', [num2str(funnelDielectricHeight), '[um]']);


            % build geometry and mesh
            model.component('comp1').geom().run();
            model.component('comp1').mesh().run();

            %% display geometry
            % figure(1)
            % mphgeom(model);
            % xlim([0 funnelRBot]*1e-6); ylim([0 funnelHt+1]*1e-6);
            % saveas(gcf, [sweepName, '/geometry.png']);

            %% run study
            model.study('std1').run;
            % mphsave(model,"adiabaticSave.mph");

            %% extracting fields
            um = 1e-6; % 1 um / 1 m

            % coordinates for field/permittivity extraction.
            raFields = study.windowR0*um; % start of r
            rbFields = (study.windowR0 + study.windowWidth)*um; % end of r
            drFields = (study.windowDr)*um; % increment of r

            zaFields = (study.windowZ0)*um; % start of z
            zbFields = (study.windowZ0 + study.windowHeight)*um; % end of z
            dzFields = (study.windowDz)*um; % increment of z

            rArrFields =(raFields:drFields:rbFields);
            zArrFields = (zaFields:dzFields:zbFields);

            [rGridFields,zGridFields] = meshgrid(rArrFields,zArrFields);
            rzCoordsFields = [rGridFields(:),zGridFields(:)]';

            % field extraction
            normE = reshape(mphinterp(model,'ewfd.normE', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));
            Er = reshape(mphinterp(model,'ewfd.Er', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));
            Ephi = reshape(mphinterp(model,'ewfd.Ephi', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));
            Ez = reshape(mphinterp(model,'ewfd.Ez', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));

            Br = reshape(mphinterp(model,'ewfd.Br', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));
            Bphi = reshape(mphinterp(model,'ewfd.Bphi', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));
            Bz = reshape(mphinterp(model,'ewfd.Bz', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));

            % % fill output arrays
            normEC(lam0Index) = {normE};
            ErC(lam0Index) = {Er};
            EphiC(lam0Index) = {Ephi};
            EzC(lam0Index) = {Ez};

            BrC(lam0Index) = {Br};
            BphiC(lam0Index) = {Bphi};
            BzC(lam0Index) = {Bz};

            % permittivity extraction
            epsRR = reshape(mphinterp(model,'material.epsilonr11', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));
            epsPhiPhi = reshape(mphinterp(model,'material.epsilonr22', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));
            epsZZ = reshape(mphinterp(model,'material.epsilonr33', 'coord', rzCoordsFields, 'Dataset','dset1'), size(rGridFields));

            epsRRC(lam0Index) = {epsRR};
            epsPhiPhiC(lam0Index) = {epsPhiPhi};
            epsZZC(lam0Index) = {epsZZ};

            %% purcel and such
            if(study.sim == 1)
                emittedArr(lam0Index) = 0.5*mphint2(model,'real(ewfd.Er*conj(ewfd.ecd1.Jer)+ewfd.Ephi*conj(ewfd.ecd1.Jephi)+ewfd.Ez*conj(ewfd.ecd1.Jez))','surface','selection','geom1_csel13_dom','Dataset','dset1','Intvolume','on');
                emittedBottomArr(lam0Index) = -0.5*mphint2(model,'real(ewfd.Er*conj(ewfd.Hphi)-ewfd.Ephi*conj(ewfd.Hr))','line','selection','geom1_csel16_bnd','Dataset','dset1','Intsurface','on');
                emittedTopArr(lam0Index) = 0.5*mphint2(model,'real(ewfd.Er*conj(ewfd.Hphi)-ewfd.Ephi*conj(ewfd.Hr))','line','selection','geom1_csel17_bnd','Dataset','dset1','Intsurface','on');
                emittedRightArr(lam0Index) = 0.5*mphint2(model,'real(ewfd.Ephi*conj(ewfd.Hz)-ewfd.Ez*conj(ewfd.Hphi))','line','selection','geom1_csel18_bnd','Dataset','dset1','Intsurface','on');
            end

            %% transmission/reflection

            % coordinates for power calculations
            raPow = 0.01*um; % start of r
            drPow = 0.01*um; % increment of r
            rbPow = (geomWd)*um; % end of r

            rArrPow = (raPow:drPow:rbPow);

            zTrans= 0.9*supHt*um; % z value for transmission
            zRefl= -0.9*subHt*um; % z value for reflection

            rzCoordsTrans = [rArrPow',zTrans+ 0*rArrPow']';
            rzCoordsRefl = [rArrPow',zRefl+ 0*rArrPow']';

            % intensity extraction
            IntTrans = mphinterp(model, 'ewfd.Poavz', 'coord', rzCoordsTrans, 'dataset', 'dset1');
            IntRefl = mphinterp(model, 'ewfd.Poavz', 'coord', rzCoordsRefl, 'dataset', 'dset1');

            % power calculations
            powTrans = sum(rArrPow.*IntTrans)*drPow*2*pi; % transmitted power
            powRefl = sum(rArrPow.*IntRefl)*drPow*2*pi; % reflected power

            % % fill output arrays
            powTransArr(lam0Index) = powTrans;
            powReflArr(lam0Index) = powRefl;

            %% Far field

            IntTransFar = mphinterp(model, 'ewfd.normEfar', 'coord', rzCoordsTrans, 'dataset', 'dset1').^2./377;
            powTransFar = sum(zTrans.*rArrPow.*IntTransFar/(rArrPow.^2 + zTrans.^2).^(3/2))*drPow*2*pi;
            powTransFarArr(lam0Index) = powTransFar;
            powReflFarArr(lam0Index) = powTransFar;
            % powTransFarSzArr(lam0Index) = powTransFar;

            theta0 = 0*pi/180;
            theta1 = 90*pi/180;
            dthetaPow = 0.01;
            thetaArr = (theta0:dthetaPow:theta1);
            IntTransFar = mphinterp(model, 'ewfd.normEfar', 'coord', um.*[sin(thetaArr).',cos(thetaArr).'].', 'dataset', 'dset1').^2./377;
            powTransFarArr(lam0Index) = sum(sin(thetaArr).*IntTransFar)*dthetaPow*2*pi;

            IntTransFarArr(lam0Index,:) = IntTransFar;



        end
        sweepTime = toc(tStart);
        disp(['    sweepTime: ', num2str(sweepTime/60), ' min']);

        %% saving
        if(numLayersArrIndex == length(numLayersArr) + 1)
            filename = 'purcelReference';
        else
            filename = makeRunName(study);
            fileNameArr(numLayersArrIndex) = filename;
        end
        save([sweepName, '/', filename, '.mat'],...
            'lam0Arr', ...
            'ErC', 'EphiC', 'EzC', 'BrC', 'BphiC', 'BzC', 'normEC',...
            'epsRRC', 'epsPhiPhiC', 'epsZZC', ...
            'rGridFields', 'zGridFields', ...
            'powTransArr', 'powReflArr', 'powTransFarArr', 'powReflFarArr', 'IntTransFarArr', 'thetaArr', ...
            'emittedArr', 'emittedBottomArr', 'emittedTopArr', 'emittedRightArr',...
            'funnelRBot', 'funnelRTop', 'funnelHt', 'numLayers', 'gold', 'goldHt','subRef',...
            'funnelCurve','funnelType',...
            'dielectricRadius',...
            'supHt', 'subHt', 'pmlSz', 'geomWd', 'meshFunnelRatio', 'meshSmallRatio',...
            'study');

    end

    save([sweepName, '/fileNameArr.mat'], 'fileNameArr', 'study');

end