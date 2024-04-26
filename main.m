clear, clc;
close all

%%%%%---Initialization---%%%%%

%%% Wavelength Sweep
study.lam0Arr = (12:0.5:18);

%%% Simulation Type
study.sim = 0; % controls simulation type, 0 is incident from below, 1 is emitter on top

%sim = 0 Port
study.m = -1; % azimuthal mode number
study.port.Ex = 1;
study.port.Ey = 1i;
study.port.Ez = 0;

%sim = 1 Emitter
study.emitter.activeHt = 200; % [nm]
study.emitter.passiveHt = 300; % [nm]
study.emitter.polarization = 'X'; % X is circular, Z is Z

%%% Curve Stuff
study.curve = 0;   % controls funnel curve, 0 is straight line, 1 is bezier curve

%curve = 0 linear

%curve = 1 bezier
study.bezR = 0.7;%0.7;
study.bezZ = 1.5;%1.5;r
study.bezW = 1;

%%% Overall Funnel geometry
study.RBot = 3;%4.6;%2.85;%2.9;%4.6;
study.RTop = 0.30/2;%0.280/2;%0.33/2;
study.numLayersArr = [10,15,20,25,30]; % Number of bilayers
% study.numLayersArr = [10,12,14,16,18,20,22,24,26,28,30]; % Number of bilayers

%%% Materials
study.funnelType = 0; % 0 for hyperbolic, 1 for hybrid, 2 for dielectric
study.lamPlasma = 12.5; % plasma wavelength to be used by loadEpsDrude if drude == true
study.gamma = 0.68e13;
study.dieEps = 10.23;
study.metEps0 = 12.15;
study.metHt = 80; % metal layer thickness [nm]
study.dieHt = 80; % dielectric layer thickness [nm]
study.subRef = 3.8;

%funnelType = 0 hyperbolic

%funnelType = 1 hybrid
study.dielectricRadius = 1.5;

%funnelType = 2 dielectric

study.lossRatio = 1; % fraction of loss used

%gold coverage
study.gold = 1; % if 0 there will be no gold
study.goldMinR = 0.95*study.RBot; % no gold below the minimum radius

%%% Saving 

% Field saving region (um)
study.windowZ0 = -0.1;
study.windowHeight = 6;
study.windowNz = 100;
study.windowDz = (study.windowHeight - study.windowZ0)/(study.windowNz-1);
study.windowR0 = 0.010;
study.windowWidth = 4.09;
study.windowNr = 200;
study.windowDr = (study.windowWidth - study.windowR0)/(study.windowNr-1);

% Path where results will be saved
study.sweepName = 'results/HMM-sweep';

%%%%%---Run COMSOL---%%%%%
model = mphload("Funnel.mph");
runStudy(study, model);

%% %%%%%---Process and Plot---%%%%%
if(study.sim == 0)
    processStudy(study.sweepName, 50, lam0=14);
else
    processEmission(study.sweepName,lam0=14)
end







