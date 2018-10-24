% This script launches the simulation used in Fig 6.1 of my PhD thesis.
% Note: I think the computation time increases approximately with the
% square of the number of users. For this specific configuration,
% computations typically last for a few minutes.
% Last modified by Raphael Lamotte, on October 24, 2018.

clear all
close all
path(genpath(cd),path);

% Define congestion mechanism
Mechanism='MFD';
Capacity=0.25;

% Define mode (agent-based simulation ('discrete') or continuum-based
% ('continuum')
Mode='discrete';

% Define the population of users
tstar=repmat(8.05:0.1:9,1,400);
population=generateSParctan(tstar,0.5*ones(size(tstar)),2.5*ones(size(tstar)),ones(size(tstar)),4*ones(size(tstar)));%

% Define some technical settings
settings.maxIter=200;
settings.display='final';

population.L=1*rand(1,population.N);
vf=1;
updateProportion=0.05;
ini=population.tstar-population.L/vf; % initial departure times
settings.additionalPlots='plotSortingPatternsDiscrete';

%% DO NOT MODIFY AFTER THIS LINE
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

congestion=generateQuadraticSpeedMFD(Capacity,population.L,vf);
[fS,hist]=runIterationsDiscrete(settings,congestion,population,updateProportion,ini);

