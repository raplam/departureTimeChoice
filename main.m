%% Presentation
% This script launches a simulation for a predefined number of days
% (iterations). The user can specify hereafter the configuration.
% Last modified by Raphael Lamotte, on October 24, 2018.

clear all
close all
path(genpath(cd),path);

% Define congestion mechanism
Mechanism='MFD';% 'MFD' or 'bottleneck'
Capacity=0.8;

% Define mode (agent-based simulation ('discrete') or continuum-based
% ('continuum')
Mode='continuum';% 'discrete' or 'continuum'

% Define the population of users
% population=generateSPabg();%
tstar=rand(10,1)-0.5;
population=generateSParctan(tstar,0.5*ones(size(tstar)),2.5*ones(size(tstar)),ones(size(tstar)),4*ones(size(tstar)));%

% Define some technical settings
settings.maxIter=200;
settings.display='on'; %'on' 'off' 'final'
% If additional things should be plotted at every iteration, a new script
% should be created and its name should be saved in
% settings.additionalPlots. Otherwise, the default value is 'off'.
settings.additionalPlots='off';

% Define some mechanism-specific parameters
switch Mechanism
    case 'MFD'
        population.L=1+2*rand(population.N,1);
        vf=30;
    case 'bottleneck'
        settings.knownEq=0; % should only set to a positive value if Iryo's method apply (morning commute with identical alpha).
end

% Define some mode-specific parameters
switch Mode
    case 'continuum'
        departureTimes=-1.5:1/60:1.5;
        ini=ones(population.N,length(departureTimes))/population.N/length(departureTimes); % initial distribution of departures
        revisionProtocol.exponent=1;
        revisionProtocol.fun=@(R,U,lambda)SmithRevisionProtocolExponent(R,U,lambda,revisionProtocol.exponent);
        revisionProtocol.rate=1/length(departureTimes);
    case 'discrete'
        updateProportion=0.05;
        ini=population.tstar; % initial departure times
end

%% DO NOT MODIFY AFTER THIS LINE
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

if strcmp(Mechanism,'bottleneck')&&strcmp(Mode,'discrete')
    error('The combination bottleneck+discrete is the only one that is not supported yet. Sorry.');
end
switch Mechanism
    case 'MFD'
        congestion=generateQuadraticSpeedMFD(Capacity,population.L,vf);
    case 'bottleneck'
        congestion=generateBottleneck(Capacity);
        if settings.knownEq>0
            Na=ceil(1000/population.N);
            congestion=computeEquilibriumBottleneck(Na,congestion,population);
        end
end
switch Mode
    case 'continuum'
        [fS,hist]=runIterationsContinuum(departureTimes,settings,congestion,population,revisionProtocol,ini);
    case 'discrete'
        [fS,hist]=runIterationsDiscrete(settings,congestion,population,updateProportion,ini);
end

