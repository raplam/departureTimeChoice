% This script launches the simulation used in Fig 3.7 of my PhD thesis.
% Last modified by Raphael Lamotte, on October 24, 2018.

clear all
close all
path(genpath(cd),path);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
LineStyleOrder={'--','-.','-.','--','-'};

col=get(gca,'colororder');
longcol=parula(50);

screensize = get( groot, 'Screensize' );
figPricing=figure;
set(figPricing,'Position',[0,0,screensize(3)*0.5,0.2*screensize(4)]);

dt=1/60;
departureTimes=-1.5:dt:1.5;
Nt=length(departureTimes);
S=1/2;
congestion = generateBottleneck(S);
settings.maxIter=200;
settings.display='off';
settings.knownEq=0;
settings.additionalPlots='off';
revisionProtocol.exponent=1;
revisionProtocol.fun=@(R,s,lambda)SmithRevisionProtocolExponent(R,s,lambda,revisionProtocol.exponent);
settingsShort=settings;
settingsShort.maxIter=5;
settingsShort.display='on';

pricingWindowRatio=[0.85,1,1.15];
lambda=[1,5,10];
    
% Homogeneous idealization
tstar=0;
uref=1;
du=1;
u1=uref-du/2;
u2=uref+du/2;
w=4;

% True population
stdStar=0;
noise=0;
if stdStar>0
    N=10;
    perturbedTstar=tstar-((1:N)-(1+N)/2)*stdStar/sqrt((N^2-1)/12);
else
    N=1;
    perturbedTstar=tstar;
end
perturb1=u1*(1+noise*(rand(size(perturbedTstar))-0.5));
perturb2=u2*(1+noise*(rand(size(perturbedTstar))-0.5));
perturbRef=uref*(1+noise*(rand(size(perturbedTstar))-0.5));
perturbw=w*(1+noise*(rand(size(perturbedTstar))-0.5));

% Preliminary run to obtain various starting points
truePopulationWOpricing=generateSParctan(perturbedTstar,perturb1,perturb2,perturbRef,perturbw);
revisionProtocol.rate=10/length(departureTimes);
[~,histWOpricing]=runIterationsContinuum(departureTimes,settingsShort,congestion,truePopulationWOpricing,revisionProtocol,[]);

figure(figPricing)
for indw=1:length(pricingWindowRatio)
    subplot(1,length(pricingWindowRatio),indw)
    populationWOpricingHomogeneous=generateSParctan(perturbedTstar,u1*ones(size(perturbedTstar)),u2*ones(size(perturbedTstar)),uref*ones(size(perturbedTstar)),w*ones(size(perturbedTstar)));
    [tollFun,intTollFun] = computeTollHomogeneous(populationWOpricingHomogeneous,pricingWindowRatio(indw)*1/S);
 
    populationPricedDestination=generateSParctanPriced(perturbedTstar,perturb1,perturb2,perturbRef,perturbw,tollFun,intTollFun,'destination');
    populationPricedOrigin=generateSParctanPriced(perturbedTstar,perturb1,perturb2,perturbRef,perturbw,tollFun,intTollFun,'origin');
    
    % Simulations
    for indlambda=1:length(lambda)
        revisionProtocol.rate=lambda(indlambda)/length(departureTimes);
        [~,hist]=runIterationsContinuum(departureTimes,settings,congestion,populationPricedDestination,revisionProtocol,squeeze(histWOpricing.R(:,:,1)));
        plot(1:settings.maxIter,hist.potGain,'Color',col(indlambda,:));
        hold on
        [~,hist]=runIterationsContinuum(departureTimes,settings,congestion,populationPricedOrigin,revisionProtocol,squeeze(histWOpricing.R(:,:,1)));
        plot(1:settings.maxIter,hist.potGain,'--','Color',0.5*col(indlambda,:)+0.5*[0.5,0.5,0.5]);
    end
    xlabel('Days');
    ylabel('Potential gain [\%]');
    title(['$\delta=', num2str(du), '$, $\sigma^*=0$, pricing duration: $', num2str(pricingWindowRatio(indw),3), '/s$']);
    
    for i=2:settingsShort.maxIter
        perturb1=1+noise*(rand(size(perturbedTstar))-0.5);
        perturb2=1+noise*(rand(size(perturbedTstar))-0.5);
        perturbRef=1+noise*(rand(size(perturbedTstar))-0.5);
        perturbw=1+noise*(rand(size(perturbedTstar))-0.5);
        populationPricedDestination=generateSParctanPriced(perturbedTstar,u1*perturb1,u2*perturb2,uref*perturbRef,w*perturbw,tollFun,intTollFun,'destination');
        populationPricedOrigin=generateSParctanPriced(perturbedTstar,u1*perturb1,u2*perturb2,uref*perturbRef,w*perturbw,tollFun,intTollFun,'origin');
        for indlambda=1:length(lambda)
            revisionProtocol.rate=lambda(indlambda)/length(departureTimes);
            [~,hist]=runIterationsContinuum(departureTimes,settings,congestion,populationPricedDestination,revisionProtocol,squeeze(histWOpricing.R(:,:,i)));
            plot(1:settings.maxIter,hist.potGain,'Color',col(indlambda,:));
            [~,hist]=runIterationsContinuum(departureTimes,settings,congestion,populationPricedOrigin,revisionProtocol,squeeze(histWOpricing.R(:,:,i)));
            plot(1:settings.maxIter,hist.potGain,'--','Color',0.5*col(indlambda,:)+0.5*[0.5,0.5,0.5]);            
        end
        pause(0.1);
    end
end
subplot(1,length(pricingWindowRatio),length(pricingWindowRatio))
legend({['$\lambda=',num2str(lambda(1)),'$, $\$_d(t)$'],...
    ['$\lambda=',num2str(lambda(1)),'$, $\$_o(t)$'],...
    ['$\lambda=',num2str(lambda(2)),'$, $\$_d(t)$'],...
    ['$\lambda=',num2str(lambda(2)),'$, $\$_o(t)$'],...
    ['$\lambda=',num2str(lambda(3)),'$, $\$_d(t)$'],...
    ['$\lambda=',num2str(lambda(3)),'$, $\$_o(t)$']});
