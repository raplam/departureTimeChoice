% This script runs the simulations corresponding to the Fig. 3.4, 3.5, 3.6
% and 3.8 of my PhD thesis.
% Last modified by Raphael Lamotte, on October 24, 2018.

clear all
close all
path(genpath(cd),path);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

col=get(gca,'colororder');
longcol=parula(10);
screensize = get( groot, 'Screensize' );
figVerif=figure;
set(figVerif,'Position',[0,0,screensize(3)/2,0.4*screensize(4)]);
figCumPlots=figure;
set(figCumPlots,'Position',[0,0,screensize(3)/2,0.8*screensize(4)]);
figDelta=figure;
set(figDelta,'Position',[0,0,screensize(3)/2,0.4*screensize(4)]);
figMetrics=figure;
set(figMetrics,'Position',[0,0,screensize(3)/2,.6*screensize(4)]);

dt=1/60;
departureTimes=-1.5:dt:1.5;
Nt=length(departureTimes);
Capacity=1/2;
congestion=generateBottleneck(Capacity);

settings.maxIter=200;
settings.display='off'; %'on' 'off'
settings.knownEq=0;

revisionProtocol.exponent=1;
revisionProtocol.fun=@(R,U,lambda)SmithRevisionProtocolExponent(R,U,lambda,revisionProtocol.exponent);


% %% Heterogeneity and response rate
du=1.5;
lambda=[1,5,10];
stdStar=[0,0.1,0.2,0.3,0.4,0.5];
% prepare cumplots
lambdaCumPlot=[1;3;2];
stdStarCumPlot=[2;5;5];
dayCumPlot=[127,200;130,200;129,200];
potGainCumPlot=zeros(3,2);
pShiftsCumPlot=zeros(3,2);
%
for indl=1:length(lambda)
    revisionProtocol.rate=lambda(indl)/Nt;
    for inds=1:length(stdStar)
        if stdStar(inds)>0
            N=10;
            tstar=-((1:N)-(1+N)/2)*stdStar(inds)/sqrt((N^2-1)/12);
        else
            tstar=0;
        end
        population=generateSParctan(tstar,(1-du/2)*ones(size(tstar)),(1+du/2)*ones(size(tstar)),ones(size(tstar)),4*ones(size(tstar)));
        [fS,hist]=runIterationsContinuum(departureTimes,settings,congestion,population,revisionProtocol,[]);
        figure(figVerif)
        subplot(2,3,indl) % Potential gain
        plot(1:settings.maxIter,hist.potGain);
        hold on
        xlabel('Days');
        ylabel('Potential gain [\%]');
        title(['$\lambda=',num2str(lambda(indl)),'$']);
        ylim([0,100]);
        subplot(2,3,3+indl)
        semilogy(1:settings.maxIter,hist.shifts);
        hold on
        xlabel('Days');
        ylim([0.001,1])
        ylabel('Proportion of shifts');
        title(['$\lambda=',num2str(lambda(indl)),'$']);
        [~,I]=ismember([indl,inds],[lambdaCumPlot,stdStarCumPlot],'rows');
        if I>0
            figure(figCumPlots)
            % global cumulative I/O
            subplot(4,3,(I-1)*3+1)
            plot(departureTimes,cumsum(sum(squeeze(hist.R(:,:,dayCumPlot(I,1))),1)),'-','Color',col(1,:));
            hold on
            plot(hist.arrivalTimes(dayCumPlot(I,1),:),cumsum(sum(squeeze(hist.R(:,:,dayCumPlot(I,1))),1)),'--','Color',col(1,:));
            plot(departureTimes,cumsum(sum(squeeze(hist.R(:,:,dayCumPlot(I,2))),1)),'-','Color',col(2,:));
            plot(hist.arrivalTimes(dayCumPlot(I,2),:),cumsum(sum(squeeze(hist.R(:,:,dayCumPlot(I,2))),1)),'--','Color',col(2,:));
            ylim([0,1]);
            xlabel('Time');
            ylabel('Cum. Input/Output');
            title(['$\lambda = ', num2str(lambda(lambdaCumPlot(I))), ', \sigma^*=', num2str(stdStar(stdStarCumPlot(I))),', \delta=', num2str(du) '$']);
            legend({['In, day ',num2str(dayCumPlot(I,1))],...
                ['Out, day ',num2str(dayCumPlot(I,1))],...
                ['In, day ',num2str(dayCumPlot(I,2))],...
                ['Out, day ',num2str(dayCumPlot(I,2))]},'location','southeast');
            % cumulative I/O per family day 1
            subplot(4,3,(I-1)*3+2)
            for indf=1:population.N
                plot(departureTimes,squeeze(hist.R(indf,:,dayCumPlot(I,1)))/dt,'-','Color',longcol(indf,:));
                hold on
            end
            xlabel('Time');
            ylabel('Departure rate');
            title(['$\lambda = ', num2str(lambda(lambdaCumPlot(I))), ', \sigma^*=', num2str(stdStar(stdStarCumPlot(I))),'$, day ', num2str(dayCumPlot(I,1)),', $\delta=', num2str(du) '$']);
            % cumulative I/O per family day 2
            subplot(4,3,(I-1)*3+3)
            for indf=1:population.N
                plot(departureTimes,squeeze(hist.R(indf,:,dayCumPlot(I,2)))/dt,'-','Color',longcol(indf,:));
                hold on
            end
            xlabel('Time');
            ylabel('Departure rate');
            title(['$\lambda = ', num2str(lambda(lambdaCumPlot(I))), ', \sigma^*=', num2str(stdStar(stdStarCumPlot(I))),'$, day ', num2str(dayCumPlot(I,2)),', $\delta=', num2str(du) '$']);
            legend({['$t^*=', num2str(tstar(1),2),'$'],...
                ['$t^*=', num2str(tstar(2),2), '$'],['$t^*=', num2str(tstar(3),2), '$']',...
                ['$t^*=', num2str(tstar(4),2), '$'],['$t^*=', num2str(tstar(5),2), '$'],...
                ['$t^*=', num2str(tstar(6),2), '$'],['$t^*=', num2str(tstar(7),2), '$'],...
                ['$t^*=', num2str(tstar(8),2), '$'],['$t^*=', num2str(tstar(9),2), '$'],['$t^*=', num2str(tstar(10),2), '$']},'location','southeast');
            %
            potGainCumPlot(I,:)=hist.potGain(dayCumPlot(I,:));
            pShiftsCumPlot(I,:)=hist.shifts(dayCumPlot(I,:));
        end
    end
end
figure(figVerif)
for i=1:3
    subplot(2,3,lambdaCumPlot(i))
    scatter(dayCumPlot(i,:),potGainCumPlot(i,:),10,col(stdStarCumPlot(i,:),:),'o');
    subplot(2,3,3+lambdaCumPlot(i))
    scatter(dayCumPlot(i,:),pShiftsCumPlot(i,:),10,col(stdStarCumPlot(i,:),:),'o');
end
subplot(2,3,6)
legend({['$\sigma^*=', num2str(stdStar(1)), '$'],...
    ['$\sigma^*=', num2str(stdStar(2)), '$'],...
    ['$\sigma^*=', num2str(stdStar(3)), '$'],...
    ['$\sigma^*=', num2str(stdStar(4)), '$'],....
    ['$\sigma^*=', num2str(stdStar(5)), '$'],....
    ['$\sigma^*=', num2str(stdStar(6)), '$']});

%% What about delta?
lambda=[1,5,10];
stdStar=0.4;
du=[-1.5,-1,-0.5,.5,1,1.5];
% cumPlot
lambdaCumPlot=10;
duCumPlot=-1.5;
dayCumPlot=[100,200];
potGainCumPlot=zeros(1,2);
pShiftsCumPlot=zeros(1,2);

for indl=1:length(lambda)
    for indu=1:length(du)
        revisionProtocol.rate=lambda(indl)/length(departureTimes);
        if stdStar>0
            N=10;
            tstar=0-((1:N)-(1+N)/2)*stdStar/sqrt((N^2-1)/12);
        else
            tstar=0;
        end
        population=generateSParctan(tstar,(1-du(indu)/2)*ones(size(tstar)),(1+du(indu)/2)*ones(size(tstar)),ones(size(tstar)),4*ones(size(tstar)));
        [fS,hist]=runIterationsContinuum(departureTimes,settings,congestion,population,revisionProtocol,[]);
        figure(figDelta)
        subplot(2,3,indl) % Potential gain
        j=find(sort(unique(abs(du)))==abs(du(indu)));
        if du(indu)>0
            plot(1:settings.maxIter,hist.potGain,'-','Color',col(j,:));
        else
            plot(1:settings.maxIter,hist.potGain,'--','Color',col(j,:));
        end
        hold on
        xlabel('Days');
        ylabel('Potential gain [\%]');
        title(['$\lambda=',num2str(lambda(indl)),'$']);
        subplot(2,3,3+indl) % Shifts
        j=find(sort(unique(abs(du)))==abs(du(indu)));
        if du(indu)>0
            semilogy(1:settings.maxIter,hist.shifts,'-','Color',col(j,:));
        else
            semilogy(1:settings.maxIter,hist.shifts,'--','Color',col(j,:));
        end
        ylim([0.001,1])
        hold on
        xlabel('Days');
        ylabel('Proportion of shifts');
        title(['$\lambda=',num2str(lambda(indl)),'$']);
        % CumPlot
        if lambda(indl)==lambdaCumPlot && du(indu)==duCumPlot
            figure(figCumPlots)
            % global cumulative I/O
            subplot(4,3,10)
            plot(departureTimes,cumsum(sum(squeeze(hist.R(:,:,dayCumPlot(1))),1)),'-','Color',col(1,:));
            hold on
            plot(hist.arrivalTimes(dayCumPlot(1),:),cumsum(sum(squeeze(hist.R(:,:,dayCumPlot(1))),1)),'--','Color',col(1,:));
            plot(departureTimes,cumsum(sum(squeeze(hist.R(:,:,dayCumPlot(2))),1)),'-','Color',col(2,:));
            plot(hist.arrivalTimes(dayCumPlot(2),:),cumsum(sum(squeeze(hist.R(:,:,dayCumPlot(2))),1)),'--','Color',col(2,:));
            ylim([0,1]);
            xlabel('Time');
            ylabel('Cum. Input/Output');
            title(['$\lambda = ', num2str(lambdaCumPlot), ', \sigma^*=', num2str(stdStar),', \delta=', num2str(du(indu)) '$']);
            legend({['In, day ',num2str(dayCumPlot(1))],...
                ['Out, day ',num2str(dayCumPlot(1))],...
                ['In, day ',num2str(dayCumPlot(2))],...
                ['Out, day ',num2str(dayCumPlot(2))]},'location','southeast');
            % cumulative I/O per family day 1
            subplot(4,3,11)
            for indf=1:population.N
                plot(departureTimes,squeeze(hist.R(indf,:,dayCumPlot(1)))/dt,'-','Color',longcol(indf,:));
                hold on
            end
            xlabel('Time');
            ylabel('Departure rate');
            title(['$\lambda = ', num2str(lambdaCumPlot), ', \sigma^*=', num2str(stdStar),'$, day ', num2str(dayCumPlot(1)),', $\delta=', num2str(du(indu)) '$']);
            % cumulative I/O per family day 2
            subplot(4,3,12)
            for indf=1:population.N
                plot(departureTimes,squeeze(hist.R(indf,:,dayCumPlot(2)))/dt,'-','Color',longcol(indf,:));
                hold on
            end
            xlabel('Time');
            ylabel('Departure rate');
            title(['$\lambda = ', num2str(lambdaCumPlot), ', \sigma^*=', num2str(stdStar),'$, day ', num2str(dayCumPlot(2)),', $\delta=', num2str(du(indu)) '$']);
            legend({['$t^*=', num2str(tstar(1),2),'$'],...
                ['$t^*=', num2str(tstar(2),2), '$'],['$t^*=', num2str(tstar(3),2), '$']',...
                ['$t^*=', num2str(tstar(4),2), '$'],['$t^*=', num2str(tstar(5),2), '$'],...
                ['$t^*=', num2str(tstar(6),2), '$'],['$t^*=', num2str(tstar(7),2), '$'],...
                ['$t^*=', num2str(tstar(8),2), '$'],['$t^*=', num2str(tstar(9),2), '$'],['$t^*=', num2str(tstar(10),2), '$']},'location','southeast');
            %
            potGainCumPlot(:)=hist.potGain(dayCumPlot);
            pShiftsCumPlot(:)=hist.shifts(dayCumPlot);
        end
    end
end
figure(figDelta)
I1=find(lambda==lambdaCumPlot);
j=find(sort(unique(abs(du)))==abs(duCumPlot));
subplot(2,3,I1)
scatter(dayCumPlot,potGainCumPlot,10,col(j,:),'o');
subplot(2,3,3+I1)
scatter(dayCumPlot,pShiftsCumPlot,10,col(j,:),'o');

subplot(2,3,4)
legend({['$\delta=', num2str(du(1)), '$'],...
    ['$\delta=', num2str(du(2)), '$'],...
    ['$\delta=', num2str(du(3)), '$'],...
    ['$\delta=', num2str(du(4)), '$'],....
    ['$\delta=', num2str(du(5)), '$'],...
    ['$\delta=', num2str(du(6)), '$']});

%% Social cost decomposition
du=[1.5,1.5,-1.5];
stdStar=[0.1,.4,.4];
lambda=1:1:10;
for ind=1:length(stdStar)
    if stdStar(ind)>0
        N=10;
        tstar=0-((1:N)-(1+N)/2)*stdStar(ind)/sqrt((N^2-1)/12);
    else
        tstar=0;
    end
    for indl=1:length(lambda)
        revisionProtocol.rate=lambda(indl)/length(departureTimes);
        population=generateSParctan(tstar,(1-du(ind)/2)*ones(size(tstar)),(1+du(ind)/2)*ones(size(tstar)),ones(size(tstar)),4*ones(size(tstar)));
        [fS,hist]=runIterationsContinuum(departureTimes,settings,congestion,population,revisionProtocol,[]);
        plotSocialCostComponents(hist.R,hist.U,departureTimes,hist.arrivalTimes,figMetrics,ind,1,longcol(indl,:))
    end
    if du(ind)>0
        % Compute Theoretical Equilibrium
        Na=1000; % the number of agents that should be considered for each tstar
        congestion = computeEquilibriumBottleneck(Na,congestion,population);
        eqUtilities=zeros(population.N,Na);
        for i=1:population.N
            eqUtilities(i,:)=(population.UO{i}(congestion.eqDepartures(:,i))+population.UD{i}(congestion.eqArrivals(:,i)))';
        end
        SC=-squeeze(sum(sum(eqUtilities,1),2))/population.N/Na;
        % Plot it
        figure(figMetrics) % scatter plots
        subplot(3,3,3*(ind-1)+3) % SC
        plot([0,100],[SC,SC],'k');
        hold on
        subplot(3,3,3*(ind-1)+2) % Total SP
        plot([0,100],(SC-sum(sum(congestion.eqArrivals-congestion.eqDepartures))/population.N/Na)*[1,1],'k');
        hold on
        subplot(3,3,3*(ind-1)+1) % Total SP
        plot([0,100],sum(sum(congestion.eqArrivals-congestion.eqDepartures))/population.N/Na*[1,1],'k');
        hold on
    end
    % add labels
    for j=3:-1:1
        subplot(length(stdStar),3,(ind-1)*3+j)
        title(['$\sigma^*=',num2str(stdStar(ind)),', \delta=',num2str(du(ind)),'$']);
        xlabel('Potential gain [\%]');
        xlim([0,100]);
        box on
        switch j
            case 1
                ylabel('Average delay');
            case 2
                ylabel('Average schedule penalty');
            case 3
                ylabel('Average congestion cost');
                ylims=ylim;
        end
        ylim([0,ylims(2)]);
    end
end
subplot(3,3,6)
legend({'$\lambda=1$','$\lambda=2$','$\lambda=3$','$\lambda=4$',...
    '$\lambda=5$','$\lambda=6$','$\lambda=7$','$\lambda=8$',...
    '$\lambda=9$','$\lambda=10$'});
