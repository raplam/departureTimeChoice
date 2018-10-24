%% Context
% This function is called by runIterationsDiscrete.m if settings.display is
% either 'on', or 'final'. It displays some indices that are rather
% standard.
% Last modified by Raphael Lamotte, on October 24, 2018.

figure(1)
% Time series of accumulation (and speed if MFD)
switch congestion.mechanism
    case 'MFD'
        subplot(322)
        plot(t,v);
        xlabel('Time');
        ylabel('Speed');
        subplot(321)
        plot([t(1),t(end)],[congestion.ncr,congestion.ncr],'--k');
        hold on
    case 'bottleneck'
        subplot(311)
end
plot(t,n);
xlabel('Time');
ylabel('Accumulation');
xlim([t(1),t(end)]);
hold off

% Cumulative input/output diagrams
subplot(312)
switch congestion.mechanism
    case 'MFD'
        plot(t(indDep),(1:population.N)/population.N,'-k');
        hold on
        plot(t(indExChrono),(1:population.N)/population.N,'--r');
    case 'bottleneck'
        % TO DO
end
axis([t(1),t(end),0,1]);
hold off
xlabel('Time');
ylabel('Cumulative number');
legend({'Departures','Arrivals'});

% Convergence
subplot(313)
semilogy(1:iter,history.potGain(1:iter));
xlabel('Iterations');
ylabel('Potential gain');

