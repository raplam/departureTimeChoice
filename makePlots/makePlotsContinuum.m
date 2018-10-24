%% Context
% This function is called by runIterationsContinuum.m if settings.display is
% either 'on', or 'final'. It displays some indices that are rather
% standard.
% Last modified by Raphael Lamotte, on October 24, 2018.

figure(1)
% Dynamics
subplot(311)
plot(t,n);
xlabel('Time');
ylabel('Accumulation');
if strcmp(congestion.mechanism,'MFD')
    hold on
    plot([t(1),t(end)],[congestion.ncr,congestion.ncr],'--k');
end
xlim([t(1),t(end)]);
hold off

subplot(312)
switch congestion.mechanism
    case 'MFD'
        plot(t,v)
        xlabel('Time');
        ylabel('Speed');
    case 'bottleneck'
        plot(t,sum(cumUsers,1));
        hold on
        plot(t+n/congestion.S,sum(cumUsers,1));
        hold off
        xlabel('Time');
        ylabel('Cumulative number of users');
        legend({'Departures','Arrivals'},'location','northwest');
        axis([t(1),t(end),0,1]);
end

% Convergence
subplot(313)
semilogy(1:iter,history.potGain(1:iter));
xlabel('Iterations');
ylabel('Potential gain');
pause(0.01)
