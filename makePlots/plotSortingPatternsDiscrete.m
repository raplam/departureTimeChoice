%% Context
% Specific plots are called within the runIterationsContinuum and runIterationsDiscrete functions,
% right after the normal plots (makePlotsContinuum and MakePlotsDiscrete) 
%
% This specific code shows how users sort in time depending on their trip
% length.
% Last modified by Raphael Lamotte, on October 24, 2018.

figure(2)
colormap jet
scatter(t(indDep),population.L(I_s2f),6,population.tstar(I_s2f),'o');
hold on
scatter(t(indExDep),population.L(I_s2f),6,population.tstar(I_s2f),'x');
xlabel('Time');
ylabel('Trip length');
legend({'Departures','Arrivals'});
xlim([4,10]);
box on
hold off
