% This is just a small script to plot the different utility functions and
% ensure that they are correct.
% Last modified by Raphael Lamotte, on October 24, 2018.
clear all
close all
t=6:0.01:10;

tstar=[7.9,8];
u1=[1,5];
uref=[1.5,2];
u2=[3,0.5];
w=[4,4];

population=generateSParctan(tstar,u1,u2,uref,w);

dt=0.001;
t=6:dt:10;
figure
for i=1:2
    subplot(2,2,2*i-1);
    plot(t,population.uO{i}(t),':k',t,population.uD{i}(t),':r',t,population.UO{i}(t),'-k',t,population.UD{i}(t),'-r');
    hold on
    plot(t,cumsum(population.uO{i}(t))*dt+population.UO{i}(t(1)),'--k',t,-cumsum(population.uD{i}(t))*dt+population.UD{i}(t(1)),'--r');
    xlabel('time');
    ylabel('u and U');
    legend({'uO','uD','UO','UD','cum(uO)','cum(uD)'});
    subplot(2,2,2*i);
    plot(t,cumsum(population.UO{i}(t))*dt+population.intUO{i}(t(1)),':k',t,cumsum(population.UD{i}(t))*dt+population.intUD{i}(t(1)),':r');
    hold on
    plot(t,population.intUO{i}(t),'-k',t,population.intUD{i}(t),'-r');
    legend({'cum(UO)','cum(UD)','intUO','intUD'});
    xlabel('time');
    ylabel('U and intU');
end