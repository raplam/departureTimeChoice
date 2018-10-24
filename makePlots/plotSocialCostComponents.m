function [] = plotSocialCostComponents(R,U,departureTimes,arrivalTimes,fig1,inds,indp,col)
% This code is used to generate the figure 3.8 of my thesis. It plots the
% social cost and its components. Warning: it considers a value of time of 1.
% Last modified by Raphael Lamotte, on October 24, 2018.


figure(fig1) % scatter plots
SC=-squeeze(sum(sum(R.*U,1),2));
potGain=100*squeeze(sum(sum(R.*(repmat(max(U,[],2),[1,size(U,2),1])-U)./abs(U),1),2));
subplot(3,3,3*(inds-1)+1) % Delay
scatter(potGain,sum(squeeze(sum(R,1)).*(arrivalTimes'-repmat(departureTimes',[1,size(R,3)])),1),2,col(indp,:),'.');
hold on
subplot(3,3,3*(inds-1)+2) % Total SP
scatter(potGain,SC'-sum(squeeze(sum(R,1)).*(arrivalTimes'-repmat(departureTimes',[1,size(R,3)])),1),2,col(indp,:),'.');
hold on
subplot(3,3,3*(inds-1)+3) % SC
scatter(potGain,SC,2,col(indp,:),'.');
hold on


