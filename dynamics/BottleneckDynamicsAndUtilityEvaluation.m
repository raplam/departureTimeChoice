function [Utilities,tin,queue,avgArrivalTimes,cumUsers]=BottleneckDynamicsAndUtilityEvaluation(departureTimes,R,S,population)
%  This function runs the dynamics (i.e. simulates traffic for some peak
% hour) given the demand (i.e. the departure times and the corresponding
% matrix of departure rates) and the bottleneck capacity S.
% Note: the departureTimes are assumed to be uniformly spaced

% We assume that users actually depart continuously over an interval of duration T
% centered around departureTimes(i) and compute the average costs.

% Variables:
% - tin: we split time in slightly more intervals than length(departureTimes).
% The reason is that a queue might disappear within one of these invervals,
% so in this case we want to treat users separately according to whether
% they experienced a queue or not. The times tin are the boundaries between
% these intervals. They are such that users complete trips at a constant rate
% during every interval [tin(i)+queue(i)/S,tin(i+1)+queue(i+1)/S].
% - queue contains the queue at times tin.
% - Utilities containt the utility for each departureTime
% - avgArrivalTime containt the average arrival time for each departureTime
% (for illustration purposes). Note: Utilities are not simply
% UO(avgArrivalTime)+UD(avgArrivalTime). They are actually computed as the
% integral of UO(t)+UD(t) over the distribution of arrival times.
% - cumUsers contains the cumulative number of departures (for illustration
% purposes).
%
% Last modified by Raphael Lamotte, on October 24, 2018.

T=unique(departureTimes(2:end)-departureTimes(1:(end-1)));
if length(T)>1
    if max(T)-min(T)>10^-14 % this is not a numerical error.
        error('Departure Times are assumed to be uniformly spaced');
    else
        T=mean(T);
    end
end
w=sum(R,1);
Nt=length(departureTimes);
tin=zeros(1,2*Nt);% The final size is unknown (it depends on the number of transitions), so the memory allocation is overestimated.
queue=zeros(1,2*Nt);% The final size is unknown (it depends on the number of transitions), so the memory allocation is overestimated.
avgArrivalTimes=zeros(1,Nt);
tin(1)=departureTimes(1)-T/2;

Correspondance_with_departureTimes=zeros(Nt,2*Nt);
Utilities=zeros(population.N,Nt);
k=2; % index within tin.

for i=1:Nt % for all DepartureTimes
    if queue(k-1)+w(i)-S*T>=0 % there is a queue all along the interval
        tin(k)=departureTimes(i)+T/2;
        queue(k)=queue(k-1)+w(i)-S*T;
        Correspondance_with_departureTimes(i,k)=1;
        avgArrivalTimes(i)=departureTimes(i)+(queue(k)+queue(k-1))/(2*S);
    elseif queue(k-1)==0 % there was no queue at all during the interval
        tin(k)=departureTimes(i)+T/2;
        queue(k)=0;
        Correspondance_with_departureTimes(i,k)=1;
        avgArrivalTimes(i)=departureTimes(i);
    else % the queue disappeared during the interval
        duration_w_queue=queue(k-1)/(S-w(i)/T);
        tin(k)=tin(k-1)+duration_w_queue;
        queue(k)=0;
        Correspondance_with_departureTimes(i,k)=duration_w_queue/T;
        k=k+1;
        tin(k)=departureTimes(i)+T/2;
        queue(k)=0;
        Correspondance_with_departureTimes(i,k)=1-duration_w_queue/T;
        avgArrivalTimes(i)=duration_w_queue/T*(tin(k-2)+queue(k-2)/S+tin(k-1))/2+...
            (1-duration_w_queue/T)*(tin(k)+tin(k-1))/2;
    end
    k=k+1;
end
tin=tin(1:k-1);
queue=queue(1:k-1);
tex=tin+queue/S;
for i=1:population.N
    avgUO=(population.intUO{i}(tin(2:end))-population.intUO{i}(tin(1:end-1)))./(tin(2:end)-tin(1:end-1));
    avgUD=zeros(size(avgUO));
    testTexDifference=(tex(2:end)-tex(1:end-1)>eps);
    ind=find(testTexDifference);
    ind2=find(~testTexDifference);
    avgUD(ind)=(population.intUD{i}(tex(ind+1))-population.intUD{i}(tex(ind)))./(tex(ind+1)-tex(ind));
    avgUD(ind2)=population.UD{i}(tex(ind2));
    Utilities(i,:)=Correspondance_with_departureTimes(:,2:k-1)*(avgUO+avgUD)';
end
cumUsers=cumsum(R*Correspondance_with_departureTimes(:,1:k-1),2);

