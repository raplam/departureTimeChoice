function congestion = computeEquilibriumBottleneck(Na,congestion,population)
% This function solves the equivalent linear optimization problem of
%
% Iryo, T., Yoshii, T., 2007. Equivalent optimization problem for finding equilibrium in the
% bottleneck model with departure time choices. In: 4th IMA International Conference
% on Mathematics in Transport. pp. 231–244.
%
% This function uses the following inputs:
%     - a vector tstar of desired arrival times tstar (size 1xNstar),
%     - a vector t of Nt regularly spaced possible arrival times (constructed here),
%     - Na, the number of agents that should be considered for each tstar.
%     - alpha the common value of time
%     - W(t,tstar) the schedule penalty at work.
% It returns two arrays of size (Nt,Nstar), containing the equilibrium
% departure and arrival times of all agents.
% Last modified by Raphael Lamotte, on October 24, 2018.

Nstar=length(population.tstar);
tmin=min(population.tstar)-1/congestion.S; 
tmax=max(population.tstar)+1/congestion.S;
dt=1/(Na*Nstar*congestion.S);
t=(tmin:dt:tmax)';

Nt=length(t);
SP=zeros(Nt,Nstar);
for i=1:Nstar
    SP(:,i)=-population.UD{i}(t)-population.UO{i}(t);
end
SP=reshape(SP,[Nt*Nstar,1]); % SP is a vector of schedule penalties, such that SP(k) contains the schedule penalty
% for an agent of type i (\in {1,...Nstar}) who arrives at t(j), where
% j=mod(k,Nt) and i=1+floor(k/Nt), such that k=Nt*(i-1)+j.
A=repmat(speye(Nt),1,Nstar); %AX<=capacity. ensures that every time is chosen by at most one person.
Aeq=kron(speye(Nstar),ones(1,Nt));% Aeq*X=Na*ones(Nstar,1) (demand constraint)
options=optimoptions('linprog','algorithm','dual-simplex');
[arrivalRates,~,~,~,lambda]=linprog(SP,A,ones(Nt,1),Aeq,Na*ones(Nstar,1),zeros(Nt*Nstar,1),[],options);
arrivalRates=reshape(arrivalRates,[Nt,Nstar]);
departureTimes=t-lambda.ineqlin;
congestion.eqArrivals=zeros(Na,Nstar);
congestion.eqDepartures=zeros(Na,Nstar);
for i=1:Nstar
    congestion.eqArrivals(:,i)=t(arrivalRates(:,i)>0.5);
    congestion.eqDepartures(:,i)=departureTimes(arrivalRates(:,i)>0.5);
end
pause(0.1);
end

