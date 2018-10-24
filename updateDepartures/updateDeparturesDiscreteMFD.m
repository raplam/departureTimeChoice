function [bestDep, bestU, varargout] = updateDeparturesDiscreteMFD(l,d,t,vf,UO,UD,currentDep,currentU)
% function [bestDep,bestU,allTestDep,allTestU] = updateDeparturesDiscreteMFD(l,d,t,vf,UO,UD,currentDep,currentU)
% This function takes as input the two vectors d ant t characterizing the
% dynamics, the characteristics of one given individual (trip length l and
% scheduling preferences UO and UD), her current departure time and utility
% (currentDep, currentU) and the free-flow speed (vf) and tries to solve the global optimization
% problem of finding the best departure time given the current conditions
% with a simple grid search approach.
%
% Outputs:
% bestDep: new optimal departure time
% bestU: utility obtained with the new departure time
% (optional) varargout{1}: vector of evaluated utility values for this individual
% (optional) varargout{2}: vector of times at which the utility has been evaluated for this individual
%
% Last modified by Raphael Lamotte, on October 24, 2018.
if ~ismember(nargout,[2,4])
    error('The expected number of outputs is either 2 or 4')
end

bestDep=currentDep;
bestU=currentU;

if nargout==4
    varargout{1}=[];
    varargout{2}=[];
end

tmin=t(1)-l/vf; % earliest time evaluated in the first round
tmax=t(end); % latest time evaluated in the first round
dt=max(l/vf/20,(tmax-tmin)/10^3); % step for the first grid search

for iround=1:2 % the number of rounds could be increased, but 2 was deemed sufficient.
    testDep=tmin:dt:tmax;
    testDep=testDep(UO(testDep)+UD(testDep)>bestU);
    testEx=NaN(size(testDep));
    prev_d=NaN; % Index of largest instance of t that is smaller than the tested departure. Kept in memory to avoid screening the entire vector t every time.
    prev_a=NaN; % Index of largest instance of t that is smaller than the arrival time corresponding to the tested departure.
    
    for i=1:length(testDep)
        ti=testDep(i);
        if isnan(prev_d)
            if ti>t(1) %first time, starts not too early
                prev_d=find(t>=ti,1,'first');
                dti=d(prev_d-1)+(ti-t(prev_d-1))/(t(prev_d)-t(prev_d-1))*(d(prev_d)-d(prev_d-1));
            else %starts too early
                dti=d(1)-vf*(t(1)-ti);
            end
        else
            prev_d=prev_d-1+find(t(prev_d:end)>=ti,1,'first');
            if prev_d==1
                dti=d(1)-vf*(t(1)-ti);
            else
                dti=d(prev_d-1)+(ti-t(prev_d-1))/(t(prev_d)-t(prev_d-1))*(d(prev_d)-d(prev_d-1));
            end
        end
        if isnan(prev_a)%first time or finishes too early or too late
            if dti+l>d(1)
                if dti+l>=d(end) %finishes too late
                    testEx(i)=t(end)+(dti+l-d(end))/vf;
                else %finishes in time
                    prev_a=find(d>=dti+l,1,'first');
                    testEx(i)=t(prev_a-1)+(l+dti-d(prev_a-1))/(d(prev_a)-d(prev_a-1))*(t(prev_a)-t(prev_a-1));
                end
            else %finishes too early
                testEx(i)=ti+l/vf;
            end
        else
            tmp=find(d(prev_a:end)>=dti+l,1,'first');
            if ~isempty(tmp)
                prev_a=prev_a-1+tmp;
                if prev_a==1
                    testEx(i)=t(prev_a)-(d(1)-dti-l)/vf;
                else
                    testEx(i)=t(prev_a-1)+(l+dti-d(prev_a-1))/(d(prev_a)-d(prev_a-1))*(t(prev_a)-t(prev_a-1));
                end
            else
                prev_a=length(d);
                testEx(i)=t(end)+(dti+l-d(end))/vf;
            end
        end
    end
    testU=UO(testDep)+UD(testEx);
    [bestUtmp,I]=max(testU);
    if bestUtmp>bestU
        bestU=bestUtmp;
        bestDep=testDep(I);
    end
    % prepare for next round
    tmin=bestDep-dt;
    tmax=bestDep+dt;
    dt=dt/100; % new time step
    if nargout==4
        varargout{1}=[varargout{1},testDep];
        varargout{2}=[varargout{2},testU];
    end
end