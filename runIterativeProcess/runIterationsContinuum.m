function [finalState,history] = runIterationsContinuum(departureTimes,settings,congestion,population,revisionProtocol,R)
% This function takes as input a population of individuals with their
% preferences (population), a congestion mechanism (congestion), the range of allowed departure times (departureTimes),
% an adjustment mechanism (revisionProtocol), an initial matrix of
% departure rates (R) and some settings.
% It then mimicks the day to day evolution.
%
% Outputs:
% - finalState is a struct variable containing informations regarding the
% state on the last iteration of the simulation.
% - history is a struct variable containing informations for all days of
% the simulation.
%
% Last modified by Raphael Lamotte, on October 24, 2018.


col=get(gca,'colororder');
Nt=length(departureTimes);

% memory allocation
history.TU=NaN(1,settings.maxIter);               % history of total utility
history.potGain=NaN(1,settings.maxIter);          % Potential Gain history
history.R=NaN(population.N,Nt,settings.maxIter);  % history of state variable
history.U=NaN(population.N,Nt,settings.maxIter);    % history of utilities
history.shifts=NaN(1,settings.maxIter);          % history of utilities
switch congestion.mechanism
    case 'MFD'
        % preliminary treatment of L
        [uniqueL,~,indfromUniqueL]=unique(population.L);
        Nl=length(uniqueL);
        A=sparse(indfromUniqueL,1:population.N,ones(population.N,1));
        dep4Dyn=reshape(repmat(departureTimes,Nl,1),[],1);
        L4Dyn=reshape(repmat(uniqueL,1,Nt),[],1);
        UtilitiesAtDep=zeros(population.N,Nt); % precomputed to save computational time in case there are many events
        for i=1:population.N
            UtilitiesAtDep(i,:)=population.UO{i}(departureTimes);
        end
        history.arrivalTimes=NaN(settings.maxIter,Nl,Nt);      % history of arrival times
    case 'bottleneck'
        history.arrivalTimes=NaN(settings.maxIter,Nt);      % history of arrival times
end


if isempty(R)
    R=ones(population.N,Nt)/population.N/Nt; %matrix indicating which proportion of each of the N classes depart at each of the Nt times. Initially random (but normalized)
end
for iter=1:settings.maxIter
    %% Update
    history.R(:,:,iter)=R;
    % Run dynamics
    switch congestion.mechanism
        case 'MFD'
            % Create vector of fictive users (one user = all users having
            % same l departing at same time)
            weights=reshape(A*R,[],1);
            % Run dynamics with fictive users
            [d,t,indDep,indExChrono,indExDep,n,v] = MFDDynamics(dep4Dyn,L4Dyn,weights,congestion.speed_fct);
            % Fill matrix of arrival times 
            history.arrivalTimes(iter,:,:)=reshape(t(indExDep),[1,Nl,Nt]);
            % Fill matrix of utilities
            Utilities=zeros(population.N,Nt);
            for i=1:population.N
                Utilities(i,:)=population.UD{i}(squeeze(history.arrivalTimes(iter,indfromUniqueL(i),:))');
            end
            Utilities=Utilities+UtilitiesAtDep;
        case 'bottleneck'
            [Utilities,t,n,history.arrivalTimes(iter,:),cumUsers]=BottleneckDynamicsAndUtilityEvaluation(departureTimes,R,congestion.S,population);
        otherwise
            error('unknown congestion mechanism');
    end
    history.potGain(iter)=100*sum(sum(R.*(repmat(max(Utilities,[],2),1,Nt)-Utilities)./abs(Utilities),2));
    history.TU(iter)=sum(sum(R.*Utilities));
    % Compute new R
    [R,history.shifts(iter)] = updateDeparturesContinuum(R,Utilities,revisionProtocol);
    % Run_script_making_plots
    if strcmp(settings.display,'on')||(strcmp(settings.display,'final')&& iter==settings.maxIter)
        fprintf('iteration %i\n',iter);
        makePlotsContinuum;
        if exist(settings.additionalPlots,'file')==2 % 2 is the code that exist returns when a file with this name exists
            run(settings.additionalPlots);
        end
    end
    % save history
    history.U(:,:,iter)=Utilities;
end
finalState.t=t;
finalState.n=n;
