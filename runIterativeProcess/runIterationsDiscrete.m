function [finalState,history] = runIterationsDiscrete(settings,congestion,population,updateProportion,dep)
% This function takes as input a population of individuals with their
% preferences (population), a congestion mechanism (congestion), their initial departure times (dep),
% the proportion that update their decision at every iteration and mimics
% the day to day evolution.
%
% Outputs:
% - finalState is a struct variable containing informations regarding the
% state on the last iteration of the simulation.
% - history is a struct variable containing informations for all days of
% the simulation.
% Last modified by Raphael Lamotte, on October 24, 2018.

history.potGain = zeros(1,settings.maxIter);            % Potential Gain history
history.dep = NaN(population.N,settings.maxIter);     % history of departure times
I_s2f = 1:population.N;
I_f2s = 1:population.N;
lastUpdate=zeros(population.N,1); % last iteration the user's decision has been updated
for iter = 1:settings.maxIter
    history.dep(:,iter)=dep;
    %% Sort users according to their departure times
    [sortedDep,I]  =  sort(dep(I_s2f)); % sorting dep(I_f2s) is faster than sorting dep, because dep(I_f2s) should be already almost sorted (except at first iteration)
    Irev(I)=1:population.N;
    I_f2s = Irev(I_f2s); % this is such that dep(I_f2s) = sortedDep; "f2s" because it goes from the "family space" to the "sorted space".
    I_s2f(I_f2s) = 1:population.N;
    
    %% Run the dynamics
    switch congestion.mechanism
        case 'MFD'
            [d,t,indDep,indExChrono,indExDep,n,v] = MFDDynamics(sortedDep,population.L(I_s2f),ones(population.N,1)/population.N,congestion.speed_fct);
        case 'bottleneck'
            % TO DO - with discrete agents
        otherwise
            error('unknown congestion mechanism');
    end
    
    %% Update the decisions of some users
    
    % Choose those who will be updated pseudo-randomly, giving more chances to those that
    % have not been updated for a long time
    indUpdate = selectPseudoRandomly(exp(-lastUpdate/10),ceil(updateProportion*population.N));
    newdeps = zeros(size(indUpdate)); %new departure times
    for j = 1:length(indUpdate)
        switch congestion.mechanism
            case 'MFD'
                % Compute prior utility
                prevDep = t(indDep(I_f2s(indUpdate(j))));
                prevArr = t(indExDep(I_f2s(indUpdate(j))));
                prevU = population.UO{indUpdate(j)}(prevDep)+population.UD{indUpdate(j)}(prevArr);
                % Update utility
                [newdeps(j),bestU] = updateDeparturesDiscreteMFD(population.L(indUpdate(j)),d',t',congestion.vf,population.UO{indUpdate(j)},population.UD{indUpdate(j)},prevDep,prevU);
                history.potGain(iter) = history.potGain(iter)+max((bestU-prevU)/abs(bestU),0)/length(indUpdate);
            case 'bottleneck'
                % TO DO - with discrete agents
        end
    end
    lastUpdate(indUpdate)=iter*ones(length(indUpdate),1);
    dep(indUpdate) = newdeps;
    % Run_script_making_plots
    if strcmp(settings.display,'on')||(strcmp(settings.display,'final')&& iter==settings.maxIter)
        fprintf('iteration %i\n',iter);
        makePlotsDiscrete;
        if exist(settings.additionalPlots,'file')==2 % 2 is the code that exist returns when a file with this name exists
            run(settings.additionalPlots);
        end
    end
end
finalState.d=d;
finalState.t=t;
finalState.n=n;
finalState.v=v;
finalState.indDep=indDep;
finalState.indExChrono=indExChrono;
finalState.indExDep=indExDep;
end

