function [d,t,ind_dep,chrono_exits,dep_exits,n,v]=MFDDynamics(dep,L,w,speed_fct)
% This function runs the dynamics (i.e. simulates traffic for some peak
% hour) given the demand (i.e. the departure times and trip lengths of all
% users). The way the network reacts to some demand is governed by a
% speed-MFD.

% Inputs:
% dep: vector of length N of departure times of all users (sorted from the earliest to the latest)
% L: vector of trip lengths of all users (the user departing at dep(i) has trip length L(i))
% w: vector of weights of all users (the user departing at dep(i) has
% weight w(i)). The weight represents how much a user contributes to the
% total accumulation. It is useful when considering a continuum of users.
% speed_fct: function giving the speed as a function of the accumulation

% Outputs:
% t: vector of size 2N*1 including the times of all events (sorted from the earliest to the latest) 
% d: vector of size 2N*1. d(i) is the distance traveled at time t(i) by a fictive user who entered at the beginning of times
% and never exits.
% n: vector of size 2N*1 containing the accumulations right AFTER each
% event.
% v: vector of size 2N*1 containing the speed right AFTER each
% event.
% ind_dep: vector of size N. Such that t(ind_dep(i))=dep(i). (i.e. gives
% the position of the departures in the vector t).
% chrono_exits: vector of size N. Such that the user the ith exit occurred at
% time t(chrono_exits(i))
% dep_exits: vector of size N. Such that the user that entered at dep(i)
% exits at t(dep_exits(i));

% Code written by Raphael Lamotte, raphael.lamotte@epfl.ch
% last updated: Dec 12, 2017

N=length(dep); %number of  users
n=zeros(1,2*N);
t=zeros(1,2*N);
d=zeros(1,2*N);
v=zeros(1,2*N);
ind_dep=zeros(1,N);
chrono_exits=zeros(1,N);
dep_exits=zeros(1,N);

%% These variables are global because they are also modified by the function
%add_exit
d_exit=zeros(N+2,1);% Vector that is filled progressively such that the user
% i must exit when d(t)=d_exit(i+1). The size is N+2 such that for any
% real user, one can always access the d_exit of the user before and after
% (this simplifies the code).
ind_exits=zeros(N+2,1);% Vector that is filled progressively such that the user
% exitting when d(t)=d_exit(i+1) is the one whose input data (dep, weights) was at indice ind_exits(i+1). The size is N+2 such that for any
% real user, one can always access the d_exit of the user before and after
% (this simplifies the code).
d_exit(1)=-1;
d_exit(N+2)=inf;
%using next_one and prev_one allows to work with tables of variable size
% without copying huge tables every time. Instead, we just change the value
% of next_one and prev_one.
next_one=zeros(N+2,1);
prev_one=zeros(N+2,1);

%initialisation
% at first, there is only 1 real user (index=2) that is surrounded by two
% fictive users (index =1 and N+2).
next_one(1)=2;
prev_one(N+2)=2;
d_exit(2)=L(1);
ind_exits(2)=1;
prev_one(2)=1;
next_one(2)=N+2;

n(1)=w(1);
v(1)=speed_fct(n(1));
t(1)=dep(1);
ind_dep(1)=1;

j_dep=2;%the next departure will be the j_dep'th
j_ex=1;%the next exit will be the j_ex'th
next_exit=2; %index of the next exit in the vector d_exit.
for j=2:2*N %for all events
    %     if mod(j,100)==0
    %         display(j);
    %     end
    if j_dep<=N %if there are still users that have to depart
        dt1=dep(j_dep)-t(j-1); %time to the next departure
        dt2=(d_exit(next_exit)-d(j-1))/v(j-1); %time to the next arrival
        dt=min(dt1,dt2);
        d(j)=d(j-1)+v(j-1)*dt;
        t(j)=t(j-1)+dt;
        if dt1<dt2 %departure
            ind_dep(j_dep)=j;
            n(j)=n(j-1)+w(j_dep);
            next_exit=add_exit(L(j_dep)+d(j),j_dep-1,next_exit,j_dep);%add one 
            %element in the vector d_exit (requires changing the values of
            % next_one, prev_one and updating d_exit
            j_dep=j_dep+1;
        else %exit
            n(j)=n(j-1)-w(ind_exits(next_exit));
            chrono_exits(j_ex)=j;
            dep_exits(ind_exits(next_exit))=j;
            j_ex=j_ex+1;
            next_exit=next_one(next_exit);
        end
    else %all vehicles have entered => exit
%         if next_exit==0
%             pause(1)
%         end
        dt=(d_exit(next_exit)-d(j-1))/v(j-1);
        d(j)=d(j-1)+v(j-1)*dt;
        t(j)=t(j-1)+dt;
        n(j)=n(j-1)-w(ind_exits(next_exit));
        chrono_exits(j_ex)=j;
        dep_exits(ind_exits(next_exit))=j;
        j_ex=j_ex+1;
        next_exit=next_one(next_exit);
    end
    v(j)=speed_fct(n(j));
end
function next_exit = add_exit(newd,i,next_exit,ind_input)
% i is the number of elements previously in the vectors (not including the
% 1st and last, which are fictive).
% newd is the value of d when that user must exit
% ind is the index in the input data of the corresponding user
% next_exit is the current index of the next exit in the vector of exits

% Code written by Raphael Lamotte, raphael.lamotte@epfl.ch
% last updated: December 4, 2017

if newd<=d_exit(next_exit)
    d_exit(i+2)=newd;
    ind_exits(i+2)=ind_input;
    first_before=prev_one(next_exit);
    first_after=next_exit;
    next_one(first_before)=i+2;
    prev_one(first_after)=i+2;
    next_one(i+2)=first_after;
    prev_one(i+2)=first_before;
    next_exit=i+2;
else
    n_start=round(sqrt(2*i));
    ind_start=2:round(i/n_start):(i+1);
    [~,ind]=min(abs(d_exit(ind_start)-newd));
    %look for the good position
    if d_exit(ind_start(ind))>newd
        first_before=prev_one(ind_start(ind));
        while d_exit(first_before)>newd
            first_before=prev_one(first_before);
        end
        first_after=next_one(first_before);
    else
        first_after=next_one(ind_start(ind));
        while d_exit(first_after)<newd
            first_after=next_one(first_after);
        end
        first_before=prev_one(first_after);
    end
    
    %update the vectors
    d_exit(i+2)=newd;
    ind_exits(i+2)=ind_input;
    next_one(first_before)=i+2;
    prev_one(first_after)=i+2;
    next_one(i+2)=first_after;
    prev_one(i+2)=first_before;
end
end
end