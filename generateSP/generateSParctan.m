function population = generateSParctan(tstar,u1,u2,uref,w)
% Generates a population with an arctan and a constant marginal utility
% rate. The arctan utility rate is given by
% u(x)=(u1+u2)/2+(u2-u1)/pi*atan(w*(x-offset(tstar,u1,u2,uref))), such that it tends towards u1
% when x -> -inf and to u2 when x -> +Inf.
% If u1<u2, the arctan marginal utility rate correspond to the utility at
% destination.
% If u1>=u2, the arctan marginal utility rate correspond to the utility at
% the origin.
% The other marginal utility rate is constant equal to uref. uref should be
% between u1 and u2 for the problem to be well-defined (otherwise people
% have no reason to travel).
% offset(tstar,u1,u2,uref) is defined such that the two marginal utility rates
% intersect for x=tstar.
% The inputs tstar, u1, u2, uref and w should all have the same number of elements.
% Last modified by Raphael Lamotte, on October 24, 2018.

population.tstar=reshape(tstar,1,[]);
population.N=length(population.tstar);
if length(unique([numel(tstar),numel(u1),numel(u2),numel(uref),numel(w)]))~=1
    error('All inputs should have the same number of elements');
else
    u1=reshape(u1,1,[]);
    u2=reshape(u2,1,[]);
    uref=reshape(uref,1,[]);
    w=reshape(w,1,[]);
end
offset=population.tstar+tan(pi/2*(u1+u2-2*uref)./(u2-u1))./w;
intAtan=@(x,umean,d,offset,w)umean*(x-offset)+d/pi*((x-offset).*atan(w*(x-offset))-log(w^2*(x-offset).^2+1)/(2*w));
int2Atan=@(x,umean,d,offset,w)umean/2*(x-offset).^2+...
    d/pi*((w^2*(x-offset).^2-1).*atan(w*(x-offset))/(2*w^2)-...
    ((x-offset).*log(w^2*(x-offset).^2+1)-(x-offset))/(2*w));
for i=1:population.N
    if u2(i)>u1(i)
        population.uO{i}=@(x)uref(i)*ones(size(x)); % marginal utility rate at origin
        population.uD{i}=@(x)(u1(i)+u2(i))/2+(u2(i)-u1(i))/pi*atan(w(i)*(x-offset(i))); % marginal utility rate at destination
        population.UO{i}=@(x)uref(i)*(x-tstar(i)); % integral of marginal utility rate at origin
        population.UD{i}=@(x)-intAtan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))+intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i)); % integral of marginal utility rate at destination, with x as lower limit
        population.intUO{i}=@(x)uref(i)*(x-tstar(i)).^2/2; % integral of UO (used to compute average UO over an interval)
        population.intUD{i}=@(x)-int2Atan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))+(x-offset(i))*intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i)); % integral of UD (used to compute average UD over an interval)
    else
        population.uD{i}=@(x)uref(i)*ones(size(x));
        population.uO{i}=@(x)(u1(i)+u2(i))/2+(u2(i)-u1(i))/pi*atan(w(i)*(x-offset(i)));
        population.UD{i}=@(x)-uref(i)*(x-tstar(i));
        population.UO{i}=@(x)intAtan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i));
        population.intUD{i}=@(x)-uref(i)*(x-tstar(i)).^2/2;
        population.intUO{i}=@(x)int2Atan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-(x-offset(i))*intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i));
    end
end
end
