function population = generateSParctanPriced(tstar,u1,u2,uref,w,tollFun,intTollFun,mode)
% Same as generateSParctan, except that users pay a toll.
% tollFun(t) gives the toll paid when departing (resp. arriving) at t when
% mode='origin' (resp. 'destination')
% intTollFun(t) is the integral of tollFun from some constant to t.
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
switch mode
    case 'origin' % In this case the toll is a function of the time at which users leave their origin.
        for i=1:population.N
            if u2(i)>u1(i)
                population.uO{i}=@(x)uref(i)*ones(size(x));
                population.uD{i}=@(x)(u1(i)+u2(i))/2+(u2(i)-u1(i))/pi*atan(w(i)*(x-offset(i)));
                population.UO{i}=@(x)uref(i)*(x-tstar(i))-tollFun(x);
                population.UD{i}=@(x)-intAtan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))+intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i));
                population.intUO{i}=@(x)uref(i)*(x-tstar(i)).^2/2-intTollFun(x);
                population.intUD{i}=@(x)-int2Atan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))+(x-offset(i))*intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i));
            else
                population.uD{i}=@(x)uref(i)*ones(size(x));
                population.uO{i}=@(x)(u1(i)+u2(i))/2+(u2(i)-u1(i))/pi*atan(w(i)*(x-offset(i)));
                population.UD{i}=@(x)-uref(i)*(x-tstar(i));
                population.UO{i}=@(x)intAtan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-tollFun(x);
                population.intUD{i}=@(x)-uref(i)*(x-tstar(i)).^2/2;
                population.intUO{i}=@(x)int2Atan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-(x-offset(i))*intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-intTollFun(x);
            end
        end
    case 'destination' % In this case the toll is a function of the time at which users leave their destination.
        for i=1:population.N
            if u2(i)>u1(i)
                population.uO{i}=@(x)uref(i)*ones(size(x));
                population.uD{i}=@(x)(u1(i)+u2(i))/2+(u2(i)-u1(i))/pi*atan(w(i)*(x-offset(i)));
                population.UO{i}=@(x)uref(i)*(x-tstar(i));
                population.UD{i}=@(x)-intAtan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))+intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-tollFun(x);
                population.intUO{i}=@(x)uref(i)*(x-tstar(i)).^2/2;
                population.intUD{i}=@(x)-int2Atan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))+(x-offset(i))*intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-intTollFun(x);
            else
                population.uD{i}=@(x)uref(i)*ones(size(x));
                population.uO{i}=@(x)(u1(i)+u2(i))/2+(u2(i)-u1(i))/pi*atan(w(i)*(x-offset(i)));
                population.UD{i}=@(x)-uref(i)*(x-tstar(i))-tollFun(x);
                population.UO{i}=@(x)intAtan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i));
                population.intUD{i}=@(x)-uref(i)*(x-tstar(i)).^2/2-intTollFun(x);
                population.intUO{i}=@(x)int2Atan(x,(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i))-(x-offset(i))*intAtan(tstar(i),(u1(i)+u2(i))/2,u2(i)-u1(i),offset(i),w(i));
            end
        end
    otherwise
        error('unknown mode of pricing');
end
end
