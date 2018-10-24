function [tollFun,intTollFun] = computeTollHomogeneous(population,duration)
% This script computes the socially optimal toll, when the bottleneck is
% forecasted to be used for a given "duration".
% This preliminary version only works for homogeneous users with symmetric
% schedule preferences (i.e. such that at equilibrium, the bottleneck is
% used over the interval [tstar-duration/2, tstar+duration/2]. We could
% circumvent this limitation in future implementations, but this was not
% needed when written.
% Last modified by Raphael Lamotte, on October 24, 2018.

if length(population.tstar)>1
    error('pricing only works with homogeneous population');
end

t0=population.tstar-duration/2;
tend=population.tstar+duration/2;

Ueq=population.UO{1}(t0)+population.UD{1}(t0);
if abs(Ueq-population.UO{1}(tend)-population.UD{1}(tend))>10^(-10)
    error('Current pricing algorithm only works with symmetric schedule preferences');
end

tollFun=@(x)(x>t0 & x<tend).*(population.UO{1}(x)+population.UD{1}(x)-Ueq);
intTollFun=@(x)(x>t0 & x<tend).*(population.intUO{1}(x)+population.intUD{1}(x)-...
    population.intUO{1}(t0)-population.intUD{1}(t0)-Ueq*(x-t0))+...
    (x>=tend).*(population.intUO{1}(tend)+population.intUD{1}(tend)-...
    population.intUO{1}(t0)-population.intUD{1}(t0)-Ueq*(tend-t0));
end