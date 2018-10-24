function [R,shifts] = updateDeparturesContinuum(R,U,revisionProtocol)
% Computes the new departure rate matrix R and the total number of users changing departure (shifts).
% Last modified by Raphael Lamotte, on October 24, 2018.

[N,Nt]=size(R);
shifts=zeros(Nt);
for i=1:N
    F=revisionProtocol.fun(R(i,:),U(i,:),revisionProtocol.rate);% F(i,j) contains the probability for each user of i to change to j
    popShifts=F.*repmat(R(i,:)',1,Nt);
    R(i,:)=R(i,:)-sum(popShifts,2)'+sum(popShifts,1);
    shifts=shifts+popShifts;
end
shifts=sum(sum(shifts));
R=round(R*10^10)/10^10;
R=diag(1./sum(R,2))*R/N;
end

