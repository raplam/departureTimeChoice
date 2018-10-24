function Probabilities = SmithRevisionProtocolExponent(R,U,lambda,p)
% This function returns the matrix of probabilities of shifting according
% to a discrete adaptation of Smith's revision protocol.
% R should be a vector of size (1,Nt) indicating current departures.
% U should be a vector of size (1,Nt) indicating current utilities.
% lambda and p are parameters of the revision protocol.
% (Probabilities)_{i,j} indicates the probabilities that users currently in
% i switch to j.
% Last modified by Raphael Lamotte, on October 24, 2018.
Nt=length(R);
Probabilities=lambda*max(zeros(Nt,Nt),repmat(U,Nt,1)-repmat(U',1,Nt)).^p;
Probabilities=Probabilities./repmat(max(sum(Probabilities,2),ones(Nt,1)),1,Nt); % makes sure that the sum of probabilities do not exceed 1 for one sending zone.
end

