function [ ind] = selectPseudoRandomly(x,n)
% Returns the indices of n values. The probability a value is picked is
% proportional to its value. Values of x are assumed to be positive.
% Last modified by Raphael Lamotte, on October 24, 2018.

S=sum(x);
i=S*sort(rand(n,1));
S_tmp=x(1);
ind_tmp=1;
ind=zeros(n,1);
for j=1:n
    while S_tmp<i(j)
        ind_tmp=ind_tmp+1;
        S_tmp=S_tmp+x(ind_tmp);
    end
    ind(j)=ind_tmp;
end
keep_ind=zeros(n,1);
for j=1:n
    if x(ind(j))~=0
        if j>1
            if ind(j)~=ind(j-1)
                keep_ind(j)=1;
            end
        else
            keep_ind(j)=1;
        end
    end
end
ind=ind(keep_ind>0);        
end

