function C=hankelconstraint(m,n)
% Construct hankel constraint matrix C such that C*vec(Y)=0
% for all m by n hankel matrices Y.
%
% Magnus Jansson 111018

N = m+n-1; 

y=1:N;
Y = hankel(y(1:m),y(m:N));

C=[];
for k=2:N-1
    i = find(Y(:)==k);
    for ii=1:length(i)-1
        v=zeros(1,m*n);
        v(i(ii))=1;v(i(ii+1))=-1;
        C=[C; v];
    end
end


