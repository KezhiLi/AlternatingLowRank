function S=hankelstructure(m,n)
% Construct a matrix S such that vec(Y) = S*theta 
% for an m by n hankel matrix Y where theta is a vector with 
% parameters from the first column and last row of Y.
%
% Magnus Jansson 111101

N = m+n-1; 

y=1:N;
Y = hankel(y(1:m),y(m:N));
S = zeros(m*n,N);
S(1,:) = [1 zeros(1,N-1)];
for k=2:N-1
    i = find(Y(:)==k);
    v = zeros(1,N);
    v(k) = 1;
    S(i,:) = repmat(v,length(i),1);
end
S(m*n,:) = [zeros(1,N-1) 1];
end