function  [X, a, c] = func_generatelowrank_hankel( X, r, prm_struct, tol_struct )
%% I/O


%% Construct low rank

%SVD truncation
[U,Sigma,V] = svds(X,r);
X0 = U*Sigma*V';


%{A,B,C} parameterization
[n,p] = size(X0);
x = zeros(n+p-1,1);
x(1:n) = X0(:,1);
x((n+1):(n+p-1)) = X0(n,2:p);

%Parameter fitting
[b,a] = prony(x,r-1,r);
%Compute A,B and C
A = [-a(2:(r+1)); eye(r-1) zeros(r-1,1)];
C = b;
B = [1;zeros(r-1,1)];
%Output parameters
c = b;
a = a(2:(r+1));

%Compute impulse respone of new system
x = zeros(n+p-1,1);
x(1) = C*B;
for i = 2:(n+p-1)
    x(i) = C*A^(i-1)*B;
end

X = hankel(x(1:n),x(n:(n+p-1)));


end