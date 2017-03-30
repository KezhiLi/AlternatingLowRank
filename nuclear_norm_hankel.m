%Nuclear norm minimization solution to the problem
%Martin Sundin, 

%y = A*h + n, i.e.
%min ||X||_*, X hankel
%A*iS*vec(X) = b
function h_hat = nuclear_norm_hankel(A,y,lambda)

[~,N] = size(A);
n = ceil(N/2);
p = N - n + 1;
iS = hankel_struct_inv(n,p);
AiS = A*iS;


%nuclear norm minimization
cvx_begin sdp quiet
    variable X(n,p) hankel;
    minimize norm_nuc(X) + lambda/2*norm(y-AiS*X(:),2);
cvx_end

h_hat = iS*X(:);

