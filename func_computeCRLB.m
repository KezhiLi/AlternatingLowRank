function [ CRLB ] = func_computeCRLB( sigma2, A, X, r )
%% I/O

%% Run
[n_dim,p_dim] = size(X);

%Obtain singular vectors
[U,Lambda,V] = svd(X);

U0 = U(:,1:r);
V0 = V(:,1:r);

U1 = U(:,r+1:end);
V1 = V(:,r+1:end);

%Construct P
%P = [ kron(V1,U0), kron(V0,U0), kron(V0,U1) ];

%CRLB check
if ((n_dim+p_dim)*r - r^2) ~= rank(A* [ kron(V1,U0), kron(V0,U0), kron(V0,U1) ] )
    disp('CRLB breakdown')
    CRLB = inf;
else
    CRLB = sigma2 * trace(inv( [ kron(V1,U0), kron(V0,U0), kron(V0,U1) ]'*(A'*A)*[ kron(V1,U0), kron(V0,U0), kron(V0,U1) ] ));
    CRLB = real(CRLB);
end

end

