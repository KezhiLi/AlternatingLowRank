function  [X, L] = func_generatelowrank_pdef( X, r, prm_struct, tol_struct )
%% I/O


%% Construct low rank
[N,dummy] = size(X);

L = randn(N,r);
X = L*L';

end