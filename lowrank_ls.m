%Simple algorithm for recovering a low-rank solution
%Martin Sundin, 2014-03-26

function Xhat = lowrank_ls(y,A,p,q,r)

Xhat = A\y;
Xhat = reshape(Xhat,p,q);
[u,s,v] = svds(Xhat,r);
Xhat = u*s*v';
