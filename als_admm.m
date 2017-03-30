%ALS-ADMM for fixed beta
%Martin Sundin, 2013-02-24

function [Xhat,L,R] = als_admm(y,A,p,q,r,beta,gamma)
%Gamma should be 0 < gamma < 1?

Xhat = A\y;
Xhat = reshape(Xhat,p,q);
[u,s,v] = svds(Xhat,r);
L = u*sqrt(s);
R = sqrt(s)*v';
%gamma = 0.9;

C = hankelconstraint(p,q);
[ll,~] = size(C);
lambda = zeros(ll,1);

AACC = (A'*A) + beta*(C'*C);
Aty = A'*y;

maxiter = 100;
mindiff = 1e-3;
diff = 1;
iter = 0;
while (iter < maxiter) && (diff > mindiff)
    iter = iter + 1
    Xold = Xhat;
    AtyCl = Aty + C'*lambda;
    L = (kron(R,eye(p,p))*AACC*kron(R',eye(p,p)))\(kron(R,eye(p,p))*AtyCl);
    L = reshape(L,p,r);
    R = (kron(eye(q,q),L')*AACC*kron(eye(q,q),L))\(kron(eye(q,q),L')*AtyCl);
    R = reshape(R,r,q);
    Xhat = L*R;
    lambda = lambda + gamma*beta*C*Xhat(:);
    
    diff = norm(Xold - Xhat,'fro');
    if iter>100
        break;
    end
end
