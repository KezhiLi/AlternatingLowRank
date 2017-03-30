%Different version of simplehankel for faster implementation
%Martin Sundin, 2014-03-19

function [Xhat,L,R] = simplehankel2(y,A,p,q,r)

%Initialize
hhat = (A*hankelstructure(p,q))\y;

%Break conditions
maxiter = 100;
iter = 0;
tol1 = 1e-3;
diff = 1;
while (iter < maxiter) && (diff > tol1)
    iter = iter + 1;
    
    %Project to lowrank
    Xhat = hankel(hhat(1:p),hhat(p:end));
    L = Xhat(:,1:r);
    R = L\Xhat;
    Xhat = L*R;
    
    %Project to hankel
    hhat = hankel_proj(Xhat);
    diff = norm([Xhat(:,1);Xhat(p,2:end)'] - hhat,inf);
end