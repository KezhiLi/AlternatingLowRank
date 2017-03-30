%Cramer-Rao Bound for Hankel structured matrices using the
%ABC-parametrization
%normalized = must be multiplied by sigma^2 to get true CRLB
%Martin Sundin 2011-10-21

function crb = crb_hankel_abc(a,C,r,A,n,p)

%Generate S matrix for n times p hankel
k = min(n,p);
S = sparse(n*p,n+p-1);
for i = 1:(k-1)
    S(i:(n-1):(n*(i-1)+1),i) = ones(i,1);
    S((n*(p-i+1)):(n-1):(n*p-i+1),n+p-i) = ones(i,1);
    %for j = 1:i
        %S(i-j+1+n*(j-1),i) = 1;
        %S(n*p-i-n*j+n+1,n+p-i) = 1;
    %end
end

if n >= p
    for i = p:n
        S(i:(n-1):(i+(p-1)*(n-1)),i) = ones(p,1);
%         for j = 1:p
%             S(i+(n-1)*(j-1),i) = 1;
%         end
    end
end
if p > n
    for i = n:p
        %size(n*i:(n-1):(n*(i-1)+1))
        S((n*(i-n+1)):(n-1):(n*(i-1)+1),i) = ones(n,1);
%         for j = 1:n
%             S(n*i+(n-1)*(j-1),i) = 1;
%         end
    end
end

%Compute dg/dalpha = [U V]
%Derivative matrix
D = zeros(r,n+p-3);
Aparam = [-a; eye(r-1) zeros(r-1,1)];
u = Aparam(:,1);
D(:,1) = u;
for i = 2:(n+p-2)
    u = Aparam*u;
    D(:,i) = u;
end

%Check D is finite
if sum(~isfinite(D(:))) ~= 0
    crb = inf; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('CRLB breakdown: infinite D')
    return
end

%Calculate V
V = sparse(n+p-1,r);
%i = 1
V(2,1) = -C(1);
V(3,1) = -C(1)*D(1,1) - C*D(:,1);
for k = 3:(n+p-2)
    s = - C(1)*D(1,k-1) - C*D(:,k-1);
    for j = 1:(k-2)
        s = s - C*D(:,j)*D(1,k-j-1);
    end
    V(k+1,1) = s;
end
%i = 2,3,...,r
for i = 2:r
    for k = i:(n+p-2)
        s = -C(1)*D(i,k-1);
        for j = 1:(k-2)
            s = s - C*D(:,j)*D(i,k-j-1);
        end
        V(k+1,i) = s;
    end
end

E1 = sparse(r,1);
E1(1) = 1;
W = [[E1 D]' V];

T1 = S*W;
T2 = A*T1;

Fisher = T2'*T2;

if det(Fisher) < 1e-6
   crb = inf;
   disp('CRLB Breakdown: Near singular Fisher')
   return
end

crb = trace(Fisher\(T1'*T1));
crb = real(crb);
