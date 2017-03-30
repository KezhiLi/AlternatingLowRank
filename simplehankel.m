function  [ X_hat, iter_flag ] = simplehankel( y, A, m,n,r,C,constr_tol,max_iter)
% Low rank reconstruction with linear constraints: 
%  y = A*vec(X) + noise 
% X is assumed to be of size m x n and rank r. Ideally X should 
% satisfy the constraints C*vec(X)=0. 
%
% 
% Stopping criteria:   
% constr_tol    max( abs( C*vec(X_hat) ) < constr_tol 
% max_iter      maximum number of iterations. If exceeded, iter_flag=1.
%
%
%
% Magnus Jansson 111107, 111111 (completely new version)


% init
% could be any linear structure here
S = hankelstructure(m,n); %vec(hankel)=S*theta
%AS = A*S;

theta = pinv(A*S)*y;% not unique if AS singular
X_hat = reshape(S*theta,m,n); % initial hankel estimate
L = X_hat(:,1:r); %first r columns
%R     = pinv(L) * X_hat;
%L     = X_hat   * pinv(R);
R     = L\X_hat;
L     = X_hat/R;
X_hat = L*R;    

iter=0; 
iter_flag=0;
while true
       cc = max(abs(C*X_hat(:)));
  if cc < constr_tol, break, end
 iter = iter+1; 
 if iter > max_iter, iter_flag = 1; break, end
% constraints not satisfied
        theta = S\X_hat(:);
        X_hat = reshape(S*theta,m,n); 
%        R     = pinv(L) * X_hat;
%        L     = X_hat   * pinv(R);
        R     = L\ X_hat;
        L     = X_hat/R;
        X_hat = L*R; 
end

end