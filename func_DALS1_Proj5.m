function [ X_hat ] = func_DALS1_Proj5( y, A, sigma2_noise, N,P,K, abs_tol, rel_tol, x_struct, rho, rho2, rho3,rho4 )
%% Alternating Least-Squares (ILS) Low-Rank Matrix Reconstruction
%Alternating solution to X = arg min || y - A(X) ||^2_,2 subject to rank(X) = K
%In addition: X can have linear structure: Hankel, Toeplitz, Circulant
% Kezhi 2013-011-30

%Input:
% y             observation Mx1
% A             sensing matrix MxNP
% N,P           dimensions of X
% K             rank of X
% abs_tol       absolute limit on squared error (not in use given sigma2_noise)
% sigma2_noise  measurement noise level
% rel_tol       relative limit on reduction of squared error
% x_struct      'none', 'hank', 'toep', 'circ'

%Output:
% X_hat         matrix of rank K  NxP

%% Initialize

%Prior matrix structure
if sum(x_struct ~= 'hank' ) == 0
    prm_struct = 1;
    disp('[LS: hankel]')
elseif sum(x_struct ~= 'toep' ) == 0
    prm_struct = 2;
    disp('[LS: toeplitz]')
elseif sum(x_struct ~= 'circ' ) == 0
    prm_struct = 3;
    disp('[LS: circulant]')
elseif sum(x_struct ~= 'pdef' ) == 0
    prm_struct = 4;
    disp('[LS: positive semidefinite]')
else
    prm_struct = 0;
    disp('[LS: no structure]')
end

%Obtain dimensions
[M,dummy]  = size(y);

%Initialize structures
X_hat   = zeros(N,P);
L       = zeros(N,K);
R       = zeros(K,P);
I_n     = speye(N,N);
I_p     = speye(P,P);

pp_m = 2;

%Project measurement onto domain of X
Py              = reshape( A'*y, N,P );

%%%%%%
HH = Py;
S = hankelstructure(N,P);
%%%%%%
%Py              = reshape( pinv(A)*y, N,P );
%[U0,Lambda0,V0] = svds( Py, K );

% %Initialize iteration
% L        = U0*sqrt(Lambda0); %orthonormal
% % R = pinv(L) * Py;
% % n_L=norm(L);
% % n_R=norm(R);
% % if n_L>n_R
% %     L=L/sqrt(n_L*n_R);
% %     R=R*sqrt(n_L*n_R);
% % else
% %     R=R/sqrt(n_L*n_R);
% %     L=L*sqrt(n_L*n_R);
% % end

% hhat = (A*hankelstructure(N,P))\y;
% Xhat = hankel(hhat(1:N),hhat(N:end));
% L = Xhat(:,1:K);
%     R = L\Xhat;
%     Xhat = L*R;

if prm_struct == 1
theta = pinv(A*S)*y;% not unique if AS singular
X_hat = reshape(S*theta,N,P); % initial hankel estimate
L = X_hat(:,1:K); %first r columns
%R     = pinv(L) * X_hat;
%L     = X_hat   * pinv(R);
R     = L\X_hat;
L     = X_hat/R;
X_hat = L*R;    
z=R;
if prm_struct == 1
    hhat = hankel_proj(X_hat);
    X_hat = hankel(hhat(1:N),hhat(N:end));
end
else
    Py              = reshape( pinv(A)*y, N,P );
    [U0,Lambda0,V0] = svds( Py, K );
    L        = U0*sqrt(Lambda0); %orthonormal
    z = sqrt(Lambda0) * V0';
     %R = pinv(L) * Py;
end

sq_error = inf;

ti=0;

s=L;


%z=R;
%z = pinv(s) * Py;


%z=zeros(K,P);
t=zeros(N,K);
u=zeros(K,P);
% rho=0.5;
% rho2=0.5;
% rho3=0;
while(true)
    
    ti=ti+1;
    
    %Store residual
    sq_error_old = sq_error;
% 

%%%%%%%%%%%%%%%%%%%%%%%%%
%     if prm_struct == 1
%         X_hat = L*R; 
%         theta = S\X_hat(:);
%         X_hat = reshape(S*theta,N,P); 
%     end
%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute R
    %-----------------------------
    %Least-squares
    B_l = A*kron(I_p,L);
    [L_1,L_2]=size(B_l);
    
    

    u=zeros(K,P);
    for pp=1:pp_m;
    %R   = reshape( pinv(B_l)*y, K, P );
    
    R = inv((B_l')*B_l+rho*eye(L_2))*((B_l')*y+rho*reshape((z-u),K*P,1));
    R = reshape(R, K, P);
    w=R+u*rho2;
    
          %Projection
    if prm_struct ~= 0
        %X_hat = func_project_matrix(s*w, prm_struct);
        %z  = pinv(s) * X_hat;
        
        z = inv(s'*s+rho4*eye(K))*(s'*X_hat+rho4*(R+u)); 

    else
        z = w;
    end 
    
        u=u+(R-z)*rho3;
    end
    
    
    %Compute L
    %-----------------------------
    %Least-squares
    
    B_r = A*kron(R',I_n);
    %L   = reshape( pinv(B_r)*y, N, K );
    [R_1,R_2]=size(B_r);
    
 
           

    t=zeros(N,K);
    for pp=1:pp_m;
    L = inv((B_r')*B_r+rho*eye(R_2))*((B_r')*y+rho*reshape((s-t),N*K,1));
    L = reshape(L, N, K);
    v=L+t*rho2;
    
           %Projection
    if prm_struct ~= 0
        %X_hat = func_project_matrix(v*z, prm_struct);
        %s     = X_hat   * pinv(z);        

       %s = X_hat/z; %s = X_hat/z;
    
        s = (X_hat*z'+rho4*(L+t))*inv(z*z'+rho4*eye(K));
  
    else 
        s = v;
    end
     t=t+(L-s)*rho3;
    end
    

X_hat = s*z;
if prm_struct ~= 0
    hhat = hankel_proj(X_hat);
    X_hat = hankel(hhat(1:N),hhat(N:end));
    %X_hat=func_project_matrix(X_hat, prm_struct);
else
    hhat=[X_hat(:,1);X_hat(N,2:end)'];
end
  %  end
    
    %Compute residual
    %-----------------------------
    resid    = norm(y - A*reshape( X_hat, N*P, 1 ),2)+...
        norm([X_hat(:,1);X_hat(N,2:end)'] - hhat,inf);
    sq_error = resid'*resid;
    disp(sq_error)
    if (sq_error < 0.1*sigma2_noise )  || (sq_error_old/sq_error < rel_tol)
    %if (ti>10)  || (sq_error < 0.1*sigma2_noise)
    %if ti>20
        sq_error = sq_error_old; 
        break; %terminate
    end

   
end

%Produce estimate
%---------------
X_hat = s*z;

% if prm_struct ~= 0
% %    
% %     
% % X_hat =HH;
% X_hat=func_project_matrix(X_hat, prm_struct);

%end




end


function [X_hat] = func_project_matrix(X, prm_struct)

if prm_struct == 1
    X_hat = func_lsqhank(X);
elseif prm_struct == 2
    X_hat = func_lsqtoep(X);
elseif prm_struct == 3
    X_hat = func_lsqcirc(X);
elseif prm_struct == 4
    X_hat = func_lsqposdef(X);
else
    disp('[LS ERROR]: unknown matrix structure')
    X_hat = X;
end

end

%--------------------------------------
% PROJECTIONS
%--------------------------------------
function [hh] = func_lsqhank(hin)
    [m,n]=size(hin);
    hin=fliplr(hin);
    tmp = zeros(m+n-1,1);
    for k=-(m-1):(n-1)
        tmp(k+m) = mean(diag(hin,k));
    end
    tmp=flipud(tmp);
    hh = hankel(tmp(1:m), tmp(m:end));
end
%%%%%%%%%%%%%%%%%%%%%
% function [Xhat] = func_lsqhank(hin)
%  hh=hankel_proj(hin);
%  N=(length(hh)+1)/2;
%  Xhat = hankel(hh(1:N),hh(N:end));
% end
%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------
function [X] = func_lsqtoep(X)

[n,m] = size(X);
k = min(n,m);
x = zeros(m,1);
y = zeros(n,1);

for i = 1:(k-1)
    sy = 0;
    sx = 0;
    for j = 1:i
        sy = sy + X(n-i+j,j);
        sx = sx + X(j,m-i+j);
    end
    y(n+1-i) = sy/i;
    x(m+1-i) = sx/i;
end


if n <= m
    for i = 1:(m-n+1)
        s = 0;
        for j = 1:n
            s = s + X(j,i+j-1);
        end
        x(i) = s/n;
    end
    y(1) = x(1);
end

if m < n
    for i = 1:(n-m+1)
        s = 0;
        for j = 1:m
            s = s + X(i+j-1,j);
        end
        y(i) = s/m;
    end
    x(1) = y(1);
end
X = toeplitz(y,x);

end


%--------------------------------------
function [X] = func_lsqcirc(X)

[n,m] = size(X);
if n ~= m
    disp('[LS ERROR]: Circulant matrix not square!');
end
for i = 1:n
	s = 0;
	for j = 1:(n-i+1)
		s = s + X(i+j-1,j);
	end
	for j = 1:(i-1)
		s = s + X(j,n-i+j+1);
	end
	c = s/n;
	for j = 1:(n-i+1)
		X(i+j-1,j) = c;
	end
	for j = 1:(i-1)
		X(j,n-i+j+1) = c;
	end
end

end

%--------------------------------------
function [X] = func_lsqposdef(X)

[n,m] = size(X);
if n ~= m
    disp('[LS ERROR]: Positive semidefinite matrix!');
end

%Symmetrization
X = (X + X')/2;

%Eigendecomposition
[V,L] = eig(X);

for j = 1:n
    %Ensure real valued
    if imag(L(j,j)) ~= 0
        L(j,j) = real(L(j,j));
    end
    
    %Ensure positive
    if real(L(j,j)) < 0
        L(j,j) = 0;
    end
end

X = V*L*V';


end
