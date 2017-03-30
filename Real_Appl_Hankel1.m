%Kezhi 2014-04-01
clear all
clc

%% Set parameters
%Monte Carlo
M_runs    = 10;
SMNR_list = 5:2.5:20;

%NxP matrix with rank K
N        = 100;
P        = 100;
K        = 4;  %% 5, or 6
x_struct = 'hank';

%Measurement dimensions and noise
rho  = 0.2;
flag_random_sensing = 1;


%% Allocate memory

sq_e_ls_struct_M_ALS   = zeros(M_runs,length(SMNR_list));


%sq_e_simple_normal   = zeros(M_runs,length(SMNR_list));
%%%

%% Algorithm parameters
%Least Squares
abs_tol = 1e-6;
rel_tol = 1.01;

%C=hankelconstraint(N,P);
%constr_tol=0.001; max_iter=100;
%% Run Monte Carlo
tic
for j = 1:M_runs

    %Generate low-rank matrix
    %-----------------------------------------
    X_true     = randn(N,P); %initial random
    tol_struct = 1e-14; %tolerance for structured matrix

    %Prior matrix structure
    if sum(x_struct ~= 'hank' ) == 0
        prm_struct = 1;
    elseif sum(x_struct ~= 'toep' ) == 0
        prm_struct = 2;
    elseif sum(x_struct ~= 'circ' ) == 0
        prm_struct = 3;
    else
        prm_struct = 0;
    end

    %Generate   
    [X_true,a_prm, c_prm] = func_generatelowrank_hankel( X_true, K, prm_struct, tol_struct );

    %Measurement properties
    M = ceil(rho*N*P);
    y = zeros(M,1);
    
    %Sensing operator
    %[A] = func_generateA(N,P,M,flag_random_sensing);
    
    A_ori = zeros(M, P);
    A = kron( A_ori, eye(N));
    
    
    
    
    
    %SMNR in [dB]
    count = 1;
    for SMNR = SMNR_list
        
        %Display
        %-----------------------------------------
        disp('-----------------------------')
        disp('Monte Carlo progress (%)')
        disp( j/M_runs * 100 )
        disp('SMNR [dB]')
        disp( SMNR )
        disp('-----------------------------')
        
        %Generate measurements
        %-----------------------------------------
        %Generate measurement
        sigma2_n_tot = norm(X_true,'fro')^2 * 10^(-SMNR/10);
        n_true       = sqrt(sigma2_n_tot/M) * randn(M,1);
        y            = A*reshape(X_true,N*P,1);
        y            = y + n_true;



        %%%%%%%%%%%%%%%%%%%

        % Kezhi
        [X_hat_ls_struct_M_ALS] = func_DALS1_Proj51( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, x_struct,...
        0.1, 0.1, 0.1,0.01);

        %%%%%%%%%%%%%%%%%%%
        
        %Insert other algorithms
        %--------------

        %Measures
        sq_e_ls_struct_M_ALS(j,count) = norm(X_true-X_hat_ls_struct_M_ALS,'fro')^2;

        %%%
        
        %Display
%         % Kezhi
         disp('Instantaneous SRER of M-ALS:')
%         disp(10*log10(sq_x(j,count)./[sq_e_ls_M_ALS(j,count) sq_e_ls_struct_M_ALS(j,count)] ))
%         %%%%%%%%%%

        
        %Add to counter
        count = count + 1;

       
    end

end
toc
%% Convert into SRER and SMNR
%%%
% Kezhi
SRER_ls_struct_M_ALS = 10*log10(  mean(sq_x,1)./mean(sq_e_ls_struct_M_ALS,1)  );

%% Plot
figure
plot( SMNR_emp, SRER_ls_struct_M_ALS, 'ro-.', 'LineWidth', 1.5 )
axis([min(SMNR_emp) max(SMNR_emp)  -5 inf]); 
xlabel('SMNR [dB]'), ylabel('SRER [dB]')
legend('ADLS-hankel')
%legend('ALS','D-ALS','simple normal','CRB','ALS-hankel','D-ALS-hankel', 'simplehankel','CRB-hankel')
%legend('ALS','D-ALS','CRB','ALS-hankel','D-ALS-hankel', 'Simple-Hankel','CRB-hankel')
 %axis([5 20  -10 40]);

%% Save data
%save MC_results_SMNR_002_hank   N P K x_struct rho M_runs SMNR_list SMNR_emp   SRER_ls SRER_ls_struct SRER_crlb SRER_crlb_struct  sq_x sq_n  sq_e_ls sq_e_ls_struct sq_e_crlb sq_e_crlb_struct
