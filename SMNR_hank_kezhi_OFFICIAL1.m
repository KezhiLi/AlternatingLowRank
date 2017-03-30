%Kezhi 2014-04-01
clear all
clc

%% Set parameters
%Monte Carlo
M_runs    = 50;
SMNR_list = 5:2.5:20

%NxP matrix with rank K
N        = 100;
P        = 100;
K        = 4;  %% 5, or 6
x_struct = 'hank';

%Measurement dimensions and noise
rho  = 0.2;
flag_random_sensing = 1;


%% Allocate memory
sq_x = zeros(M_runs,length(SMNR_list));
sq_n = zeros(M_runs,length(SMNR_list));

sq_e_ls          = zeros(M_runs,length(SMNR_list));
sq_e_ls_struct   = zeros(M_runs,length(SMNR_list));

sq_e_crlb        = zeros(M_runs,length(SMNR_list));
sq_e_crlb_struct = zeros(M_runs,length(SMNR_list));

%%%%
% Kezhi
sq_e_ls_M_ALS          = zeros(M_runs,length(SMNR_list));
sq_e_ls_struct_M_ALS   = zeros(M_runs,length(SMNR_list));

sq_e_struct_simplehankel   = zeros(M_runs,length(SMNR_list));
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
    [A] = func_generateA(N,P,M,flag_random_sensing);
    
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

        %Compute Cramér Rao Bound
        %-----------------------------------------
        CRLB        = func_computeCRLB( sigma2_n_tot/M, A, X_true, K );
        CRLB_struct = (sigma2_n_tot/M) * crb_hankel_abc(a_prm,c_prm,K,A,N,P);
        
        
        % Reconstruction
        %-----------------------------------------
        
        %Least-Squares
        %--------------
        if isfinite(CRLB)
            resid_tol = sigma2_n_tot;
        else
            resid_tol = 1e-1;
            disp('Warning: CRLB = -inf')
        end
        
%                 %%%%%%%%%%%%%%%%
%         % Kezhi
%         if CRLB_struct>1e8
%             CRLB_struct=norm(X_true,'fro')^2;
%         end
%         if abs(CRLB)>1e8
%              CRLB =norm(X_true,'fro')^2;
%         end
%         %%%%%%%%%%%%%%%
        
        
        [X_hat_ls] = func_LS_lowrankrec_proj( y, A, resid_tol, N,P,K, abs_tol, rel_tol, 'none' );
        disp('LS [s]')
        %disp(toc)
        
                %%%%%%%%%%%%%%%%%%
        % Kezhi
        %[X_hat_ls_M_ALS] = func_LS_lowrankrec_proj_kezhi_finalH3( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, 'none',...
        %0.5, 0.5, 0.5);
        [X_hat_ls_M_ALS] = func_DALS1_Proj5( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, 'none',...
        0.5, 0.5, 0.5,0);
    
        [ X_hat_simple_normal ] = lowrank_ls( y, A, N,P,K);
        %%%%%%%%%%%%%%%%%%%
        
        %Least-Squares structure
        %--------------
        %tic
        [X_hat_ls_struct] = func_LS_lowrankrec_proj( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, x_struct );
        disp('LS structure [s]')
        %disp(toc)

                %%%%%%%%%%%%%%%%%%
        % Kezhi
        [X_hat_ls_struct_M_ALS] = func_DALS1_Proj51( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, x_struct,...
        0.1, 0.1, 0.1,0.01);
    
        %[ X_hat_struct_simplehankel, iter_flag ] = simplehankel( y, A, N,P,K,C,constr_tol,max_iter);
       [ X_hat_struct_simplehankel, iter_flag ] = simplehankel2( y, A, N,P,K);
        %%%%%%%%%%%%%%%%%%%
        
        %Insert other algorithms
        %--------------

        %Measures
        %-----------------------------------------
        sq_x(j,count) = norm(X_true,'fro')^2;
        sq_n(j,count) = norm(n_true,'fro')^2;
        
        sq_e_ls(j,count)        = norm(X_true-X_hat_ls,'fro')^2;
        sq_e_ls_struct(j,count) = norm(X_true-X_hat_ls_struct,'fro')^2;

        sq_e_crlb(j,count)        = CRLB;
        sq_e_crlb_struct(j,count) = CRLB_struct;
        %%%
        % Kezhi
        sq_e_ls_M_ALS(j,count)        = norm(X_true-X_hat_ls_M_ALS,'fro')^2;
        sq_e_ls_struct_M_ALS(j,count) = norm(X_true-X_hat_ls_struct_M_ALS,'fro')^2;
        
        sq_e_struct_simplehankel(j,count) = norm(X_true-X_hat_struct_simplehankel,'fro')^2;
        %sq_e_simple_normal(j,count) = norm(X_true-X_hat_simple_normal,'fro')^2;

        %%%
        
        %Display
        %-----------------------------------------
        disp('Instantaneous SMNR:')
        disp([ SMNR 10*log10(sq_x(j,count)/sq_n(j,count)) ])
        
        disp('Instantaneous SRER:')
        disp(10*log10(sq_x(j,count)./[sq_e_ls(j,count) sq_e_ls_struct(j,count)]))
 
                %%%%%%%%%%%
        % Kezhi
        disp('Instantaneous SRER of M-ALS:')
        disp(10*log10(sq_x(j,count)./[sq_e_ls_M_ALS(j,count) sq_e_ls_struct_M_ALS(j,count)] ))
        %%%%%%%%%%
        
        disp('CRB:')
        disp(10*log10(sq_x(j,count)./[sq_e_crlb(j,count) sq_e_crlb_struct(j,count)]))
        
        %Add to counter
        count = count + 1;

       
    end

end
toc
%% Convert into SRER and SMNR

SMNR_emp         = 10*log10(  mean(sq_x,1)./mean(sq_n,1)  );

SRER_ls          = 10*log10(  mean(sq_x,1)./mean(sq_e_ls,1)  );
SRER_ls_struct   = 10*log10(  mean(sq_x,1)./mean(sq_e_ls_struct,1)  );

%%%
% Kezhi
SRER_ls_M_ALS        = 10*log10(  mean(sq_x,1)./mean(sq_e_ls_M_ALS,1)  );
SRER_ls_struct_M_ALS = 10*log10(  mean(sq_x,1)./mean(sq_e_ls_struct_M_ALS,1)  );

SRER_struct_simplehankel = 10*log10(  mean(sq_x,1)./mean(sq_e_struct_simplehankel,1)  );
%SRER_simple_normal = 10*log10(  mean(sq_x,1)./mean(sq_e_simple_normal,1)  );

%%%

SRER_crlb        = 10*log10(  mean(sq_x,1)./mean(sq_e_crlb,1)  );
SRER_crlb_struct = 10*log10(  mean(sq_x,1)./mean(sq_e_crlb_struct,1)  );

%% Plot
% figure(1)
% plot( SMNR_emp, SRER_ls,        'bo--', 'LineWidth', 1.5 ), grid on, hold on, box on
% plot( SMNR_emp, SRER_ls_struct, 'r^-.', 'LineWidth', 1.5 )
% plot( SMNR_emp, SRER_crlb,      'k', 'LineWidth', 1.5 )
% plot( SMNR_emp, SRER_crlb_struct,      'k--', 'LineWidth', 1.5 )
% 
% xlabel('SMNR [dB]'), ylabel('SRER [dB]')
% legend('ALS','ALS-hankel','CRB','CRB-hankel')
figure
plot( SMNR_emp, SRER_ls,        'b^--', 'LineWidth', 1.5 ), grid on, hold on, box on
plot( SMNR_emp, SRER_ls_M_ALS,        'ro--', 'LineWidth', 1.5 ),
%plot( SMNR_emp, SRER_simple_normal,        'mo--', 'LineWidth', 1.5 ),
plot( SMNR_emp, SRER_crlb,      'k--', 'LineWidth', 1.5 )
plot( SMNR_emp, SRER_ls_struct, 'b^-.', 'LineWidth', 1.5 )
plot( SMNR_emp, SRER_ls_struct_M_ALS, 'ro-.', 'LineWidth', 1.5 )
plot( SMNR_emp, SRER_struct_simplehankel, 'm*-.', 'LineWidth', 1.5 )
plot( SMNR_emp, SRER_crlb_struct, 'k-', 'LineWidth', 1.5 )
axis([min(SMNR_emp) max(SMNR_emp)  -5 inf]); 
xlabel('SMNR [dB]'), ylabel('SRER [dB]')
legend('ALS','ADLS','CRB','ALS-hankel','ADLS-hankel', 'ALE-hankel','CRB-hankel')
%legend('ALS','D-ALS','simple normal','CRB','ALS-hankel','D-ALS-hankel', 'simplehankel','CRB-hankel')
%legend('ALS','D-ALS','CRB','ALS-hankel','D-ALS-hankel', 'Simple-Hankel','CRB-hankel')
 %axis([5 20  -10 40]);

%% Save data
%save MC_results_SMNR_002_hank   N P K x_struct rho M_runs SMNR_list SMNR_emp   SRER_ls SRER_ls_struct SRER_crlb SRER_crlb_struct  sq_x sq_n  sq_e_ls sq_e_ls_struct sq_e_crlb sq_e_crlb_struct
