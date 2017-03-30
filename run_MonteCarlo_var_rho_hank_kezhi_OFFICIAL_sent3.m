%Kezhi 2013-03-19
clear all
clc

%% Set parameters
%Monte Carlo
M_runs   = 50; %% M_runs   = 500;
rho_list = [0.1 0.15 0.2 0.25 0.3 0.35]

%NxP matrix with rank K
N        = 80;
P        = 80;
K        = 4; % 80/80/5
x_struct = 'hank';

%Measurement dimensions and noise
SMNR = 10;
flag_random_sensing = 1;


%% Allocate memory
sq_x = zeros(M_runs,length(rho_list));
sq_n = zeros(M_runs,length(rho_list));

sq_e_ls          = zeros(M_runs,length(rho_list));
sq_e_ls_struct   = zeros(M_runs,length(rho_list));
%%%%
% Kezhi
sq_e_ls_M_ALS          = zeros(M_runs,length(rho_list));
sq_e_ls_struct_M_ALS   = zeros(M_runs,length(rho_list));

sq_e_struct_simplehankel   = zeros(M_runs,length(rho_list));

% sq_e_struct_simplehankel   = zeros(M_runs,length(rho_list));
%%%


sq_e_crlb        = zeros(M_runs,length(rho_list));
sq_e_crlb_struct = zeros(M_runs,length(rho_list));

%% Additional information
store_rank_x           = zeros(M_runs,length(rho_list));
store_rank_x_ls        = zeros(M_runs,length(rho_list));
store_rank_x_ls_struct = zeros(M_runs,length(rho_list));

store_minsvalue_x           = zeros(M_runs,length(rho_list));
store_minsvalue_x_ls        = zeros(M_runs,length(rho_list));
store_minsvalue_x_ls_struct = zeros(M_runs,length(rho_list));

store_bias_x_ls        = zeros(M_runs,length(rho_list));
store_bias_x_ls_struct = zeros(M_runs,length(rho_list));

%%%%
% Kezhi
store_rank_x_ls_M_ALS        = zeros(M_runs,length(rho_list));
store_minsvalue_x_ls_M_ALS        = zeros(M_runs,length(rho_list));
store_bias_x_ls_M_ALS        = zeros(M_runs,length(rho_list));

store_rank_x_ls_struct_M_ALS = zeros(M_runs,length(rho_list));
store_minsvalue_x_ls_struct_M_ALS = zeros(M_runs,length(rho_list));
store_bias_x_ls_struct_M_ALS = zeros(M_runs,length(rho_list));

store_rank_x_struct_simplehankel = zeros(M_runs,length(rho_list));
store_minsvalue_x_struct_simplehankel = zeros(M_runs,length(rho_list));
store_bias_x_struct_simplehankel = zeros(M_runs,length(rho_list));

% store_rank_x_struct_simplehankel = zeros(M_runs,length(rho_list));
% store_minsvalue_x_struct_simplehankel = zeros(M_runs,length(rho_list));
% store_bias_x_struct_simplehankel = zeros(M_runs,length(rho_list));
%%%%

%% Algorithm parameters
%Least Squares
abs_tol = 1e-6;
rel_tol = 1.005;

% C=hankelconstraint(N,P);
% constr_tol=0.001; max_iter=100;

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


   
    count = 1;
    for rho = rho_list
        
        %Display
        %-----------------------------------------
        disp('-----------------------------')
        disp('Monte Carlo progress (%)')
        disp( j/M_runs * 100 )
        disp('rho')
        disp( rho )
        disp('-----------------------------')
        
        %Generate
        [X_true,a_prm, c_prm] = func_generatelowrank_hankel( X_true, K, prm_struct, tol_struct );

        
        %Measurement properties
        M = ceil(rho*N*P);
        y = zeros(M,1);
        
        %Sensing operator
        [A] = func_generateA(N,P,M,flag_random_sensing);
        
        %Generate measurements
        %-----------------------------------------
        %Generate measurement
        sigma2_n_tot = norm(X_true,'fro')^2 * 10^(-SMNR/10);
        n_true       = sqrt(sigma2_n_tot/M) * randn(M,1);
        y            = A*reshape(X_true,N*P,1);
        y            = y + n_true;

        %Compute Cramér Rao Bound
        %-----------------------------------------
        CRLB = func_computeCRLB( sigma2_n_tot/M, A, X_true, K );
        CRLB_struct = (sigma2_n_tot/M)* crb_hankel_abc(a_prm,c_prm,K,A,N,P);
        

        
        % Reconstruction
        %-----------------------------------------
        
        %Least-Squares
        %--------------
        if isfinite(CRLB)
            resid_tol = sigma2_n_tot;
        else
            resid_tol = -2e-1;
            disp('Warning: CRLB = -inf')
        end
        
        
        %%%%%%%%%%%%%%%%
        % Kezhi
        if CRLB_struct>1e8
            CRLB_struct=norm(X_true,'fro')^2;
        end
        if abs(CRLB)>1e8
             CRLB =norm(X_true,'fro')^2;
        end
        %%%%%%%%%%%%%%%
        [X_hat_ls] = func_LS_lowrankrec_proj( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, 'none' );
        disp('LS [s]')
        %disp(toc)
        
        %%%%%%%%%%%%%%%%%%
        % Kezhi
        [X_hat_ls_M_ALS] = func_LS_lowrankrec_proj_kezhi_finalH3( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, 'none',...
        0.5, 0.5, 0.5);
       
        %[ X_hat, iter_flag ] = simplehankel( y, A, m,n,r,C,constr_tol,max_iter);
        %%%%%%%%%%%%%%%%%%%
        
        
        %Least-Squares structure
        %--------------
        %tic
        [X_hat_ls_struct] = func_LS_lowrankrec_proj( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, x_struct );
        disp('LS structure [s]')
        %disp(toc)
        
        %%%%%%%%%%%%%%%%%%
        % Kezhi
        [X_hat_ls_struct_M_ALS] = func_LS_lowrankrec_proj_kezhi_finalH4( y, A, sigma2_n_tot, N,P,K, abs_tol, rel_tol, x_struct,...
        0.1, 0.1, 0.1);
    

        %[ X_hat_struct_simplehankel, iter_flag ] = simplehankel( y, A, N,P,K,C,constr_tol,max_iter);
        [ X_hat_struct_simplehankel, iter_flag ] = simplehankel2( y, A, N,P,K);
        %%%%%%%%%%%%%%%%%%%

        %Store additional information
        %-----------------------------------------
        store_rank_x(j,count)           = rank(X_true);
        store_rank_x_ls(j,count)        = rank(X_hat_ls);
        store_rank_x_ls_struct(j,count) = rank(X_hat_ls_struct);
        
          
        store_minsvalue_x(j,count)           = min(svds(X_true, rank(X_true)));
        store_minsvalue_x_ls(j,count)        = min(svds(X_hat_ls, rank(X_hat_ls)));
        store_minsvalue_x_ls_struct(j,count) = min(svds(X_hat_ls_struct, rank(X_hat_ls_struct)));
        
        store_bias_x_ls(j,count)         = sum( X_true(:) - X_hat_ls(:) );
        store_bias_x_ls_struct(j,count)  = sum( X_true(:) - X_hat_ls_struct(:) );
        
        %%%%%%%%
        % Kezhi
        store_rank_x_ls_M_ALS(j,count)        = rank(X_hat_ls_M_ALS);
        store_minsvalue_x_ls_M_ALS(j,count)        = min(svds(X_hat_ls_M_ALS, rank(X_hat_ls_M_ALS)));
        store_bias_x_ls_M_ALS(j,count)         = sum( X_true(:) - X_hat_ls_M_ALS(:) );
        
        store_rank_x_ls_struct_M_ALS(j,count) = rank(X_hat_ls_struct_M_ALS);
        store_minsvalue_x_ls_struct_M_ALS(j,count) = min(svds(X_hat_ls_struct_M_ALS, rank(X_hat_ls_struct_M_ALS)));
        store_bias_x_ls_struct_M_ALS(j,count)  = sum( X_true(:) - X_hat_ls_struct_M_ALS(:) );
        
        store_rank_x_struct_simplehankel(j,count) = rank(X_hat_struct_simplehankel);
        store_minsvalue_x_struct_simplehankel(j,count) = min(svds(X_hat_struct_simplehankel, rank(X_hat_struct_simplehankel)));
        store_bias_x_struct_simplehankel(j,count)  = sum( X_true(:) - X_hat_struct_simplehankel(:) );
        
%         store_rank_x_struct_simplehankel(j,count) = rank(X_hat_struct_simplehankel);
%         store_minsvalue_x_struct_simplehankel(j,count) = min(svds(X_hat_struct_simplehankel, rank(X_hat_struct_simplehankel)));
%         store_bias_x_struct_simplehankel(j,count)  = sum( X_true(:) - X_hat_struct_simplehankel(:) );
        %%%%%%%%%%
        
        
        %Measures
        %-----------------------------------------
        sq_x(j,count) = norm(X_true,'fro')^2;
        sq_n(j,count) = norm(n_true,'fro')^2;
        
        sq_e_ls(j,count)        = norm(X_true-X_hat_ls,'fro')^2;
        sq_e_ls_struct(j,count) = norm(X_true-X_hat_ls_struct,'fro')^2;
        %%%
        % Kezhi
        sq_e_ls_M_ALS(j,count)        = norm(X_true-X_hat_ls_M_ALS,'fro')^2;
        sq_e_ls_struct_M_ALS(j,count) = norm(X_true-X_hat_ls_struct_M_ALS,'fro')^2;
        
        sq_e_struct_simplehankel(j,count) = norm(X_true-X_hat_struct_simplehankel,'fro')^2;
%         sq_e_struct_simplehankel(j,count) = norm(X_true-X_hat_struct_simplehankel,'fro')^2;
        %%%
        
        sq_e_crlb(j,count)         = CRLB;
        sq_e_crlb_struct(j,count)  = CRLB_struct;
        
        %Display
        %-----------------------------------------
        disp('Instantaneous SMNR:')
        disp([ SMNR 10*log10(sq_x(j,count)/sq_n(j,count)) ])
        
        disp('Instantaneous SRER:')
        disp(10*log10(sq_x(j,count)./[sq_e_ls(j,count) sq_e_ls_struct(j,count)] ))
        
        %%%%%%%%%%%
        % Kezhi
        disp('Instantaneous SRER of M-ALS:')
        disp(10*log10(sq_x(j,count)./[sq_e_ls_M_ALS(j,count) sq_e_ls_struct_M_ALS(j,count)] ))
        
        %disp('Instantaneous SRER of SimpleHankel:')
        %disp(10*log10(sq_x(j,count)./[sq_e_ls_M_ALS(j,count) sq_e_struct_simplehankel] ))
        %%%%%%%%%%
        
        disp('CRB:')
        disp(10*log10(sq_x(j,count)./[sq_e_crlb(j,count)  sq_e_crlb_struct(j,count)]))
        
        %Add to counter
        count = count + 1;

       
    end

end
toc
%% Convert into SRER and SMNR

SMNR_emp       = 10*log10(  mean(sq_x,1)./mean(sq_n,1)  )

SRER_ls        = 10*log10(  mean(sq_x,1)./mean(sq_e_ls,1)  );
SRER_ls_struct = 10*log10(  mean(sq_x,1)./mean(sq_e_ls_struct,1)  );
%%%
% Kezhi
SRER_ls_M_ALS        = 10*log10(  mean(sq_x,1)./mean(sq_e_ls_M_ALS,1)  );
SRER_ls_struct_M_ALS = 10*log10(  mean(sq_x,1)./mean(sq_e_ls_struct_M_ALS,1)  );

SRER_struct_simplehankel = 10*log10(  mean(sq_x,1)./mean(sq_e_struct_simplehankel,1)  );

%%%

SRER_crlb        = 10*log10(  mean(sq_x,1)./mean(sq_e_crlb,1)  );
SRER_crlb_struct = 10*log10(  mean(sq_x,1)./mean(sq_e_crlb_struct,1)  );

%% Plot
figure
plot( rho_list, SRER_ls,        'b^--', 'LineWidth', 1.5 ), grid on, hold on, box on
plot( rho_list, SRER_ls_M_ALS,        'ro--', 'LineWidth', 1.5 ),
plot( rho_list, SRER_crlb,      'k--', 'LineWidth', 1.5 )
plot( rho_list, SRER_ls_struct, 'b^-.', 'LineWidth', 1.5 )
plot( rho_list, SRER_ls_struct_M_ALS, 'ro-.', 'LineWidth', 1.5 )
plot( rho_list, SRER_struct_simplehankel, 'm*-.', 'LineWidth', 1.5 )
plot( rho_list, SRER_crlb_struct, 'k-', 'LineWidth', 1.5 )

%%%

%%%

xlabel('\xi'), ylabel('SRER [dB]')
 legend('ALS','D-ALS','CRB', 'ALS-hankel','D-ALS-hankel', 'Simple-Hankel','CRB-hankel')
%legend('ALS','D-ALS','CRB', 'ALS-hankel','D-ALS-hankel', 'CRB-hankel')
%% Save data
%save MC_results_rho_001_hank_kezhi3   N P K x_struct rho M_runs rho_list SMNR_emp SRER_ls SRER_ls_struct SRER_crlb SRER_crlb_struct  sq_x sq_n  sq_e_ls sq_e_ls_struct sq_e_crlb sq_e_crlb_struct



