function [A] = func_generateA(N,P,M,prm_sensing)
%% I/O declaration
% Generate sensing matrix
%Dave Zachariah 2011-10-10

%% Initialize
if prm_sensing == 0
   A = func_matrixcompletion(N,P,M);
   disp('[Generate A: Matrix completion]') 
elseif prm_sensing == 1
   A = func_rip(N,P,M);
   disp('[Generate A: Random matrix RIP]')  
elseif prm_sensing == 2
   A = func_sensehankel(N,P,M);
   disp('[Generate A: Sense Hankel]')  
elseif prm_sensing == 3
   A =  func_sensehermitian(N,P,M);
   disp('[Generate A: Sense Hermitian]')   
else
   A = [];
   disp('[Generate A: ERROR! Unknown sensing structure]')   
end




end


%------------------------------------------------
%------------------------------------------------
function [A] = func_matrixcompletion(N,P,M)

%Allocate memory
A = zeros(M,N*P);

%Randomize location
for m = 1:M
    A_tmp   = zeros(N,P);
    n_coord = round(rand*(N-1) + 1);
    s_coord = round(rand*(P-1) + 1);
    A_tmp(n_coord,s_coord) = 1;
    
    A(m,:) = reshape( A_tmp, N*P, 1 )';
end

end



function [A] = func_rip(N,P,M)

%Allocate memory
A = zeros(M,N*P);

%Random linear combination
A = sqrt(1/M) * randn(M,N*P);

end



function [A] = func_sensehankel(N,P,M)

%Allocate memory
A = zeros(M,N*P);

%Number of samples
T     = N+P-1; 

%Locations of samples in matrix
order        = randperm(T);
t_loc        = zeros(1,T);
t_loc(1:M)   = 1:M;
t_loc(order) = t_loc;
H_loc        = hankel( t_loc(1:N), t_loc(N:N+P-1) );

%TEMP
%disp(t_loc)

%Randomize location
for m = 1:M
    
    %Find observed indices in vec(X)
    idx = find( H_loc == m );
    
    %Linear operator
    a         = zeros(1,N*P);
    a(1,idx) = ones(1,length(idx));
    A(m,:)    = a/sum(a);   
    
end

end



function [A] = func_sensehermitian(N,P,M)

%Allocate memory
A = zeros(M,N*P);

%Randomize location
i_set = [];
j_set = [];
count = 1;



%Randomize location
for m = 1:M
    A_tmp   = zeros(N,P);
    
    i_coord = round(rand*(N-1) + 1);
    j_coord = round(rand*(i_coord-1) + 1);
    
    A_tmp(i_coord,j_coord) = 1/2;
    A_tmp(j_coord,i_coord) = 1/2;
    
    A(m,:) = reshape( A_tmp, N*P, 1 )';
end



% while (length(i_set)<M) && (length(j_set)<M)
%     
%     %Randomize i and j
%     i_coord = round(rand*(N-1) + 1);
%     j_coord = round(rand*(i_coord-1) + 1);
%     
%     %disp([ i_coord  j_coord ])
%     %pause
%     
%     %Check set membership
%     %if (~ismember(i_coord,i_set)) && (~ismember(j_coord,j_set))
%     if count <= M   
%         %Add to sets
%         i_set = union(i_set,i_coord);
%         j_set = union(j_set,j_coord);
%         
%         %Create measurement
%         A_tmp   = zeros(N,P);
%         A_tmp(i_coord,j_coord) = 1/2;
%         A_tmp(j_coord,i_coord) = 1/2;
%         
%         A(count,:) = reshape( A_tmp, N*P, 1 )';
%        
%         
%         count = count + 1;
%     end
%     
% end

end

