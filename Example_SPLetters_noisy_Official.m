clear;
clc
clear all;
  A = [1 1 0 1 1 1 ; 1 1 1 0 0 0 ; 0 1 1 0 0 0; 1 0 0 1 0 0 ; 1 0 0 0 1 1 ; 1 0 0 0 1 1];
 Det = diag(sum(A,2));
 L = Det - A;
 
 E = 1/6;
 Discrete_L = eye(6)-(1/6)*L;
 N=size(Discrete_L,1);
 M=N;
 

 X_0 = [1.3389 2.0227 1.9872 6.0379 2.7219 1.9881]';


 X_k = X_0;
 X_List = X_0;

 End_Value = 30;
  for i = 1:End_Value
     X_k = Discrete_L * X_k;
     X_List = [X_List X_k];
  end
        nn=6;
        r=1;
        i=nn;
        HankelFirstColumn = X_List(r,1:i);
        HankelLastRow = X_List(r,i:(2*i-1));
        Hankels{r,i} = hankel(HankelFirstColumn,HankelLastRow);
         
        T_HankelFirstColumn = X_List(r,(1+1):(i+1));
        T_HankelLastRow = X_List(r,(i+1):((2*i-1)+1));
        T_Hankels{r,i} = hankel(T_HankelFirstColumn,T_HankelLastRow);
       
        SubstructionHankel{r,i} = -Hankels{r,i}+T_Hankels{r,i};
        
 real_sub= SubstructionHankel{r,i};
 

 
        r=1;
        i=6;
        x_struct = 'hank';
        A=eye(nn*nn);
        A(3,4)=1;
        A(4,3)=1;
        A(8,9)=1;
        A(9,8)=1;
        A(13,14)=1;
        A(14,13)=1;
        A(13,19)=1;

AA=A(1,:);
kk=1;
for ii=2:nn^2-1
    if sum(abs(AA(kk,:)-A(ii,:)))~=0
        kk=kk+1;
        AA=[AA;A(ii,:)];
    end
end   


times=100;
tt=8;

result_ALS=zeros(3,tt);
result_D_ALS=zeros(3,tt);


 for j=1:tt;
j
      SMNR=24+j*2;
 sigma2_n_X0 = norm(X_0)^2/N * 10^(-SMNR/10);
 
norm_1=zeros(3,times);
norm_2=zeros(3,times);


count=1;
while (count<times+1);
%count
%for count=1:times; 
X_k = X_0;
 X_List = X_0;
noise       = sqrt(sigma2_n_X0) * randn(6,End_Value+1);
 
  for ii = 1:End_Value
%      noise       = sqrt(sigma2_n_X0) * randn(N,1);
%      X_k = Discrete_L * X_k +noise;
     X_k = Discrete_L * X_k ;
     X_List = [X_List X_k];
  end
 
 X_List=X_List+ noise;
%  xx=X_List(:,end);
%  ckmean=mean(xx);
 

        HankelFirstColumn = X_List(r,1:i);
        HankelLastRow = X_List(r,i:(2*i-1));
        Hankels{r,i} = hankel(HankelFirstColumn,HankelLastRow);
         
        T_HankelFirstColumn = X_List(r,(1+1):(i+1));
        T_HankelLastRow = X_List(r,(i+1):((2*i-1)+1));
        T_Hankels{r,i} = hankel(T_HankelFirstColumn,T_HankelLastRow);
       
        SubstructionHankel{r,i} = -Hankels{r,i}+T_Hankels{r,i};
 

        
        P=nn;
        K=5;
        abs_tol = 1e-6;
        rel_tol = 1.01;
        %y=reshape(SubstructionHankel{r,i},nn*nn,1);
        y=A*reshape(SubstructionHankel{r,i},nn*nn,1);
        yy=AA*reshape(SubstructionHankel{r,i},nn*nn,1);
        %%%%%%%
        sigma2_n_X0_temp=0.000001;
        %%%%%%
        
        tic; tstart = tic;
        [X_hat_ls_struct] = func_LS_lowrankrec_proj2_letters( y, A, sigma2_n_X0_temp, nn,P,K, abs_tol, rel_tol, x_struct );
         %[X_hat_ls_struct] = func_LS_lowrankrec_proj( y, A, sigma2_n_X0_temp, nn,P,K, abs_tol, rel_tol, x_struct );
        telapsed = toc(tstart);
         norm_1(3,count)=telapsed;

         tic; tstart = tic;
        [X_hat_ls_struct_D_ALS] = func_LS_lowrankrec_proj_kezhi_final3_letters( yy, AA, sigma2_n_X0_temp, nn,P,K, abs_tol, rel_tol, x_struct,...
        0.1, 0.1, 0.1);
        telapsed = toc(tstart);
         norm_2(3,count)=telapsed;
    
        S1 = null(X_hat_ls_struct(1:6,1:6));
        S2 = null(X_hat_ls_struct_D_ALS(1:6,1:6));
        %S3 = null(X_hat_ls_nuclear);
        
        norm_1(1,count)=norm(real_sub-X_hat_ls_struct,'fro')^2/norm(real_sub,'fro')^2;
        norm_2(1,count)=norm(real_sub-X_hat_ls_struct_D_ALS,'fro')^2/norm(real_sub,'fro')^2;
        
         norm_1(2,count)=mean(X_List(1,5:10)*S1/sum(S1));
         norm_2(2,count)=mean(X_List(1,5:10)*S2/sum(S2));


        if (norm_1(1,count)>10-j)||(norm_2(1,count)>10-j)||(abs(norm_1(2,count)-2.6828)/2.6828>1)...
                ||(abs(norm_2(2,count)-2.6828)/2.6828>1)
             count=count-1;
        end

count=count+1;

 end
 
 result_ALS(1,j)=mean(norm_1(1,:));
 result_ALS(2,j)=mean(abs(norm_1(2,:)-2.6828)/2.6828);
 result_ALS(3,j)=mean(norm_1(3,:));
 
 result_D_ALS(1,j)=mean(norm_2(1,:));
 result_D_ALS(2,j)=mean(abs(norm_2(2,:)-2.6828)/2.6828);
 result_D_ALS(3,j)=mean(norm_2(3,:));

 end

        
% xx=[26:2:40];
% figure,
% plot(xx,result_ALS(1,:),'bo--');
% hold on
% plot(xx,result_D_ALS(1,:),'ro-');
% plot(xx,result_ALS(2,:),'b*--');
% plot(xx,result_D_ALS(2,:),'r*-');
% grid
% xlabel('SMNR');
% ylabel('NMSE');        
% legend('error: ALS','error: D-ALS','ALS','D-ALS')     
% 
% 
% figure,
% semilogy(xx,result_ALS(1,:),'bo--');
% hold on
% semilogy(xx,result_D_ALS(1,:),'ro-');
% semilogy(xx,result_ALS(2,:),'b*--');
% semilogy(xx,result_D_ALS(2,:),'r*-');
% grid
% xlabel('SMNR');
% ylabel('NMSE');        
% legend('error: ALS','error: D-ALS','ALS','D-ALS')        
        
result_ALS
result_D_ALS

%save all        
%save all 66noisy        
      