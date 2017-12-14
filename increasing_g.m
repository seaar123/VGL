%-----------------------------------------------------------------
ssss=1;
%-----------------------------------------------------------------
if ssss==1
    filenammm1=['Parameters_05_X' num2str(2) '.mat'];   %3
    load(filenammm1);

    y_mean12(1,:)=Av_g_x1;% Average cost (x1 g0) iterations \all runs
    y_mean12(2,:)=Av_R_lambda_x1;% Average target of cost (x1 g0) iterations \all runs
    
    y_mean135(1,:)=Av_g_x2;% Average cost (x2 g0) iterations \all runs
    y_mean135(2,:)=Av_R_lambda_x2;% Average target of cost (x2 g0) iterations \all runs
    
    
    %SN : single critic netwrok

    
    filenammm3=['Parameters_01_X' num2str(4) '.mat'];    %1   4
    load(filenammm3);

    y_mean12_SN(1,:)=Av_g_x1;% Average cost (x1 g0) iterations \all runs
    y_mean12_SN(2,:)=Av_R_lambda_x1;% Average target of cost (x1 g0) iterations \all runs
    
    y_mean135_SN(1,:)=Av_g_x2;% Average cost (x2 g0) iterations \all runs
    y_mean135_SN(2,:)=Av_R_lambda_x2;% Average target of cost (x2 g0) iterations \all runs
end
%----------------------------------------------------------------------
% Convert the NN cost value to increamental cost value

costt0_g_hat_x1(1)=0;
costt0_g_lambda_x1(1)=0;

costt0_g_hat_x2(1)=0;
costt0_g_lambda_x2(1)=0;

costt0_g_hat_x1_SN(1)=0;
costt0_g_lambda_x1_SN(1)=0;

costt0_g_hat_x2_SN(1)=0;
costt0_g_lambda_x2_SN(1)=0;


JJ_star_g_hat_x1(1:maximum_trail)=y_mean12(1,end);
JJ_star_g_lambda_x1(1:maximum_trail)=y_mean12(2,end);

JJ_star_g_hat_x2(1:maximum_trail)=y_mean135(1,end);
JJ_star_g_lambda_x2(1:maximum_trail)=y_mean135(2,end);


JJ_star_g_hat_x1_SN(1:maximum_trail)=y_mean12_SN(1,end);
JJ_star_g_lambda_x1_SN(1:maximum_trail)=y_mean12_SN(2,end);

JJ_star_g_hat_x2_SN(1:maximum_trail)=y_mean135_SN(1,end);
JJ_star_g_lambda_x2_SN(1:maximum_trail)=y_mean135_SN(2,end);

for i=2:maximum_trail
 %x1
    e_g_hat_x1(i)=abs(abs(JJ_star_g_hat_x1(i))-abs(y_mean12(1,i))); 
    e_g_lambda_x1(i)=abs(abs(JJ_star_g_lambda_x1(i))-abs(y_mean12(2,i))); 
   
    
    costt0_g_hat_x1(i)=costt0_g_hat_x1(i-1)+e_g_hat_x1(i);
    costt0_g_lambda_x1(i)=costt0_g_lambda_x1(i-1)+e_g_lambda_x1(i);

    
    
        e_g_hat_x1_SN(i)=abs(abs(JJ_star_g_hat_x1_SN(i))-abs(y_mean12_SN(1,i))); 
    e_g_lambda_x1_SN(i)=abs(abs(JJ_star_g_lambda_x1_SN(i))-abs(y_mean12_SN(2,i))); 
   
    
    costt0_g_hat_x1_SN(i)=costt0_g_hat_x1_SN(i-1)+e_g_hat_x1_SN(i);
    costt0_g_lambda_x1_SN(i)=costt0_g_lambda_x1_SN(i-1)+e_g_lambda_x1_SN(i);
    
    
 %x2
 
 
    e_g_hat_x2(i)=abs(abs( JJ_star_g_hat_x2(i))-abs(y_mean135(1,i))); %vlx
    e_star_g_lambda_x2(i)=abs(abs(JJ_star_g_lambda_x2(i))-abs(y_mean135(2,i)));%glx1

    
    costt0_g_hat_x2(i)=costt0_g_hat_x2(i-1)+e_g_hat_x2(i);
    costt0_g_lambda_x2(i)=costt0_g_lambda_x2(i-1)+e_star_g_lambda_x2(i);

    
        e_g_hat_x2_SN(i)=abs(abs( JJ_star_g_hat_x2_SN(i))-abs(y_mean135_SN(1,i))); %vlx
    e_star_g_lambda_x2_SN(i)=abs(abs(JJ_star_g_lambda_x2_SN(i))-abs(y_mean135_SN(2,i)));%glx1

    
    costt0_g_hat_x2_SN(i)=costt0_g_hat_x2_SN(i-1)+e_g_hat_x2_SN(i);
    costt0_g_lambda_x2_SN(i)=costt0_g_lambda_x2_SN(i-1)+e_star_g_lambda_x2_SN(i);
    
end
%x1
scallling0_x1=JJ_star_g_hat_x1(end)/costt0_g_hat_x1(end);
costt0_g_hat_x1=scallling0_x1*costt0_g_hat_x1;
scalllingl_x=JJ_star_g_lambda_x1(end)/costt0_g_lambda_x1(end);
costt0_g_lambda_x1=scalllingl_x*costt0_g_lambda_x1;

figure;
subplot(2,2,1)
plot(0:maximum_trail-1,costt0_g_lambda_x1,0:maximum_trail-1,costt0_g_hat_x1,'LineWidth',2)
hold on
%plot(0:maximum_trail-1,JJ_star_v0_x1,'g',0:maximum_trail-1,JJ_star_vl_x1,'g','LineWidth',2)
legend('g(x_1(.))','g^\lambda(x_1(.))')
title({'The Convergence of The Average Gradient Cost Functions for ';' The Two Velocities by Using TSRFN-VGL(\lambda=0.98) )';' and TSRFN-NSVGL(\lambda=0.98 with 0.001 Random Noise of one run ';' ';' The First State Gradient Cost Functions of ';' The TSRFN-VGL(0.98) '})
xlabel('Iteratons')
%ylabel('The First State Gradient Cost Functions')
grid on



scallling0_x1_SN=JJ_star_g_hat_x1_SN(end)/costt0_g_hat_x1_SN(end);
costt0_g_hat_x1_SN=scallling0_x1_SN*costt0_g_hat_x1_SN;
scalllingl_x_SN=JJ_star_g_lambda_x1_SN(end)/costt0_g_lambda_x1_SN(end);
costt0_g_lambda_x1_SN=scalllingl_x_SN*costt0_g_lambda_x1_SN;

subplot(2,2,2)
plot(0:maximum_trail-1,costt0_g_hat_x1_SN,0:maximum_trail-1,costt0_g_lambda_x1_SN,'LineWidth',2)
hold on
%plot(0:maximum_trail-1,JJ_star_v0_x1,'g',0:maximum_trail-1,JJ_star_vl_x1,'g','LineWidth',2)
legend('g(x_1(.))','g^\lambda(x_1(.))')
title({'The First State Gradient Cost Function of ';' The TSRFN-SNVGL(0.98)'})
xlabel('Iteratons')
%ylabel('The First State Gradient Cost Function')
grid on




%x2
scallling0_x2=JJ_star_g_hat_x2(end)/costt0_g_hat_x2(end);
costt0_g_hat_x2=abs(scallling0_x2*costt0_g_hat_x2);
scalllingl_x2=JJ_star_g_lambda_x2(end)/costt0_g_lambda_x2(end);
costt0_g_lambda_x2=abs(scalllingl_x2*costt0_g_lambda_x2);

subplot(2,2,3)
plot(0:maximum_trail-1,costt0_g_lambda_x2-0.06,0:maximum_trail-1,costt0_g_hat_x2-0.06,'LineWidth',2)
hold on
%plot(0:maximum_trail-1,JJ_star_v0_x1,'g',0:maximum_trail-1,JJ_star_vl_x1,'g','LineWidth',2)
legend('g(x_2(.))','g^\lambda(x_2(.))')
title({'The Second State Gradient Cost Function of ';' The TSRFN-VGL(0.98)'})
xlabel('Iteratons')
%ylabel('The Second State Gradient Cost Function')
grid on


scallling0_x2_SN=JJ_star_g_hat_x2_SN(end)/costt0_g_hat_x2_SN(end);
costt0_g_hat_x2_SN=abs(scallling0_x2_SN*costt0_g_hat_x2_SN);
scalllingl_x2_SN=JJ_star_g_lambda_x2_SN(end)/costt0_g_lambda_x2_SN(end);
costt0_g_lambda_x2_SN=abs(scalllingl_x2_SN*costt0_g_lambda_x2_SN);

subplot(2,2,4)
plot(0:maximum_trail-1,costt0_g_hat_x2_SN,0:maximum_trail-1,costt0_g_lambda_x2_SN,'LineWidth',2)
hold on
%plot(0:maximum_trail-1,JJ_star_v0_x1,'g',0:maximum_trail-1,JJ_star_vl_x1,'g','LineWidth',2)
legend('g(x_2(.))','g^\lambda(x_2(.))')
title({'The Second State Gradient Cost Function of ';' The TSRFN-SNVGL(0.98)'})
xlabel('Iteratons')
%ylabel('The Second State Gradient Cost Function')
grid on





fdfd=fdfdg