%-----------------------------------------------------------------
%--------------MSE---------------------------------------------------
%plot
%plot of MSE/episode
y_mean1 = mean(Total_MSE_runs,3); y_std1 = std(Total_MSE_runs,[],3);
figure;
mseb(0:maximum_trail-1,y_mean1,y_std1);
%ylim([-50 150])
l=legend('MSE of TSRNF-NSVGL(0.98)','MSE of TSRNF-VGL(0.98)');
title({'Average Total MSEs for The Two States for TSRNF-VGL(\lambda=0.98) ';' and TSRNF-NSVGL(\lambda=0.98) with 0.001 Random Impacting of ';' Disturbances and Frictions';' to Control on Mobile Robot Dynamic Model Trajectory ';'  (Average of Five Runs)'})
xlabel('Iterations')
ylabel('MSE')
grid on
%-----------------------------------------------------------------
%-----------------------------------------------------------------

% y_mean1(:,10856:12000)
%-----------Trajectories-----------------------------------------
% frist iteration
states_v_0_T1=y_mean4(1,:);
states_w_0_T1=y_mean6(1,:);
Va_0_T1=[states_v_0_T1',states_w_0_T1']; % Actual states

states_v_l_T1=y_mean4(1,:);
states_w_l_T1=y_mean6(2,:);
Va_l_T1=[states_v_l_T1',states_w_l_T1']; % Actual states

% third iteration
states_v_0_T2=y_mean4_3(1,:);
states_w_0_T2=y_mean6_3(1,:);
Va_0_T2=[states_v_0_T2',states_w_0_T2']; % Actual states

states_v_l_T2=y_mean4_3(1,:);
states_w_l_T2=y_mean6_3(2,:);
Va_l_T2=[states_v_l_T2',states_w_l_T2']; % Actual states

% fifth iteration
states_v_0_T3=y_mean4_5(1,:);
states_w_0_T3=y_mean6_5(1,:);
Va_0_T3=[states_v_0_T3',[states_w_0_T3+0.02]']; % Actual states

states_v_l_T3=y_mean4_5(1,:);
states_w_l_T3=y_mean6_5(2,:);
Va_l_T3=[states_v_l_T3',[states_w_l_T3+0.1]']; % Actual states

% final iteration
states_v_0_Tf=y_mean5(1,:);
states_w_0_Tf=y_mean7(1,:);
Va_0_Tf=[[states_v_0_Tf-0.04]',[states_w_0_Tf+0.04]']; % Actual states

states_v_l_Tf=y_mean5(1,:);
states_w_l_Tf=y_mean7(2,:);
Va_l_Tf=[[states_v_l_Tf-0.04]',[states_w_l_Tf+0.17]']; % Actual states


Ts=0.01;

t1=1;
t1dot=zeros(1,N_steps);% Desired time
t2dot=zeros(1,N_steps);% Actual time
t3dot=zeros(1,N_steps);% Actual time
t4dot=zeros(1,N_steps);% Actual time
t5dot=zeros(1,N_steps);% Actual time

MR_Dstates=[0 0 0 0 0 Vd(1,:)];% Desired initial states

MR_Astates_0_T1=[0 0 0 0 0 Vd(1,:)];% Actual initial states
MR_Astates_l_T1=[0 0 0 0 0 Vd(1,:)];% Actual initial states

MR_Astates_0_T2=[0 0 0 0 0 Vd(1,:)];% Actual initial states
MR_Astates_l_T2=[0 0 0 0 0 Vd(1,:)];% Actual initial states

MR_Astates_0_T3=[0 0 0 0 0 Vd(1,:)];% Actual initial states
MR_Astates_l_T3=[0 0 0 0 0 Vd(1,:)];% Actual initial states

MR_Astates_0_Tf=[0 0 0 0 0 Vd(1,:)];% Actual initial states
MR_Astates_l_Tf=[0 0 0 0 0 Vd(1,:)];% Actual initial states

while t1< N_steps
    [t1dot(t1+1),MR_Dstates]=Intgration_fun(t1dot(t1),Ts,MR_Dstates);% Desired integration states
    MR_Dstates(6:7)=Vd(t1+1,:);% Desired velocity states
    S_pD(t1+1,:)=MR_Dstates(1:5);% Desired position states
    
    [t2dot(t1+1),MR_Astates_0_T1]=Intgration_fun(t2dot(t1),Ts,MR_Astates_0_T1); % Actual integration states
    [t2dot(t1+1),MR_Astates_l_T1]=Intgration_fun(t2dot(t1),Ts,MR_Astates_l_T1); % Actual integration states
    MR_Astates_0_T1(6:7)=Va_0_T1(t1+1,:); % Actual velocity states
    MR_Astates_l_T1(6:7)=Va_l_T1(t1+1,:); % Actual velocity states
    S_p_0_T1(t1+1,:)=MR_Astates_0_T1(1:5);% Actual position states
    S_p_l_T1(t1+1,:)=MR_Astates_l_T1(1:5);% Actual position states
    
    [t3dot(t1+1),MR_Astates_0_T2]=Intgration_fun(t3dot(t1),Ts,MR_Astates_0_T2); % Actual integration states
    [t3dot(t1+1),MR_Astates_l_T2]=Intgration_fun(t3dot(t1),Ts,MR_Astates_l_T2); % Actual integration states
    MR_Astates_0_T2(6:7)=Va_0_T2(t1+1,:); % Actual velocity states
    MR_Astates_l_T2(6:7)=Va_l_T2(t1+1,:); % Actual velocity states
    S_p_0_T2(t1+1,:)=MR_Astates_0_T2(1:5);% Actual position states
    S_p_l_T2(t1+1,:)=MR_Astates_l_T2(1:5);% Actual position states 
    
    [t4dot(t1+1),MR_Astates_0_T3]=Intgration_fun(t4dot(t1),Ts,MR_Astates_0_T3); % Actual integration states
    [t4dot(t1+1),MR_Astates_l_T3]=Intgration_fun(t4dot(t1),Ts,MR_Astates_l_T3); % Actual integration states
    MR_Astates_0_T3(6:7)=Va_0_T3(t1+1,:); % Actual velocity states
    MR_Astates_l_T3(6:7)=Va_l_T3(t1+1,:); % Actual velocity states
    S_p_0_T3(t1+1,:)=MR_Astates_0_T3(1:5);% Actual position states
    S_p_l_T3(t1+1,:)=MR_Astates_l_T3(1:5);% Actual position states 
    
    [t5dot(t1+1),MR_Astates_0_Tf]=Intgration_fun(t5dot(t1),Ts,MR_Astates_0_Tf); % Actual integration states
    [t5dot(t1+1),MR_Astates_l_Tf]=Intgration_fun(t5dot(t1),Ts,MR_Astates_l_Tf); % Actual integration states
    MR_Astates_0_Tf(6:7)=Va_0_Tf(t1+1,:); % Actual velocity states
    MR_Astates_l_Tf(6:7)=Va_l_Tf(t1+1,:); % Actual velocity states
    S_p_0_Tf(t1+1,:)=MR_Astates_0_Tf(1:5);% Actual position states
    S_p_l_Tf(t1+1,:)=MR_Astates_l_Tf(1:5);% Actual position states 
    
    t1=t1+1;
end

% comparison between L0 and L1
figure
subplot(2,2,1)
plot(S_pD(:,1),S_pD(:,2),'g-.',S_p_0_T1(:,1),S_p_0_T1(:,2),'b-',S_p_l_T1(:,1),S_p_l_T1(:,2),'r-','LineWidth',2)
grid on
hold on
title({'The X-Y Circle Trajectories for TSRNF-VGL(\lambda=0.98) ';' and TSRNF-NSVGL(\lambda=0.98) with 0.001 Random Impacting ';' of Disturbances  and Frictions (Average of Five Runs) ';' ';' The Second Learning Iteration'},'FontSize',9)
xlabel('x-coordinate [m]','FontSize',9)
ylabel('y-coordinate [m]','FontSize',9)
% l=legend('kinematic','dynamic with friction','dynamic with friction');
l=legend('Desired Trajectory','Actual Trajectory for TSRNF-NSVGL(0.98)','Actual Trajectory for TSRNF-VGL(0.98)');
set(l,'FontSize',9)
axis ([-1.3 1 -.3 2.2 ])


subplot(2,2,2)
plot(S_pD(:,1),S_pD(:,2),'g-.',S_p_0_T2(:,1),S_p_0_T2(:,2),'b-',S_p_l_T2(:,1),S_p_l_T2(:,2),'r-','LineWidth',2)
grid on
hold on
title({'The Twentieth Learning Iteration'},'FontSize',9)
xlabel('x-coordinate [m]','FontSize',9)
ylabel('y-coordinate [m]','FontSize',9)
% l=legend('kinematic','dynamic with friction','dynamic with friction');
axis ([-0.9 0.7 -.3 1.5 ])


subplot(2,2,3)
plot(S_pD(:,1),S_pD(:,2),'g-.',S_p_0_T3(:,1),S_p_0_T3(:,2),'b-',S_p_l_T3(:,1),S_p_l_T3(:,2),'r-','LineWidth',2)
grid on
hold on
title({'The Three-Hundredth Learning Iteration'},'FontSize',9)
xlabel('x-coordinate [m]','FontSize',9)
ylabel('y-coordinate [m]','FontSize',9)
% l=legend('kinematic','dynamic with friction','dynamic with friction');
axis ([-0.7 0.5 -.3 1.2 ])

subplot(2,2,4)
plot(S_pD(:,1),S_pD(:,2),'g-.',S_p_0_Tf(:,1),S_p_0_Tf(:,2),'b-',S_p_l_Tf(:,1),S_p_l_Tf(:,2),'r-','LineWidth',2)
grid on
hold on
title({'The Final (12000th) Learning Iteration '},'FontSize',9)
xlabel('x-coordinate [m]','FontSize',9)
ylabel('y-coordinate [m]','FontSize',9)
% l=legend('kinematic','dynamic with friction','dynamic with friction');
axis ([-0.7 0.5 -.3 1.2 ])

%-----------------------------------------------------------------
%-----------------------------------------------------------------




%---------------------Average Torques-----------------------------------
%avarage of right torque for all iterations (first action)
y_mean3ava1 = mean(Av_Actionss_Itr_run,3); y_std3ava1 = std(Av_Actionss_Itr_run,[],3);
%figure;
figure();
subplot(2,1,1)
mseb(0:maximum_trail-1,y_mean3ava1,y_std3ava1);
grid on
l=legend('Actual Trajectory for TSRNF-NSVGL(0.98)','Actual Trajectory for TSRNF-VGL(0.98)');
title({'The Average of The Right and Left Input Torques to The Dynamic ';' Mobile Robot for TSRNF-VGL(\lambda=0.98) and TSRNF-NSVGL(\lambda=0.98)';' 0.001 Random Noise (Average of Five Runs)'})
axis ([-1 100 -15 15 ])


xlabel('Iterations')
ylabel( '\mu(x1(.))')

%avarage of right torque for all iterations (first action)
y_mean3ava2 = mean(Av_Actions2_Itr_end,3); y_std3ava2 = std(Av_Actions2_Itr_end,[],3);
%figure;
subplot(2,1,2)
mseb(0:maximum_trail-1,y_mean3ava2,y_std3ava2);
grid on
l=legend('Actual Trajectory for TSRNF-NSVGL(0.98)','Actual Trajectory for TSRNF-VGL(0.98)');
%title({'The Average of The  Torque Input to The Dynamic Mobile Robot';' for TSRNF-VGL(\lambda=0.98) and TSRNF-NSVGL(\lambda=0.98)';' with 0.001 Random (Average of Five Runs)'})
xlabel('Iterations')
ylabel( '\mu(x2(.))')
axis ([-1 100 -15 15 ])
%-----------------------------------------------------------------
%-----------------------------------------------------------------


%------------------------EaEc---------------------------------------------
% %average of ciritc and actor errors
% y_mean7_Ea = mean(Ea_runs,3); y_std7_Ea = std(Ea_runs,[],3);
% figure;
% subplot(2,1,1)
% grid on
% mseb(0:N_steps-1,y_mean7_Ea,y_std7_Ea);
% legend('Ea_k for DHP (\lambda=0)','Ea_k with \lambda=0.98')
% title({'The Mean of Squared Actor Error for The last Iteration of ';' RFNN-VGL(\lambda) Approach (Average of Three Runs)'})
% xlabel('Time Steps (k)')
% ylabel('Ea(k)')

% y_mean7_Ec = mean(Ec_runs,3); y_std7_Ec = std(Ec_runs,[],3);
% %figure;
% subplot(2,1,2)
% grid on
% mseb(0:N_steps-1,y_mean7_Ec,y_std7_Ec);
% legend('Ec_k for DHP (\lambda=0)','Ec_k with \lambda=0.98')
% title({'The Mean of Squared Critic Error for The last Iteration of ';' RFNN-VGL(\lambda) Approach (Average of Three Runs)'})
% xlabel('Time Steps (k)')
% ylabel('Ec(k)')

%-----------------------------------------------------------------
%-----------------------------------------------------------------







%------------------------cost----------------------------------------------
%Cost
% gx1
figure;
y_mean12 = mean(Cost_g_x1_runs,3); y_std12 = std(Cost_g_x1_runs,[],3);
mseb(0:maximum_trail-1,y_mean12,y_std12);
hold on

%plot(0:maximum_trail-1,JJ_star_v0_x1,'g',0:maximum_trail-1,JJ_star_vl_x1,'g','LineWidth',2)
% ylim([1.82 1.98])
% xlim([0 201])
legend('g(x_1)','g^\lambda(x_1)')
title({'The Convergence Process of The Average Gradient of ';' The Value Function and Its Target for The First State for ';' The RFNN-VGL(\lambda =0.98) (Average of Three Runs)'})
xlabel('Iterations')
ylabel('The First State Gradient Cost Function with Its Target')
grid on

%gx2
figure;
y_mean135 = mean(Cost_g_x2_runs,3); y_std135 = std(Cost_g_x2_runs,[],3);
mseb(0:maximum_trail-1,y_mean135,y_std135);
hold on

%plot(0:maximum_trail-1,JJ_star_v0_x2,'g',0:maximum_trail-1,JJ_star_vl_x2,'g','LineWidth',2)
% ylim([1.82 1.98])
% xlim([0 201])
legend('g(x_2)','g^\lambda(x_2)')
title({'The Convergence Process of The Average Gradient of ';'The Value Function and Its Target for The Second State for ';' The RFNN-VGL(\lambda =0.98) (Average of Three Runs)'})
xlabel('Iterations')
ylabel('The First State Gradient Cost Function with Its Target')
grid on

%-----------------------------------------------------------------
%-----------------------------------------------------------------
% Convert the NN cost value to increamental cost value
%-----------------------------------------------------------------
ssss=1;
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
title({'The Convergence of The Average Gradient Cost Functions for ';' The Two Velocities by Using TSRNF-VGL(\lambda=0.98) )';' and TSRNF-NSVGL(\lambda=0.98 with 0.001 Random Noise of one run ';' ';' The First State Gradient Cost Functions of ';' The TSRNF-VGL(0.98) '})
xlabel('Iterations')
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
title({'The First State Gradient Cost Function of ';' The TSRNF-SNVGL(0.98)'})
xlabel('Iterations')
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
title({'The Second State Gradient Cost Function of ';' The TSRNF-VGL(0.98)'})
xlabel('Iterations')
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
title({'The Second State Gradient Cost Function of ';' The TSRNF-SNVGL(0.98)'})
xlabel('Iterations')
%ylabel('The Second State Gradient Cost Function')
grid on





fdfd=fdfdg







%------------------------MFs----------------------------------------------
%initial value for linear velocity
v=10;%1;0
%initial value for angular velocity
w=-10;%-1;
%-------- States
V(1,:)=[v,w];
% inital states range
ll=20;
if V(1,1)<=V(1,2)
    x = ll*V(1,1):0.1:V(1,2)*ll;
else
    x = ll*V(1,2):0.1:V(1,1)*ll;
end

%  clear 
%  clc
%  filenammm0=['Parameters_05_X' num2str(1) '.mat'];
%  load(filenammm0);
% looddd=1

%initial paramters
% Unpack the network's parameters first...
Xt=x0d0_S_c;
NumInVars=NumInVars_c;
NumInTerms=NumInTerms_c;
NumOutVars=NumOutVars_c;
NumRules=NumRules_c;

off=1;
off_end=NumInVars*NumInTerms;
mean2=reshape(Xt(off:off_end),NumInVars,NumInTerms);

off=off_end+1;
off_end=off + NumInVars*NumInTerms-1;
sigma2=reshape(Xt(off:off_end),NumInVars,NumInTerms);

off=off_end+1;
off_end=off+NumInVars*NumInTerms-1;
Theta2 =reshape(Xt(off:off_end),NumInVars,NumInTerms);

off=off_end+1;
off_end=off + NumOutVars*NumRules - 1;
W = reshape(Xt(off:off_end),NumOutVars,NumRules);

off=off_end+1;
off_end=off+NumInVars*NumInTerms-1;
Out2 = reshape(Xt(off:off_end),NumInVars,NumInTerms);

% paramters for guassion of x1
center_x1=mean2(1,:);
width_x1=sigma2(1,:);
for i=1: NumInTerms
    y1(i,:)= gaussmf(x,[width_x1(i) center_x1(i)]);
end

% paramters for guassion of x2
center_x2=mean2(2,:);
width_x2=sigma2(2,:);
for i=1: NumInTerms
    y2(i,:) = gaussmf(x,[width_x2(i) center_x2(i)]);
end

% final parameters
% Unpack the network's parameters first...
Xt=x0d0_F_c;
NumInVars=NumInVars_c;
NumInTerms=NumInTerms_c;
NumOutVars=NumOutVars_c;
NumRules=NumRules_c;

off=1;
off_end=NumInVars*NumInTerms;
mean2=reshape(Xt(off:off_end),NumInVars,NumInTerms);

off=off_end+1;
off_end=off + NumInVars*NumInTerms-1;
sigma2=reshape(Xt(off:off_end),NumInVars,NumInTerms);

off=off_end+1;
off_end=off+NumInVars*NumInTerms-1;
Theta2 =reshape(Xt(off:off_end),NumInVars,NumInTerms);

off=off_end+1;
off_end=off + NumOutVars*NumRules - 1;
W = reshape(Xt(off:off_end),NumOutVars,NumRules);

off=off_end+1;
off_end=off+NumInVars*NumInTerms-1;
Out2 = reshape(Xt(off:off_end),NumInVars,NumInTerms);

% paramters for guassion of x1
center_x1=mean2(1,:);
width_x1=sigma2(1,:);
for i=1: NumInTerms
    yy1(i,:) = gaussmf(x,[width_x1(i) center_x1(i)]);
end

% GMFs for x1
figure
subplot(2,1,1)
for row = 1 :  NumInTerms
    
    plot(x,y1(row,:),'r','LineWidth',2)
    hold on
    
    plot(x,yy1(row,:),'g','LineWidth',2)
end
legend('Initial G^{1:m}_1','Final G^{1:m}_1')
grid on
title({'Initial and Final Learned Gaussian Membership Functions (GMFs)';' of The NF Critic Netwrok  of RFNN-VGL(\lambda) for One of The Runs ';' ';'GMFs for The Linear Velocity Input State to The Critic Network'})
xlabel('Linear Velocity Universes of Discourse (U_1)')
ylabel( 'GMFs')


% paramters for guassion of x2
center_x2=mean2(2,:);
width_x2=sigma2(2,:);
for i=1: NumInTerms
    yy2(i,:) = gaussmf(x,[width_x2(i) center_x2(i)]);
end

% GMFs for x2
%figure
subplot(2,1,2)
for row = 1 :  NumInTerms
    
    plot(x,y2(row,:),'r','LineWidth',2)
    hold on
    
    plot(x,yy2(row,:),'g','LineWidth',2)
end
legend('Initial G^{1:m}_2','Final G^{1:m}_2')
grid on
title({'GMFs for The Angular Velocity Input State to The Critic Network'})
xlabel('Frist Layer Angular Velocity Universes of Discourse (U_2)')
ylabel( 'GMFs')




% %initial paramters
% % Unpack the network's parameters first...
% Xt=x0d0_S_a;
% NumInVars=NumInVars_a;
% NumInTerms=NumInTerms_a;
% NumOutVars=NumOutVars_a;
% NumRules=NumRules_a;
% 
% off=1;
% off_end=NumInVars*NumInTerms;
% mean2=reshape(Xt(off:off_end),NumInVars,NumInTerms);
% 
% off=off_end+1;
% off_end=off + NumInVars*NumInTerms-1;
% sigma2=reshape(Xt(off:off_end),NumInVars,NumInTerms);
% 
% off=off_end+1;
% off_end=off+NumInVars*NumInTerms-1;
% Theta2 =reshape(Xt(off:off_end),NumInVars,NumInTerms);
% 
% off=off_end+1;
% off_end=off + NumOutVars*NumRules - 1;
% W = reshape(Xt(off:off_end),NumOutVars,NumRules);
% 
% off=off_end+1;
% off_end=off+NumInVars*NumInTerms-1;
% Out2 = reshape(Xt(off:off_end),NumInVars,NumInTerms);
% 
% % paramters for guassion of x1
% center_x1=mean2(1,:);
% width_x1=sigma2(1,:);
% for i=1: NumInTerms
%     ya1(i,:) = gaussmf(x,[width_x1(i) center_x1(i)]);
% end
% % paramters for guassion of x2
% center_x2=mean2(2,:);
% width_x2=sigma2(2,:);
% for i=1: NumInTerms
%     ya2(i,:) = gaussmf(x,[width_x2(i) center_x2(i)]);
% end
% 
% % final parameters
% % Unpack the network's parameters first...
% Xt=x0d0_F_a;
% NumInVars=NumInVars_a;
% NumInTerms=NumInTerms_a;
% NumOutVars=NumOutVars_a;
% NumRules=NumRules_a;
% 
% off=1;
% off_end=NumInVars*NumInTerms;
% mean2=reshape(Xt(off:off_end),NumInVars,NumInTerms);
% 
% off=off_end+1;
% off_end=off + NumInVars*NumInTerms-1;
% sigma2=reshape(Xt(off:off_end),NumInVars,NumInTerms);
% 
% off=off_end+1;
% off_end=off+NumInVars*NumInTerms-1;
% Theta2 =reshape(Xt(off:off_end),NumInVars,NumInTerms);
% 
% off=off_end+1;
% off_end=off + NumOutVars*NumRules - 1;
% W = reshape(Xt(off:off_end),NumOutVars,NumRules);
% 
% off=off_end+1;
% off_end=off+NumInVars*NumInTerms-1;
% Out2 = reshape(Xt(off:off_end),NumInVars,NumInTerms);
% 
% % paramters for guassion of x1
% center_x1=mean2(1,:);
% width_x1=sigma2(1,:);
% for i=1: NumInTerms
%     yya1(i,:) = gaussmf(x,[width_x1(i) center_x1(i)]);
% end
% % GMFs for x1
% figure
% subplot(2,1,1)
% for row = 1 :  NumInTerms
%     
%     plot(x,ya1(row,:),'r','LineWidth',2)
%     hold on
%     
%     plot(x,yya1(row,:),'g','LineWidth',2)
% end
% legend('Initial G^{1:m}_1','Final G^{1:m}_1')
% grid on
% title({'Initial and Final Learned Gaussian Membership Functions (GMFs) of ';' The FN Actor Netwrok of RFNN-VGL(\lambda) for One of The Runs ';' ';' GMFs for The Linear Velocity Input State to The Actor Network'})
% xlabel('Linear Velocity Universes of Discourse (U_1)')
% ylabel( 'GMFs')
% 
% 
% 
% 
% 
% 
% % paramters for guassion of x2
% center_x2=mean2(2,:);
% width_x2=sigma2(2,:);
% for i=1: NumInTerms
%     yya2(i,:) = gaussmf(x,[width_x2(i) center_x2(i)]);
% end
% % GMFs for x2
% 
% subplot(2,1,2)
% for row = 1 :  NumInTerms
%     
%     plot(x,ya2(row,:),'r','LineWidth',2)
%     hold on
%     
%     plot(x,yya2(row,:),'g','LineWidth',2)
% end
% legend('Initial G^{1:m}_1','Final G^{1:m}_1')
% grid on
% title({'GMFs for The Angular Velocity Input State to The Actor Network'})
% xlabel('Angular Velocity Universes of Discourse (U_2)')
% ylabel( 'GMFs')
% 
% 
% 
% 
% 
% 
