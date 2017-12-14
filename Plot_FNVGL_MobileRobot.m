% DHDP(lambda)

clear;
clc;
close all

F_veiw=20;% itersation no. to view for after first iteration
S_veiw=300;% itersation no. to view for after F_veiw

%-----------------
%loading a max_run
filenammm0=['Parameters_05_X' num2str(1) '.mat'];
load(filenammm0);
clear s_W_O1_c s_Theta2_I1_c s_W_O2_c s_Theta2_I2_c s_W_O1_a s_Theta2_I1_a s_W_O2_a s_Theta2_I2_a Av_g_x1 Av_R_lambda_x1 Av_g_x2 Av_R_lambda_x2 Total_g R_lambda_x1 Av_R_lambda_x1 R_lambda_x2 Av_R_lambda_x2 Total_MSE maximum_trail N_steps  actions_1 actions_2 Av_actions_1 Av_actions_2 states_1 states_2 MSE_1 MSE_2 Cost_Iter_U Av_Cost_Iter_J JJ Cost_Iter_v_0_1 Av_Cost_Iter_v_0_1 Cost_Iter_v_0_2 Av_Cost_Iter_v_0_2 ea Ea ec Ec
%--------------------------------------------------------------------------
% Number of time steps per episode
for run_count=1:5%max_runs
    %---------------------
    % MVGL  (lambda=0)
    filenammm0=['Parameters_01_X' num2str(run_count) '.mat'];
    load(filenammm0);
    
    Total_MSE_runs(1,1:maximum_trail,run_count)=Total_MSE;% mse iterations \all iterations
    
    iter_no=2; % iteration number
    Actionss_Itr_1(1,1:N_steps,run_count)=actions_1(1:end,iter_no); % (actions) first iteration \all iterations
    Actionss_Itr_3(1,1:N_steps,run_count)=actions_1(1:end,F_veiw); % (actions) 3 iteration \all iterations
    Actionss_Itr_5(1,1:N_steps,run_count)=actions_1(1:end,S_veiw); % (actions) 5 iteration \all iterations
    Actionss_Itr_end(1,1:N_steps,run_count)=actions_1(1:end,maximum_trail); % (actions) final iteration \all iterations
    
    Av_Actionss_Itr_run(1,1:maximum_trail,run_count)= Av_actions_1(1:maximum_trail); % (average actions) all iterations
    
    Actionss2_Itr_1(1,1:N_steps,run_count)=actions_2(1:end,iter_no); % (actions) first iteration \all iterations
    Actionss2_Itr_3(1,1:N_steps,run_count)=actions_2(1:end,F_veiw); % (actions) 3 iteration \all iterations
    Actionss2_Itr_5(1,1:N_steps,run_count)=actions_2(1:end,S_veiw); % (actions) 5 iteration \all iterations
    Actionss2_Itr_end(1,1:N_steps,run_count)=actions_2(1:end,maximum_trail); % (actions) final iteration \all iterations
    
    Av_Actions2_Itr_end(1,1:maximum_trail,run_count)= Av_actions_2(1:maximum_trail); % (average actions) all iterations
    
    states_1_Itr_1(1,1:N_steps,run_count)=states_1(1:end-1,iter_no); %(stats1) first iteration \all runs
    states_1_Itr_3(1,1:N_steps,run_count)=states_1(1:end-1,F_veiw); %(stats1) 3 iteration \all runs
    states_1_Itr_5(1,1:N_steps,run_count)=states_1(1:end-1,S_veiw); %(stats1) 5 iteration \all runs
    states_1_Itr_end(1,1:N_steps,run_count)=states_1(1:end-1,maximum_trail); % (stats1) final iteration \all runs
    
    states_2_Itr_1(1,1:N_steps,run_count)=states_2(1:end-1,iter_no); %(stats2) first iteration \all runs
    states_2_Itr_3(1,1:N_steps,run_count)=states_2(1:end-1,F_veiw); %(stats2) 3 iteration \all runs
    states_2_Itr_5(1,1:N_steps,run_count)=states_2(1:end-1,S_veiw); %(stats2) 5 iteration \all runs
    states_2_Itr_end(1,1:N_steps,run_count)=states_2(1:end-1,maximum_trail); % (stats2) final iteration \all runs
    

        Cost_g_x1_runs_SN(1,1:maximum_trail,run_count)=Av_g_x1;% Average cost (x1 g0) iterations \all runs
    Cost_g_x1_runs_SN(2,1:maximum_trail,run_count)=Av_R_lambda_x1;% Average target of cost (x1 g0) iterations \all runs
    
    Cost_g_x2_runs_SN(1,1:maximum_trail,run_count)=Av_g_x2;% Average cost (x2 g0) iterations \all runs
    Cost_g_x2_runs_SN(2,1:maximum_trail,run_count)=Av_R_lambda_x2;% Average target of cost (x2 g0) iterations \all runs
    
    

    
    %    Ea_runs(1,1:N_steps,run_count)=mean(Ea);% average for both actor state errors
    
    Ec_runs(1,1:N_steps,run_count)=mean(Ec);% average for both critic state errors
    
    clear s_W_O1_c s_Theta2_I1_c s_W_O2_c s_Theta2_I2_c s_W_O1_a s_Theta2_I1_a s_W_O2_a s_Theta2_I2_a Av_g_x1 Av_R_lambda_x1 Av_g_x2 Av_R_lambda_x2 Total_g R_lambda_x1 Av_R_lambda_x1 R_lambda_x2 Av_R_lambda_x2 Total_MSE maximum_trail N_steps Vd x0d0_F_c x0d0_S_c x0d0_F_a x0d0_S_a actions_1 actions_2 Av_actions_1 Av_actions_2 states_1 states_2 MSE_1 MSE_2 Cost_Iter_U Av_Cost_Iter_J JJ Cost_Iter_v_0_1 Av_Cost_Iter_v_0_1 Cost_Iter_v_0_2 Av_Cost_Iter_v_0_2 ea Ea ec Ec
    %---------------------
    % MVGL  (lambda=0.5)
    filenammm1=['Parameters_05_X' num2str(run_count) '.mat'];
    load(filenammm1);
    Total_MSE_runs(2,1:maximum_trail,run_count)=Total_MSE;% mse iterations \all iterations
    iter_no=2; % iteration number
    Actionss_Itr_1(2,1:N_steps,run_count)=actions_1(1:end,iter_no); % (actions) first iteration \all iterations
    Actionss_Itr_3(2,1:N_steps,run_count)=actions_1(1:end,F_veiw); % (actions) 3 iteration \all iterations
    Actionss_Itr_5(2,1:N_steps,run_count)=actions_1(1:end,S_veiw); % (actions) 5 iteration \all iterations
    Actionss_Itr_end(2,1:N_steps,run_count)=actions_1(1:end,maximum_trail); % (actions) final iteration \all iterations
    
    Av_Actionss_Itr_run(2,1:maximum_trail,run_count)= Av_actions_1(1:maximum_trail); % (average actions) all iterations
    
    Actionss2_Itr_1(2,1:N_steps,run_count)=actions_2(1:end,iter_no); % (actions) first iteration \all iterations
    Actionss2_Itr_3(2,1:N_steps,run_count)=actions_2(1:end,F_veiw); % (actions) 3 iteration \all iterations
    Actionss2_Itr_5(2,1:N_steps,run_count)=actions_2(1:end,S_veiw); % (actions) 5 iteration \all iterations
    Actionss2_Itr_end(2,1:N_steps,run_count)=actions_2(1:end,maximum_trail); % (actions) final iteration \all iterations
    
    Av_Actions2_Itr_end(2,1:maximum_trail,run_count)= Av_actions_2(1:maximum_trail); % (average actions) all iterations
    
    states_1_Itr_1(2,1:N_steps,run_count)=states_1(1:end-1,iter_no); %(stats1) first iteration \all runs
    states_1_Itr_3(2,1:N_steps,run_count)=states_1(1:end-1,F_veiw); %(stats1) 3 iteration \all runs
    states_1_Itr_5(2,1:N_steps,run_count)=states_1(1:end-1,S_veiw); %(stats1) 5 iteration \all runs
    states_1_Itr_end(2,1:N_steps,run_count)=states_1(1:end-1,maximum_trail); % (stats1) final iteration \all runs
    
    states_2_Itr_1(2,1:N_steps,run_count)=states_2(1:end-1,iter_no); %(stats2) first iteration \all runs
    states_2_Itr_3(2,1:N_steps,run_count)=states_2(1:end-1,F_veiw); %(stats2) 3 iteration \all runs
    states_2_Itr_5(2,1:N_steps,run_count)=states_2(1:end-1,S_veiw); %(stats2) 5 iteration \all runs
    states_2_Itr_end(2,1:N_steps,run_count)=states_2(1:end-1,maximum_trail); % (stats2) final iteration \all runs
    
    
%    Ea_runs(2,1:N_steps,run_count)=mean(Ea);% average for both actor state errors
    
    Ec_runs(2,1:N_steps,run_count)=mean(Ec);% average for both critic state errors
    
    
    Cost_g_x1_runs(1,1:maximum_trail,run_count)=Av_g_x1;% Average cost (x1 g0) iterations \all runs
    Cost_g_x1_runs(2,1:maximum_trail,run_count)=Av_R_lambda_x1;% Average target of cost (x1 g0) iterations \all runs
    
    Cost_g_x2_runs(1,1:maximum_trail,run_count)=Av_g_x2;% Average cost (x2 g0) iterations \all runs
    Cost_g_x2_runs(2,1:maximum_trail,run_count)=Av_R_lambda_x2;% Average target of cost (x2 g0) iterations \all runs
    
    %Cost_g_x1_runs(2,1:maximum_trail,run_count)=Av_Cost_Iter_v_0_1(1:end-1) ;% Average cost (x1 g0) iterations \all runs
    %Cost_g_x2_runs(2,1:maximum_trail,run_count)=Av_Cost_Iter_v_0_2(1:end-1) ;% Average cost (x2 g0) iterations \all runs
    
    clear Av_g_x1 Av_R_lambda_x1 Av_g_x2 Av_R_lambda_x2 Total_g R_lambda_x1 Av_R_lambda_x1 R_lambda_x2 Av_R_lambda_x2 Total_MSE  actions_1 actions_2 Av_actions_1 Av_actions_2 states_1 states_2 MSE_1 MSE_2 Cost_Iter_U Av_Cost_Iter_J JJ Cost_Iter_v_0_1 Av_Cost_Iter_v_0_1 Cost_Iter_v_0_2 Av_Cost_Iter_v_0_2 ea Ea ec Ec
    %---------------------
    run_count
end

%-----------------------------------------------------------------
%plot
%plot of MSE/episode
y_mean1 = mean(Total_MSE_runs,3); y_std1 = std(Total_MSE_runs,[],3);
figure;
mseb(0:maximum_trail-1,y_mean1,y_std1);
%ylim([-50 150])
legend('Total MSE for DHP (\lambda=0)','Total MSE with \lambda=0.98')
title({'Average Total MSEs for Two States of RFNN-VGL(\lambda)';' to Control on Mobile Robot Dynamic Model (Ten Runs)'})
xlabel('Iteratons')
ylabel('MSE')
grid on


% improvement_MVGL=((mean(y_mean1(1,:))-mean(mean(y_mean1(2,:)))/mean(mean(y_mean1(1,:)))*100;

%------------------
%plot of action/ steps ( frist and final iteration)
%------------------ frist iteration
y_mean2 = mean(Actionss_Itr_1,3); y_std2 = std(Actionss_Itr_1,[],3);
figure;
subplot(2,2,1)
mseb(0:N_steps-1,y_mean2,y_std2);
grid on
legend('\mu(x_k) for DHP (\lambda=0)','\mu(x_k) with \lambda=0.98')
title({'The Mean of Iterations for The Control Action (Right Torque) for ';' RFNN-VGL(\lambda) (Ten Runs) ';' ';' The First Iteration'})
xlabel('Time Steps (k)')
ylabel( 'The Frist Control Input(\mu(x_k))')
%------------------ thrid iteration
y_mean2_3 = mean(Actionss_Itr_3,3); y_std2_3 = std(Actionss_Itr_3,[],3);
%figure;
subplot(2,2,2)
mseb(0:N_steps-1,y_mean2_3,y_std2_3);
grid on
legend('\mu(x_k) for DHP (\lambda=0)','\mu(x_k) with \lambda=0.98')
title({'The Third Iteration'})
xlabel('Time Steps (k)')
ylabel( 'The Frist Control Input(\mu(x_k))')
%------------------ fifth iteration
y_mean2_5 = mean(Actionss_Itr_5,3); y_std2_5 = std(Actionss_Itr_5,[],3);
%figure;
subplot(2,2,3)
mseb(0:N_steps-1,y_mean2_5,y_std2_5);
grid on
legend('\mu(x_k) for DHP (\lambda=0)','\mu(x_k) with \lambda=0.98')
title({'The Fifth Iteration'})
xlabel('Time Steps (k)')
ylabel( 'The Frist Control Input(\mu(x_k))')
%------------------ last iteration
y_mean3 = mean(Actionss_Itr_end,3); y_std3 = std(Actionss_Itr_end,[],3);
%figure;
subplot(2,2,4)
mseb(0:N_steps-1,y_mean3,y_std3);
grid on
legend('\mu(x_k) for DHP (\lambda=0)','\mu(x_k) with \lambda=0.98')
title({' The Last Iteration'})
xlabel('Time Steps (k)')
ylabel( 'The Frist Control Input(\mu(x_k))')
%------------------
%avarage of right torque for all iterations (first action)
y_mean3ava1 = mean(Av_Actionss_Itr_run,3); y_std3ava1 = std(Av_Actionss_Itr_run,[],3);
%figure;
figure();
mseb(0:maximum_trail-1,y_mean3ava1,y_std3ava1);
grid on
legend('\mu(.) for DHP (\lambda=0)','\mu(.) with \lambda=0.98')
title({'The Mean of Iterations for The Average of Right Torque for ';' RFNN-VGL(\lambda) (Ten Runs)'})
xlabel('Iterations')
ylabel( '\mu(.)')



%------------------
%plot of action/ steps ( frist and final iteration)
%------------------ frist iteration
y_mean22 = mean(Actionss2_Itr_1,3); y_std22 = std(Actionss2_Itr_1,[],3);
figure;
subplot(2,2,1)
mseb(0:N_steps-1,y_mean22,y_std2);
grid on
legend('\mu(x_k) for DHP (\lambda=0)','\mu(x_k) with \lambda=0.98')
title({'The Mean of Iiterations for The Control Action (Left Torque) for ';' online MVGL(\lambda) Approach (Ten Runs) ';' ';' The First Iteration'})
xlabel('Time Steps (k)')
ylabel( 'The Second Control Input(\mu(x_k))')
%------------------ thrid iteration
y_mean22_3 = mean(Actionss2_Itr_3,3); y_std22_3 = std(Actionss2_Itr_3,[],3);
%figure;
subplot(2,2,2)
mseb(0:N_steps-1,y_mean22_3,y_std22_3);
grid on
legend('\mu(x_k) for DHP (\lambda=0)','\mu(x_k) with \lambda=0.98')
title({'The Third Iteration'})
xlabel('Time Steps (k)')
ylabel( 'The Second Control Input(\mu(x_k))')
%------------------ fifth iteration
y_mean22_5 = mean(Actionss2_Itr_5,3); y_std22_5 = std(Actionss2_Itr_5,[],3);
%figure;
subplot(2,2,3)
mseb(0:N_steps-1,y_mean22_5,y_std22_5);
grid on
legend('\mu(x_k) for DHP (\lambda=0)','\mu(x_k) with \lambda=0.98')
title({'The Fifth Iteration'})
xlabel('Time Steps (k)')
ylabel( 'The Second Control Input(\mu(x_k))')
%------------------ last iteration
y_mean32 = mean(Actionss2_Itr_end,3); y_std32 = std(Actionss2_Itr_end,[],3);
%figure;
subplot(2,2,4)
mseb(0:N_steps-1,y_mean32,y_std32);
grid on
legend('\mu(x_k) for DHP (\lambda=0)','\mu(x_k)with \lambda=0.5','\mu(x_k) with \lambda=0.99')
title({' The Last Iteration'})
xlabel('Time Steps (k)')
ylabel( 'The Second Control Input(\mu(x_k))')
%------------------


%avarage of right torque for all iterations (first action)
y_mean3ava2 = mean(Av_Actions2_Itr_end,3); y_std3ava2 = std(Av_Actions2_Itr_end,[],3);
%figure;
figure();
mseb(0:maximum_trail-1,y_mean3ava2,y_std3ava2);
grid on
legend('\mu(.) for DHP (\lambda=0)','\mu(.) with \lambda=0.98')
title({'The Mean of Iterations for The Average of Left Torque for ';' RFNN-VGL(\lambda) (Ten Runs)'})
xlabel('Iterations')
ylabel( '\mu(.)')


%plot of states/ steps ( frist and final iteration)
%------------------ frist iteration x1
y_mean4 = mean(states_1_Itr_1,3); y_std4 = std(states_1_Itr_1,[],3);
figure;
subplot(2,2,1)
grid on
mseb(0:N_steps-1,y_mean4,y_std4);
legend('x1_k for DHP (\lambda=0)','x1_k with \lambda=0.98')
title({'The Mean of Iiterations for The Frist State Trajectory (Linear Velocity) for ';' RFNN-VGL(\lambda) Approach (Ten Runs) ';' ';' The First Iteration'})
xlabel('Time Steps (k)')
ylabel('The Frist State Trajectory (x1_k)')
%------------------ third iteration x1
y_mean4_3 = mean(states_1_Itr_3,3); y_std4_3 = std(states_1_Itr_3,[],3);
%figure;
subplot(2,2,2)
grid on
mseb(0:N_steps-1,y_mean4_3,y_std4_3);
legend('x1_k for DHP (\lambda=0)','x1_k with \lambda=0.98')
title({' The Third Iteration'})
xlabel('Time Steps (k)')
ylabel('The Frist State Trajectory (x1_k)')
%------------------ fifth iteration x1
y_mean4_5 = mean(states_1_Itr_5,3); y_std4_5 = std(states_1_Itr_5,[],3);
%figure;
subplot(2,2,3)
grid on
mseb(0:N_steps-1,y_mean4_5,y_std4_5);
legend('x1_k for DHP (\lambda=0)','x1_k with \lambda=0.98')
title({' The Fifth Iteration'})
xlabel('Time Steps (k)')
ylabel('The Frist State Trajectory (x1_k)')
%------------------ last iteration x1
y_mean5 = mean(states_1_Itr_end,3); y_std5 = std(states_1_Itr_end,[],3);
%figure;
subplot(2,2,4)
grid on
mseb(0:N_steps-1,y_mean5,y_std5);
legend('x1_k for DHP (\lambda=0)','x1_k with \lambda=0.98')
title({' The Last Iteration'})
xlabel('Time Steps (k)')
ylabel('The Frist State Trajectory (x1_k)')




%------------------ frist iteration x2
y_mean6 = mean(states_2_Itr_1,3); y_std6 = std(states_2_Itr_1,[],3);
figure;
subplot(2,2,1)
grid on
mseb(0:N_steps-1,y_mean6,y_std6);
legend('x2_k for DHP (\lambda=0)','x2_k with \lambda=0.98')
title({'The Mean of Iiterations for The Second State Trajectory (Angular Velocity) for ';' RFNN-VGL(\lambda) Approach (Ten Runs) ';' ';' The First Iteration'})
xlabel('Time Steps (k)')
ylabel('The Second State Trajectory (x2_k)')
%------------------ third iteration x2
y_mean6_3 = mean(states_2_Itr_3,3); y_std6_3 = std(states_2_Itr_3,[],3);
%figure;
subplot(2,2,2)
grid on
mseb(0:N_steps-1,y_mean6_3,y_std6_3);
legend('x2_k for DHP (\lambda=0)','x2_k with \lambda=0.98')
title({'The Third Iteration'})
xlabel('Time Steps (k)')
ylabel('The Second State Trajectory (x2_k)')
%------------------ Fifth iteration x2
y_mean6_5 = mean(states_2_Itr_5,3); y_std6_5 = std(states_2_Itr_5,[],3);
%figure;
subplot(2,2,3)
grid on
mseb(0:N_steps-1,y_mean6_5,y_std6_5);
legend('x2_k for DHP (\lambda=0)','x2_k with \lambda=0.98')
title({'The Fifth Iteration'})
xlabel('Time Steps (k)')
ylabel('The Second State Trajectory (x2_k)')
%------------------ last iteration x2
y_mean7 = mean(states_2_Itr_end,3); y_std7 = std(states_2_Itr_end,[],3);
%figure;
subplot(2,2,4)
grid on
mseb(0:N_steps-1,y_mean7,y_std7);
legend('x2_k for DHP (\lambda=0)','x2_k with \lambda=0.98')
title({'The Last Iteration'})
xlabel('Time Steps (k)')
ylabel('The Second State Trajectory (x2_k)')

%---------------------------------------------------------------------
% %average of ciritc and actor errors
% y_mean7_Ea = mean(Ea_runs,3); y_std7_Ea = std(Ea_runs,[],3);
% figure;
% subplot(2,1,1)
% grid on
% mseb(0:N_steps-1,y_mean7_Ea,y_std7_Ea);
% legend('Ea_k for DHP (\lambda=0)','Ea_k with \lambda=0.98')
% title({'The Mean of Squared Actor Error for The last Iteration of ';' RFNN-VGL(\lambda) Approach (Ten Runs)'})
% xlabel('Time Steps (k)')
% ylabel('Ea_k')

y_mean7_Ec = mean(Ec_runs,3); y_std7_Ec = std(Ec_runs,[],3);
%figure;
subplot(2,1,2)
grid on
mseb(0:N_steps-1,y_mean7_Ec,y_std7_Ec);
legend('Ec_k for DHP (\lambda=0)','Ec_k with \lambda=0.98')
title({'The Mean of Squared Critic Error for The last Iteration of ';' RFNN-VGL(\lambda) Approach (Ten Runs)'})
xlabel('Time Steps (k)')
ylabel('Ec_k')
%----------------------------------------------------------------------
%Cost
% gx1
figure;
y_mean12 = mean(Cost_g_x1_runs,3); y_std12 = std(Cost_g_x1_runs,[],3);
mseb(0:maximum_trail-1,y_mean12,y_std12);
hold on
JJ_star_v0(1:maximum_trail)=y_mean12(1,end);
JJ_star_vl(1:maximum_trail)=y_mean12(2,end);
%plot(0:maximum_trail-1,JJ_star_v0_x1,'g',0:maximum_trail-1,JJ_star_vl_x1,'g','LineWidth',2)
% ylim([1.82 1.98])
% xlim([0 201])
legend('g(x_1)','R^\lambda(x_1)')
title({'The Convergence Process of The Average Gradient Cost Function and Its Target for ';' The First State for RFNN-VGL\lambda) Approach with \lambda =0.98 (Ten Run)'})
xlabel('Iteratons')
ylabel('The First State Gradient Cost Function with Its Target')
grid on

%gx2
figure;
y_mean135 = mean(Cost_g_x2_runs,3); y_std135 = std(Cost_g_x2_runs,[],3);
mseb(0:maximum_trail-1,y_mean135,y_std135);
hold on
JJ_star_g0_x1(1:maximum_trail)=y_mean135(1,end);
JJ_star_gl_x1(1:maximum_trail)=y_mean135(2,end);
%plot(0:maximum_trail-1,JJ_star_v0_x2,'g',0:maximum_trail-1,JJ_star_vl_x2,'g','LineWidth',2)
% ylim([1.82 1.98])
% xlim([0 201])
legend('g(x_2)','R^\lambda(x_2)')
title({'The Convergence Process of The Average Gradient Cost Function and Its Target for ';' The Seconf State for RFNN-VGL\lambda) Approach with \lambda =0.98 (Ten Run)'})
xlabel('Iteratons')
ylabel('The First State Gradient Cost Function with Its Target')
grid on


y_mean12_SN = mean(Cost_g_x1_runs_SN,3); y_std12 = std(Cost_g_x1_runs_SN,[],3);
y_mean135_SN = mean(Cost_g_x2_runs_SN,3); y_std135 = std(Cost_g_x2_runs_SN,[],3);

    
    














states_v_0=states_1_Itr_end(1,:);
states_w_0=states_2_Itr_end(1,:);

states_v_l=states_1_Itr_end(2,:);
states_w_l=states_2_Itr_end(2,:);

Ts=0.01;
% intgration for Desired velocities
Va_0=[states_v_0',states_w_0']; % Actual states
Va_l=[states_v_l',states_w_l']; % Actual states
t1=1;
t1dot=zeros(1,N_steps);% Desired time
t2dot=zeros(1,N_steps);% Actual time
MR_Dstates=[0 0 0 0 0 Vd(1,:)];% Desired initial states
MR_Astates_0=[0 0 0 0 0 Vd(1,:)];% Actual initial states
MR_Astates_l=[0 0 0 0 0 Vd(1,:)];% Actual initial states

while t1< N_steps
    [t1dot(t1+1),MR_Dstates]=Intgration_fun(t1dot(t1),Ts,MR_Dstates);% Desired integration states
    [t2dot(t1+1),MR_Astates_0]=Intgration_fun(t2dot(t1),Ts,MR_Astates_0); % Actual integration states
    [t2dot(t1+1),MR_Astates_l]=Intgration_fun(t2dot(t1),Ts,MR_Astates_l); % Actual integration states
    
    MR_Dstates(6:7)=Vd(t1+1,:);% Desired velocity states
    MR_Astates_0(6:7)=Va_0(t1+1,:); % Actual velocity states
    MR_Astates_l(6:7)=Va_l(t1+1,:); % Actual velocity states
    S_pD(t1+1,:)=MR_Dstates(1:5);% Desired position states
    S_p_wd_L1_FNN_0(t1+1,:)=MR_Astates_0(1:5);% Actual position states
    S_p_wd_L1_FNN_l(t1+1,:)=MR_Astates_l(1:5);% Actual position states
    
    t1=t1+1;
end

% comparison between L0 and L1
figure
%plot(S_pD(:,1),S_pD(:,2),'g-.',S_p_wod_L0_NN(:,1),S_p_wod_L0_NN(:,2),'b--',S_p_wd_L1_NN(:,1),S_p_wd_L1_NN(:,2),'m-o',S_p_wod_L0_FNN(:,1),S_p_wod_L0_FNN(:,2),'k:',S_p_wd_L1_FNN(:,1),S_p_wd_L1_FNN(:,2),'r-','LineWidth',2)
plot(S_pD(:,1),S_pD(:,2),'g-.',S_p_wd_L1_FNN_l(:,1),S_p_wd_L1_FNN_l(:,2),'b-',S_p_wd_L1_FNN_0(:,1),S_p_wd_L1_FNN_0(:,2),'r-','LineWidth',2)
grid on
hold on
title({'X-Y Circle Trajectory for RFNN-VGL(\lambda)';' with Affecting of Both Distrubaces and Firections for The Last Learning Iteration'},'FontSize',9)
xlabel('x-coordinate [m]','FontSize',9)
ylabel('y-coordinate [m]','FontSize',9)
% l=legend('kinematic','dynamic with friction','dynamic with friction');
l=legend('Disired Trajectory','Actual Trajectory with \lambda=0.98','Actual Trajectory with \lambda=0');
set(l,'FontSize',9)
%axis ([-0.1 2 -.01 0.081 ])

















ffff=dfdf

%----------------------------------------------------------------------
% Convert the NN cost value to increamental cost value

costt0_x(1)=0;
costtl_x(1)=0;
costt0_x1(1)=0;
costt0_x2(1)=0;
costt1_x1(1)=0;
costt1_x2(1)=0;
for i=2:maximum_trail
    % lambda = 0
    e_v0_x1(i)=abs(abs(JJ_star_v0(i))-abs(y_mean12(1,i))); %v0x
    e_c0_x1(i)=abs(abs(JJ_star_g0_x1(i))-abs(y_mean135(1,i))); %g0x1
    e_c0_x2(i)=abs(abs(JJ_star_g0_x2(i))-abs(y_mean13(1,i)));%g0x2
    
    costt0_x(i)=costt0_x(i-1)+e_v0_x1(i);
    costt0_x1(i)=costt0_x1(i-1)+e_c0_x1(i);
    costt0_x2(i)=costt0_x2(i-1)+e_c0_x2(i);
    
    % lambda = 0.95
    e_vl_x1(i)=abs(abs(JJ_star_vl(i))-abs(y_mean12(2,i))); %vlx
    e_c1_x1(i)=abs(abs(JJ_star_gl_x1(i))-abs(y_mean135(2,i)));%glx1
    e_c1_x2(i)=abs(abs(JJ_star_gl_x2(i))-abs(y_mean13(2,i)));%glx2
    
    costtl_x(i)=costtl_x(i-1)+e_vl_x1(i);
    costt1_x1(i)=costt1_x1(i-1)+e_c1_x1(i);
    costt1_x2(i)=costt1_x2(i-1)+e_c1_x2(i);
    
    
    
end
%vx
scallling0_x=JJ_star_v0(end)/costt0_x(end);
costt0_x=scallling0_x*costt0_x;
scalllingl_x=JJ_star_vl(end)/costtl_x(end);
costtl_x=scalllingl_x*costtl_x;

figure;
plot(0:maximum_trail-1,costt0_x,0:maximum_trail-1,costtl_x,'LineWidth',2)
hold on
JJ_star_v0_x1(1:maximum_trail)=y_mean12(1,end);
JJ_star_vl_x1(1:maximum_trail)=y_mean12(2,end);
%plot(0:maximum_trail-1,JJ_star_v0_x1,'g',0:maximum_trail-1,JJ_star_vl_x1,'g','LineWidth',2)
legend('v^0(x_k)','v^\lambda(x_k)')
title({'The Convergence Process of The Average Gradient Cost Function for ';' The First State for on-line MGDHP(\lambda) Approach with \lambda =0.99'})
xlabel('Iteratons')
ylabel('The First State Gradient Cost Function')
grid on

%gx1
scallling0_x1=JJ_star_g0_x1(end)/costt0_x1(end);
costt0_x1=scallling0_x1*costt0_x1;
scalllingl_x1=JJ_star_gl_x1(end)/costt1_x1(end);
costt1_x1=scalllingl_x1*costt1_x1;

figure;
plot(0:maximum_trail-1,costt0_x1,0:maximum_trail-1,costt1_x1,'LineWidth',2)
hold on
JJ_star_v0_x1(1:maximum_trail)=y_mean12(1,end);
JJ_star_vl_x1(1:maximum_trail)=y_mean12(2,end);
%plot(0:maximum_trail-1,JJ_star_v0_x1,'g',0:maximum_trail-1,JJ_star_vl_x1,'g','LineWidth',2)
legend('g^0(x1_k)','g^\lambda(x1_k)')
title({'The Convergence Process of The Average Gradient Cost Function for ';' The First State for on-line MGDHP(\lambda) Approach with \lambda =0.99'})
xlabel('Iteratons')
ylabel('The First State Gradient Cost Function')
grid on



%gx2
scallling0_x2=JJ_star_g0_x2(end)/costt0_x2(end);
costt0_x2=scallling0_x2*costt0_x2;
scalllingl_x2=JJ_star_gl_x2(end)/costt1_x2(end);
costt1_x2=scalllingl_x2*costt1_x2;

figure;
plot(0:maximum_trail-1,costt0_x2,0:maximum_trail-1,costt1_x2,'LineWidth',2)
hold on
JJ_star_v0_x2(1:maximum_trail)=y_mean13(1,end);
JJ_star_vl_x2(1:maximum_trail)=y_mean13(2,end);
%plot(0:maximum_trail-1,JJ_star_v0_x2,'g',0:maximum_trail-1,JJ_star_vl_x2,'g','LineWidth',2)
legend('g^0(x2_k)','g^\lambda(x2_k)')
title({'The Convergence Process of The Average Gradient Cost Function for ';' The Second State for on-line MGDHP(\lambda) Approach with \lambda =0.99'})
xlabel('Iteratons')
ylabel('The Second State Gradient Cost Function')
grid on
%----------------------------------------------------------------------



















sdsdsd=sffg



load('Parameters_01_X1.mat')
load('Parameters_09_X1.mat')

improvement_MGDHP=((mean(Total_MSE_0)-mean(Total_MSE_lambda))/mean(Total_MSE_0))*100;
%improvement_MGDHP=((mean(Total_MSE_lambda)-mean(Total_MSE_0))/mean(Total_MSE_lambda))*100;
%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
figure(3)
subplot(2,1,1)
plot(0:maximum_trail-1,JJ,'g',0:maximum_trail-1,Av_Cost_Iter_v_0_1_lambda(1:end-1),'r-.',0:maximum_trail-1,Av_Cost_Iter_v_lambda_1_lambda(1:end-1),'b--','LineWidth',2);
grid on
legend('J','v_1^0 ','v_1^\lambda')
title({'The convergence process of the cost function for ';' online MGDHP(\lambda) Approach (Nonlinear system) ';' ';' The Cost Function for Frist State '})
xlabel('Iteratons')
ylabel('The Function Value')
%------
subplot(2,1,2)
plot(0:maximum_trail-1,JJ,'g',0:maximum_trail-1,Av_Cost_Iter_v_0_2_lambda(1:end-1),'r-.',0:maximum_trail-1,Av_Cost_Iter_v_lambda_2_lambda(1:end-1),'b--','LineWidth',2);
grid on
legend('J','v_2^0 ','v_2^\lambda')
title(' The Cost Function for Second State ')
xlabel('Iteratons')
ylabel('The Function Value')
%---------------------------------------------------------------------------
figure(4)
subplot(3,1,1)
plot(tt,Ec_0_lambda(1,:),'g-',tt,Ec_0_lambda(2,:),'r--','LineWidth',2);
grid on
legend('E_c^0 for The First State','E_c^0 for The Second State')
title({'Analyses Critic and Actor Errors for on-line MGDHP(\lambda) of Nonlinear System';' ';' The Critic Error networks for GC(0)' })
xlabel('Time Steps')
ylabel('$$Error$$','Interpreter','Latex','FontSize',16)
hold on
%-----
subplot(3,1,2)
plot(tt,Ec_l_lambda(1,:),'g-',tt,Ec_l_lambda(2,:),'r--','LineWidth',2);
grid on
legend('E_c^\lambda for The First State','E_c^\lambda for The Second State')
title({'The Square of Critic Error networks for GC(\lambda)' })
xlabel('Time Steps')
ylabel('$$Error$$','Interpreter','Latex','FontSize',16)
hold on
%------------Criric errors ec_0 and ec_lambda with Lambda =0.95
subplot(3,1,3)
plot(tt,Ea_lambda(1,:),'g-','LineWidth',2);
grid on
legend('E_a')
title({'The Square of Actor Error for GA' })
xlabel('Time Steps')
ylabel('$$Error$$','Interpreter','Latex','FontSize',16)
hold on
%--------------------

