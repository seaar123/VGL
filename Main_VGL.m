close all
clear all
clc
%---------------------------By_CR------------------------------------------
% Initial Values
%--------------------------------------------------------------------------
maximum_passes=100;% max iterations
Terminated_state=21;%final step
D_F=0.9;%Discount Factor
Lambda=1;
Omega=[1 0 ;0 1];
FR=0; % =1 with friction =0 without friction
%-------initial states x(t)
%initial value for linear velocity
v=2;
%initial value for angular velocity
w=-0.9;
%-------- States
V(1,:)=[v,w];
%--------------------
%Robot parameters
mT=10;  %Total mass
mW=2;   %Mass of 1 wheel
r=0.05; %Wheel radius
b=0.4;  %half the robot width
d=0.1;  %CG offset from rear axle
Iyy=1;  %Wheel moment of Iertia
IT=5;   %Platform total moment of inertia
Fv=.5 *FR;  %Coefficient of Fiscous friction ( *0 means ignoring friction )
Fd=.8 *FR;  %Coefficient of Colomb friction
%----------
% initial values for MR
x0=[0 0 0 0 0 V(1,:)];  %robot initial conditions [x,y,theta, phiR, phiL, v,w]
x0=[x0 mT mW r b d Iyy IT Fv Fd]; %pass in constants xdot=0 for these
%-------- Desired States
Vd(:,:)=zeros(Terminated_state,2);
Vd(1,:)=V(1,:);
%--------------------------------------------------------------------------
% parameters for Neural Netwrok
m_in=2; % number of input in input vector (torques)
n_op=2; % number of output in disired output vector (states)
N=12;% number of neurons in the hidden layers
n=N+n_op+1;% no. of neurons ( incluse inputs )
if N<m_in
    disp('ERROR: N less than m')
    return
end
%--------------------------------------------------------------------------
%Critic Network (CN)
%-------------------------------------------------------------------------
% initialize the weights
knd_c=1;% mode for activity function (1(linear),2(sigmoid),3(tanh),4(gaussian))
a_r = -0.3; %minimum range for wieght values
b_r =  0.3; %maximam range for wieght values
W=zeros(n-1,n);
for i=m_in+1+1:n
    for j=1:i-1
        W(j,i)=(b_r-a_r).*rand(1,1)+ a_r;
    end
end
WL0=W;
LR_c=0.01;%learning_rate for critic network
%--------------------------------------------------------------------------
%Action Network (AN)
%--------------------------------------------------------------------------
% initialize the weights
knd_a=3;% mode for activity function (1(linear),2(sigmoid),3(tanh),4(gaussian))
a_r = -0.3; %minimum range for wieght values
b_r =  0.3; %maximam range for wieght values
Z=zeros(n-1,n);
for i=m_in+1+1:n
    for j=1:i-1
        Z(j,i)=(b_r-a_r).*rand(1,1)+ a_r;
    end
end
ZL0=Z;
LR_a=0.01;%learning_rate for actor network
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Iteration
for Iter=1:maximum_passes
    t=1;
    error_v(Iter)=0;
    error_w(Iter)=0;
    % { Unroll trajectory ---}
    while t<Terminated_state
        %-------- Actions
        [a,s, u(t,:)]=NET(V(t,:),Z,m_in,N,n_op,knd_a);
        [xtt,df_dx,df_du]=Dynamic_MRobot(t,x0,u(t,:));
        V(t+1,:)=xtt(6:7);
        x0(1:7)=xtt;
        t=t+1;
    end
    [a,s, u(t,:)]=NET(V(t,:),Z,m_in,N,n_op,knd_a);
    F=t;
    [U,dU_dx,dU_du]=Cost_MRobot(V(t,:),u(t,:));
    
    p=dU_dx;
    [sss,DimV2]=Matrix_Converter(W,W,1);
    Wv=zeros(DimV2,n_op);
    Zv=zeros(DimV2,n_op);
    Delta_W=zeros(DimV2,1);
    Delta_Z=zeros(DimV2,1);
    % { Backwards pass ---}
    for t=F-1:-1:1
        %-----mr
        [U,dU_dx,dU_du]=Cost_MRobot(V(t,:),u(t,:));
        [xtt,df_dx,df_du]=Dynamic_MRobot(t,x0,u(t,:));
        %-----nn Actor
        [a_a,s_a,u(t,:)]=NET(V(t,:),Z,m_in,N,n_op,knd_a);
        for i=1:n_op
            F_Yhat_a(t,i)=1;
            % BW
            [dA_dw,dA_dx]=F_NET(Z,a_a,s_a,F_Yhat_a(t,:),m_in,N,i,knd_a); % eqs 10-12
            %**********Adaptive the Matrix Dimensions
            [Zv(:,i),DimV1]=Matrix_Converter(Z,Z,1);
            [sss,DimV2]=Matrix_Converter(dA_dw,Z,1);
            [dA_dw_v(1:DimV2,i),DimV2]=Matrix_Converter(dA_dw,Z,1);
            Zv(1:DimV2,i)=Zv(1:DimV2,i)+LR_a*dA_dw_v(1:DimV2,i);
            %[Z,DimV]=Matrix_Converter(Zv(:,i),Z,2);
        end
        %-----nn critic
        [a_c1,s_c1, G_telda_t(:,t)]=NET(V(t,:),W,m_in,N,n_op,knd_c);
        [a_c2,s_c2, G_telda_t(:,t+1)]=NET(V(t+1,:),W,m_in,N,n_op,knd_c);
        for i=1:n_op
            F_Yhat_c(t,i)=1;%Yhat_pattents(t,i)-Y(t,i); % error
            % BW
            [dW_dw,dW_dx]=F_NET(W,a_c1,s_c1,F_Yhat_c(t,:),m_in,N,i,knd_c); % eqs 10-12
            %**********Adaptive the Matrix Dimensions
            [Wv(:,i),DimV1]=Matrix_Converter(W,W,1);
            [sss,DimV2]=Matrix_Converter(dW_dw,Z,1);
            [dW_dw_v(1:DimV2,i),DimV2]=Matrix_Converter(dW_dw,W,1);
            Wv(1:DimV2,i)=Wv(1:DimV2,i)+LR_c*dW_dw_v(1:DimV2,i);
            %[W,DimV]=Matrix_Converter(Wv(:,i),W,2);
        end
        %-----Target G value
        TG(:,t)=dU_dx + (D_F * df_dx * p) + dA_dx *(dU_du + D_F * df_du * p);
        %------------------------------------------------------------------
        Delta_W=Delta_W + dW_dw_v*Omega*(TG(:,t)- G_telda_t(:,t));
        %errorrr=(TG(:,t)- G_telda_t(:,t))
        Delta_Z=Delta_Z - dA_dw_v*(dU_du + D_F * df_du * G_telda_t(:,t+1));
        p=Lambda * TG(:,t) + (1-Lambda)*G_telda_t(:,t);
    end
    %  Updating Weights
    [mDelta_W,DimV]=Matrix_Converter(Delta_W,W,2);
    [mDelta_Z,DimV]=Matrix_Converter(Delta_Z,Z,2);
    W= W + LR_c * mDelta_W;
    Z= Z + LR_a * mDelta_Z;
    %---------------
    % Storage actions and states
    action_TorqR(:,Iter)=u(:,1);
    action_TorqL(:,Iter)=u(:,2);
    
    states_v(:,Iter)=V(:,1);
    states_w(:,Iter)=V(:,2);
    %---------------
    % Error Linear velocity state
    for j=1:Terminated_state
        error_v(Iter)= error_v(Iter)+0.5*(states_v(j,Iter)-Vd(j,1))^2;
    end
    % Error angular velocity state
    for j=1:Terminated_state
        error_w(Iter)= error_w(Iter)+0.5*(states_w(j,Iter)-Vd(j,2))^2;
    end
    
    error_vwL1(Iter)= (error_v(Iter)+ error_w(Iter))/2;
    
    
end
TS=Terminated_state-1;

% %Torques plot ( actions )
% for Iter=2:20
%     figure(1)
%     plot(0:TS,action_TorqR(:,1),'r',0:TS,action_TorqR(:,Iter),'y',0:TS,action_TorqR(:,maximum_passes),'g')
%     grid on
%     figure(2)
%     plot(0:TS,action_TorqL(:,1),'r',0:TS,action_TorqL(:,Iter),'y',0:TS,action_TorqL(:,maximum_passes),'g')
%     grid on
%     hold on
% end

% velocities plot ( states )
for Iter1=2:20
    figure(3)
    plot(0:TS,states_v(:,1),'r',0:TS,states_v(:,Iter1),'y',0:TS,states_v(:,maximum_passes),'g','LineWidth',2);
    grid on
    legend('First Iteration','Other Iterations', 'Final Iteration')
    title('The Trajectories for Linear Velocity with Friction {\lambda} =1')
    xlabel('Steps')
    ylabel('Lenear Velocity')
    hold on
    figure(4)
    plot(0:TS,states_w(:,1),'r',0:TS,states_w(:,Iter1),'y',0:TS,states_w(:,maximum_passes),'g','LineWidth',2);
    grid on
    legend('First Iteration','Other Iterations', 'Final Iteration')
    title('The Trajectories for Angualar Velocity with Friction {\lambda} =1')
    xlabel('Steps')
    ylabel('Aangular Velocity')
    grid on
    hold on
end

figure(5)
plot(0:maximum_passes-1,error_v,'r',0:maximum_passes-1,error_w,'k','LineWidth',2);
grid on
legend('Error for Linear Velocity State ','Error for Angular Velocity State')
title('The Errors for Trajectories with Friction at {\lambda} =1')
xlabel('Iteratons')
ylabel('Error')
%--------------------------------------------------------------------------
