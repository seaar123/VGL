fgfgfg=hfhf
clear
clc
%loading a max_run
filenammm0=['Parameters_05_X' num2str(3) '.mat'];
load(filenammm0);

%% ciritc rfnn
%-------initial states x(t)
%initial value for linear velocity
v=5;%1;0
%initial value for angular velocity
w=5;%-1;
%-------- States
V(1,:)=[v,w];
% inital states range
ll=11;
if V(1,1)<=V(1,2)
    x = ll*V(1,1):0.1:V(1,2)*ll;
else
    x = ll*V(1,2):0.1:V(1,1)*ll;
end

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
title({'Initial and Final Learned Gaussian Membership Functions (GMFs) of Critic Netwrok ';' of RFNN-VGL(\lambda) for One of The Runs ';' ';'Frist Layer GMFs for The Linear Velocity Input State to The Critic Network'})
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




%% actor rfnn

%initial paramters
% Unpack the network's parameters first...
Xt=x0d0_S_a;
NumInVars=NumInVars_a;
NumInTerms=NumInTerms_a;
NumOutVars=NumOutVars_a;
NumRules=NumRules_a;

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
    ya1(i,:) = gaussmf(x,[width_x1(i) center_x1(i)]);
end
% paramters for guassion of x2
center_x2=mean2(2,:);
width_x2=sigma2(2,:);
for i=1: NumInTerms
    ya2(i,:) = gaussmf(x,[width_x2(i) center_x2(i)]);
end

% final parameters
% Unpack the network's parameters first...
Xt=x0d0_F_a;
NumInVars=NumInVars_a;
NumInTerms=NumInTerms_a;
NumOutVars=NumOutVars_a;
NumRules=NumRules_a;

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
    yya1(i,:) = gaussmf(x,[width_x1(i) center_x1(i)]);
end
% GMFs for x1
figure
subplot(2,1,1)
for row = 1 :  NumInTerms
    
    plot(x,ya1(row,:),'r','LineWidth',2)
    hold on
    
    plot(x,yya1(row,:),'g','LineWidth',2)
end
legend('Initial G^{1:m}_1','Final G^{1:m}_1')
grid on
title({'Initial and Final Learned Gaussian Membership Functions (GMFs) of Actor Netwrok ';' of RFNN-VGL(\lambda) for One of The Runs ';' ';' Frist Layer GMFs for The Linear Velocity Input State to The Actor Network'})
xlabel('Linear Velocity Universes of Discourse (U_1)')
ylabel( 'GMFs')






% paramters for guassion of x2
center_x2=mean2(2,:);
width_x2=sigma2(2,:);
for i=1: NumInTerms
    yya2(i,:) = gaussmf(x,[width_x2(i) center_x2(i)]);
end
% GMFs for x2

subplot(2,1,2)
for row = 1 :  NumInTerms
    
    plot(x,ya2(row,:),'r','LineWidth',2)
    hold on
    
    plot(x,yya2(row,:),'g','LineWidth',2)
end
legend('Initial G^{1:m}_1','Final G^{1:m}_1')
grid on
title({'Frist Layer GMFs for The Angular Velocity Input State to The Actor Network'})
xlabel('Angular Velocity Universes of Discourse (U_2)')
ylabel( 'GMFs')





%% drowing weights
figure
%weight for output 1
[nn,mm,kk]=size(s_W_O1_c);
Bd1 = reshape(s_W_O1_c,[nn mm*kk]);
subplot(2,1,1)
plot(Bd1','LineWidth',1)
grid on
title({'Deviation of Weights from Initialized Values of Critic Netwrok ';' of RFNN-VGL(\lambda) for One of The Runs ';' ';'Connected Output Layer Weights (W_{R,1}) to g(x_1) of The Critic Network'})
xlabel('Total Training Steps')
ylabel( 'W_{R,1}')



%weight for output 2
[nn,mm,kk]=size(s_W_O2_c);
Bd2 = reshape(s_W_O2_c,[nn mm*kk]);
subplot(2,1,2)
plot(Bd2','LineWidth',1)
grid on
title({'Connected Output Layer Weights (W_{R,2}) to g(x_2) of The Critic Network'})
xlabel('Total Training Steps')
ylabel('W_{R,2}')



figure
%theta for input 1
[nn,mm,kk]=size( s_Theta2_I1_c);
Bd3 = reshape( s_Theta2_I1_c,[nn mm*kk]);
subplot(2,1,1)
plot(Bd3','LineWidth',2)
grid on
title({'Deviation of \theta from Initialized Values of Critic Netwrok ';' of RFNN-VGL(\lambda) for One of The Runs ';' ';'Recurrent Parampters for Frist Input State for GMFs of The Critic Network'})
xlabel('Total Training Steps')
ylabel( '\theta^{1:m}_1')

%theta for input 2
[nn,mm,kk]=size( s_Theta2_I2_c);
Bd4 = reshape( s_Theta2_I2_c,[nn mm*kk]);
subplot(2,1,2)
plot(Bd4','LineWidth',2)
grid on
title({' Recurrent Parampters for Second Input State for GMFs of The Critic Network'})
xlabel('Total Training Steps')
ylabel( '\theta^{1:m}_2')


%%
figure
% storing weights for actor
[nn,mm,kk]=size(s_W_O1_a);
Bd11 = reshape(s_W_O1_a,[nn mm*kk]);
subplot(2,1,1)
plot(Bd11','LineWidth',1)
grid on
title({'Deviation of Weights from Initialized Values of The Actor Netwrok ';' of RFNN-VGL(\lambda) for One of The Runs ';' ';'Connected Output Layer Weights (W_{R,1}) to g(x_1) of The Actor Network'})
xlabel('Total Training Steps')
ylabel( 'W_{R,1}')

%weight for output 2
[nn,mm,kk]=size(s_W_O2_a);
Bd21 = reshape(s_W_O2_a,[nn mm*kk]);
subplot(2,1,2)
plot(Bd21','LineWidth',1)
grid on
title({'Connected Output Layer Weights (W_{R,2}) to g(x_2) of The Actor Network'})
xlabel('Total Training Steps')
ylabel('W_{R,2}')


figure
%theta for input 1
[nn,mm,kk]=size( s_Theta2_I1_a);
Bd31 = reshape( s_Theta2_I1_a,[nn mm*kk]);
subplot(2,1,1)
plot(Bd31','LineWidth',2)
grid on
title({'Deviation of \theta from Initialized Values of The Actor Netwrok ';' of RFNN-VGL(\lambda) for One of The Runs ';' ';'Recurrent Parampters for Frist Input State for GMFs of The Actor Network'})
xlabel('Total Training Steps')
ylabel( '\theta^{1:m}_1')


%theta for input 2
[nn,mm,kk]=size( s_Theta2_I2_a);
Bd41 = reshape( s_Theta2_I2_a,[nn mm*kk]);
subplot(2,1,2)
plot(Bd41','LineWidth',2)
grid on
title({' Recurrent Parampters for Second Input State for GMFs of The Actor Network'})
xlabel('Total Training Steps')
ylabel( '\theta^{1:m}_2')
