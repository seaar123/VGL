
clear
clc
close all
%loading a max_run
filenammm0=['Parameters_05_X' num2str(1) '.mat'];
load(filenammm0,'s_Theta2_I2_a','s_Theta2_I1_a','s_Theta2_I2_c','s_Theta2_I1_c','NumOutVars_a','NumRules_a','NumInTerms_a','NumInVars_c','NumInTerms_c','NumOutVars_c','x0d0_S_c','NumRules_c','x0d0_S_a','NumInVars_a','x0d0_F_c','x0d0_F_a');

%% ciritc rfnn
%-------initial states x(t)
%initial value for linear velocity
S=-200;%1;0
%initial value for angular velocity
F=300;%-1;
%-------- States

% inital states range

    x = S:0.1:F;

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
    
    plot(x*0.1,y1(row,:)*0.1,'r','LineWidth',2)
    hold on
    
    plot(x*0.1,yy1(row,:)*0.1,'g','LineWidth',2)
end
legend('Initial G^{1:m}_1','Final G^{1:m}_1')
grid on
title({'Initial and Final Learned Gaussian Membership Functions (GMFs) of ';' The Critic Netwrok of RNF-VGL(\lambda) for One of The Runs ';' ';'GMFs for The Linear Velocity Input State to The Critic Network'})
xlabel('Linear Velocity Universes of Discourse')
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
    
    plot(x*0.1,y2(row,:)*0.1,'r','LineWidth',2)
    hold on
    
    plot(x*0.1,yy2(row,:)*0.1,'g','LineWidth',2)
end
legend('Initial G^{1:m}_2','Final G^{1:m}_2')
grid on
title({'GMFs for The Angular Velocity Input State to The Critic Network'})
xlabel('Angular Velocity Universes of Discourse')
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
    
    plot(x*0.1,ya1(row,:)*0.1,'r','LineWidth',2)
    hold on
    
    plot(x*0.1,yya1(row,:)*0.1,'g','LineWidth',2)
end
legend('Initial G^{1:m}_1','Final G^{1:m}_1')
grid on
title({'Initial and Final Learned Gaussian Membership Functions (GMFs) of ';' The Actor Netwrok of RNF-VGL(\lambda) for One of The Runs ';' ';'GMFs for The Linear Velocity Input State to The Actor Network'})
xlabel('Linear Velocity Universes of Discourse')
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
    
    plot(x*0.1,ya2(row,:)*0.1,'r','LineWidth',2)
    hold on
    
    plot(x*0.1,yya2(row,:)*0.1,'g','LineWidth',2)
end
legend('Initial G^{1:m}_1','Final G^{1:m}_1')
grid on
title({'GMFs for The Angular Velocity Input State to The Actor Network'})
xlabel('Angular Velocity Universes of Discourse')
ylabel('GMFs')





%% drowing weights
% figure
% %weight for output 1
% [nn,mm,kk]=size(s_W_O1_c);
% Bd1 = reshape(s_W_O1_c,[nn mm*kk]);
% subplot(2,1,1)
% plot(Bd1','LineWidth',1)
% grid on
% title({'Deviation of Weights from Initialized Values of Critic Netwrok ';' of RFNN-VGL(\lambda) for One of The Runs ';' ';'Connected Output Layer Weights (W_{R,1}) to g(x_1) of The Critic Network'})
% xlabel('Total Training Steps')
% ylabel( 'W_{R,1}')



%weight for output 2
% [nn,mm,kk]=size(s_W_O2_c);
% Bd2 = reshape(s_W_O2_c,[nn mm*kk]);
% subplot(2,1,2)
% plot(Bd2','LineWidth',1)
% grid on
% title({'Connected Output Layer Weights (W_{R,2}) to g(x_2) of The Critic Network'})
% xlabel('Total Training Steps')
% ylabel('W_{R,2}')



figure
%theta for input 1
[nn,mm,kk]=size( s_Theta2_I1_c);
Bd3 = reshape( s_Theta2_I1_c*0.1,[nn mm*kk]);
subplot(2,1,1)
plot(Bd3','LineWidth',2)
grid on
title({'Deviation of \theta from Initialized Values of The Critic ';' Netwrok of RNF-VGL(\lambda) for One of The Runs ';' ';'Recurrent Parampters for Frist State of The Critic Network'})
xlabel('Total Training Steps')
ylabel( '\theta^{1:m}_1')

%theta for input 2
[nn,mm,kk]=size( s_Theta2_I2_c);
Bd4 = reshape( s_Theta2_I2_c*0.1,[nn mm*kk]);
subplot(2,1,2)
plot(Bd4','LineWidth',2)
grid on
title({' Recurrent Parampters for Second State of The Critic Network'})
xlabel('Total Training Steps')
ylabel( '\theta^{1:m}_2')


% %%
% figure
% % storing weights for actor
% [nn,mm,kk]=size(s_W_O1_a);
% Bd11 = reshape(s_W_O1_a,[nn mm*kk]);
% subplot(2,1,1)
% plot(Bd11','LineWidth',1)
% grid on
% title({'Deviation of Weights from Initialized Values of The Actor Netwrok ';' of RFNN-VGL(\lambda) for One of The Runs ';' ';'Connected Output Layer Weights (W_{R,1}) to g(x_1) of The Actor Network'})
% xlabel('Total Training Steps')
% ylabel( 'W_{R,1}')
% 
% %weight for output 2
% [nn,mm,kk]=size(s_W_O2_a);
% Bd21 = reshape(s_W_O2_a,[nn mm*kk]);
% subplot(2,1,2)
% plot(Bd21','LineWidth',1)
% grid on
% title({'Connected Output Layer Weights (W_{R,2}) to g(x_2) of The Actor Network'})
% xlabel('Total Training Steps')
% ylabel('W_{R,2}')


figure
%theta for input 1
[nn,mm,kk]=size( s_Theta2_I1_a);
Bd31 = reshape( s_Theta2_I1_a*0.1,[nn mm*kk]);
subplot(2,1,1)
plot(Bd31','LineWidth',2)
grid on
title({'Deviation of \theta from Initialized Values of The Actor ';' Netwrok of RNF-VGL(\lambda) for One of The Runs ';' ';'Recurrent Parampters for Frist State of The Actor Network'})
xlabel('Total Training Steps')
ylabel( '\theta^{1:m}_1')


%theta for input 2
[nn,mm,kk]=size( s_Theta2_I2_a);
Bd41 = reshape( s_Theta2_I2_a*0.1,[nn mm*kk]);
subplot(2,1,2)
plot(Bd41','LineWidth',2)
grid on
title({' Recurrent Parampters for Second State of The Actor Network'})
xlabel('Total Training Steps')
ylabel( '\theta^{1:m}_2')
