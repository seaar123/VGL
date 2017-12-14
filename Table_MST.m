% DHDP(lambda)

clear;
clc;
close all


% Number of time steps per episode
for run_count=1:1%max_runs
    %---------------------
    % MVGL  (lambda=0)
    filenammm0=['Parameters_01_X' num2str(run_count) '.mat'];
    load(filenammm0,'maximum_trail','Total_MSE');
    
    Total_MSE_runs(1,1:maximum_trail,run_count)=Total_MSE;% mse iterations \all iterations
    
    
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

Mean_All=mean(y_mean1)
Mean_Final=y_mean1(end)