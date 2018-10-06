clear all;close all;clc

turbineData=importdata('postProcessing/turbines/0/turbine.csv',',',1).data;
time=turbineData(:,1);
cp=turbineData(:,4);

plot(time,cp)
grid on;
title(['tFinal = ' num2str(max(time)) 's']);
xlabel('time'); ylabel('Cp');
pause
