clc;close all;clear;
load("sem_table.mat");

% 1. 准备数据
% 创建网格坐标：x对应列(1-32)，y对应行(1-21)
[X, Y] = meshgrid(-10:20, 1:20);

Z=sem_table;
figure;
surf(X, Y, Z);

% 3. 美化图形
title('The semantic similarity for DeepSC');
xlabel('ξn,m');
ylabel('SNR, γn,n(dB)');
zlabel('ξn,m');

colorbar;          % 保留颜色栏
shading interp;    % 平滑着色
grid on;           % 显示网格

view(35, 30);      % 调整视角，让图形展示更直观
colormap(parula);  % 优化颜色映射（parula色表更柔和）
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman'); % 统一字体和大小
