function fig4c_transform()
clc; clear;

% 压缩比到变换因子的转换
compression_ratio = 0.30:0.05:1;
t_factor = 5*8*compression_ratio; % 变换因子

% 从之前的仿真结果获取基准值
% 这些值需要根据实际仿真结果调整
semantic_baseline = 1.17779;    % 语义通信系统基准S-SE
ideal_baseline = 0.79299*40;    % 理想系统基准
G4_baseline = 0.46323*40;       % 4G系统基准  
G5_baseline = 0.55966*40;       % 5G系统基准

% 计算不同变换因子下的S-SE
semantic = ones(1, length(compression_ratio)) * semantic_baseline;
ideal = ideal_baseline ./ t_factor;
G4 = G4_baseline ./ t_factor;
G5 = G5_baseline ./ t_factor;

% 绘图
figure('Position', [100, 100, 800, 600]);
plot(t_factor, semantic, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Semantic');
hold on;
plot(t_factor, ideal, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Ideal');
plot(t_factor, G4, 'g-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '4G System');
plot(t_factor, G5, 'm-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '5G');

xlabel('Transforming Factor $\mu$ (bits/word)', 'Interpreter', 'latex');
ylabel('S-SE,Φ(suts/s/Hz)');
title('The S-SE versus the transforming factor.');
legend('show', 'Location', 'northeast');
grid on;

% 添加关键点标记
critical_point = 19; % 19 bits/word
yl = ylim;
line([critical_point, critical_point], yl, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1,'DisplayName', 'μ=19');
text(critical_point+0.5, yl(1)+0.1*(yl(2)-yl(1)), '\mu=19', 'FontSize', 10);

critical_point2 = 27; % 27 bits/word
line([critical_point2, critical_point2], yl, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1,'DisplayName', 'μ=27');
text(critical_point2+0.5, yl(1)+0.15*(yl(2)-yl(1)), '\mu=27', 'FontSize', 10);

% 保存图像
saveas(gcf, 'Fig4c_Transform.png');
fprintf('图4(c)已保存为 Fig4c_Transform.png\n');

end