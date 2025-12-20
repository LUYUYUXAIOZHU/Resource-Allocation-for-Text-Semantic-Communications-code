function fig4a_channels()
clc; clear; 

% 参数设置
n_devices = 5;
radius = 500;
p_noise = 180000 * 10^(-17.4);
shadow_factor = 6;
p = 10; % dBm
mento = 1000;

% 语义通信参数
n_sym = 1:1:20;
snr_range = -10:1:20;
load('sem_table.mat');
f_th = 0.9;
t_factor = 40;
s_th = 1/t_factor;

% 存储结果
SE_semantic = [];
SE_ideal = [];
SE_4G = [];
SE_5G = [];

fprintf('开始计算图4(a)...\n');

for n_channels = 1:15
    fprintf('计算信道数: %d/10\n', n_channels);
    
    % 语义通信系统
    SE_mento_sem = [];
    for sim_times = 1:mento
        [H_large, H_small] = generate_channels(n_devices, n_channels, radius, shadow_factor);
        SE_matrix = calculate_semantic_SE(n_devices, n_channels, H_large, H_small, p, p_noise, ...
                                         n_sym, snr_range, sem_table, f_th, s_th);
        [~, cost] = Hungarian(-SE_matrix);
        SE_mento_sem = [SE_mento_sem, -cost];
    end
    SE_semantic = [SE_semantic, mean(SE_mento_sem)];
    
    % 理想系统
    SE_mento_ideal = [];
    for sim_times = 1:mento
        [H_large, H_small] = generate_channels(n_devices, n_channels, radius, shadow_factor);
        SE_matrix = calculate_ideal_SE(n_devices, n_channels, H_large, H_small, p, p_noise);
        [~, cost] = Hungarian(-SE_matrix);
        SE_mento_ideal = [SE_mento_ideal, -cost/t_factor];
    end
    SE_ideal = [SE_ideal, mean(SE_mento_ideal)];
    
    % 4G系统
    SE_mento_4G = [];
    for sim_times = 1:mento
        [H_large, H_small] = generate_channels(n_devices, n_channels, radius, shadow_factor);
        SE_matrix = calculate_4G_SE(n_devices, n_channels, H_large, H_small, p, p_noise);
        [~, cost] = Hungarian(-SE_matrix);
        SE_mento_4G = [SE_mento_4G, -cost/t_factor];
    end
    SE_4G = [SE_4G, mean(SE_mento_4G)];
    
    % 5G系统
    SE_mento_5G = [];
    for sim_times = 1:mento
        [H_large, H_small] = generate_channels(n_devices, n_channels, radius, shadow_factor);
        SE_matrix = calculate_5G_SE(n_devices, n_channels, H_large, H_small, p, p_noise);
        [~, cost] = Hungarian(-SE_matrix);
        SE_mento_5G = [SE_mento_5G, -cost/t_factor];
    end
    SE_5G = [SE_5G, mean(SE_mento_5G)];
end

% 绘图
figure('Position', [100, 100, 800, 600]);
plot(1:15, SE_semantic, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Semantic');
hold on;
plot(1:15, SE_ideal, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Ideal');
plot(1:15, SE_4G, 'g-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '4G');
plot(1:15, SE_5G, 'm-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '5G');

xlabel('Number of Channels,M');
ylabel('S-SE,Φ(suts/s/Hz)');
title('The S-SE versus the number of channels','FontSize',12);
legend('show', 'Location', 'northwest');
grid on;

% 保存图像
saveas(gcf, 'Fig4a_Channels.png');
fprintf('图4(a)已保存为 Fig4a_Channels.png\n');

end
function SE_matrix = calculate_semantic_SE(n_devices, n_channels, H_large, H_small, p, p_noise, ...
                                          n_sym, snr_range, sem_table, f_th, s_th)
% 计算语义SE矩阵
SE_matrix = zeros(n_devices, n_channels);

for i = 1:n_devices
    for j = 1:n_channels
        snr = calculate_SNR(H_large(i), H_small(i,j), p, p_noise);
        
        % 穷举搜索最优语义符号数
        best_SE = 0;
        for k_idx = 1:length(n_sym)
            k = n_sym(k_idx);
            similarity = get_semantic_similarity(k, snr, snr_range, sem_table, f_th);
            SE_temp = similarity / k;
            if SE_temp >= s_th && SE_temp > best_SE
                best_SE = SE_temp;
            end
        end
        SE_matrix(i,j) = best_SE;
    end
end
end

function SE_matrix = calculate_ideal_SE(n_devices, n_channels, H_large, H_small, p, p_noise)
% 计算理想系统SE矩阵
SE_matrix = zeros(n_devices, n_channels);
SE_th = 1; % 1 bit/s/Hz

for i = 1:n_devices
    for j = 1:n_channels
        snr_linear = 10^(p/10) * H_large(i) * H_small(i,j) / p_noise;
        SE = log2(1 + snr_linear);
        if SE < SE_th
            SE = 0;
        end
        SE_matrix(i,j) = SE;
    end
end
end

function SE_matrix = calculate_4G_SE(n_devices, n_channels, H_large, H_small, p, p_noise)
% 计算4G系统SE矩阵
SE_matrix = zeros(n_devices, n_channels);
SE_th = 1; % SE阈值

% 4G系统的SNR-CQI映射表
se_actual_RE = [0.15, 0.23, 0.38, 0.6, 0.88, 1.18, 1.48, 1.91, 2.41, 2.73, 3.32, 3.9, 4.52, 5.12, 5.55];
se_actual = se_actual_RE * 14 / 15; % 转换为bits/s/Hz
snr_th = 10*log10((2.^(se_actual*1.25)-1)/0.66);

for i = 1:n_devices
    for j = 1:n_channels
        snr_temp = 10^(p/10) * (H_large(i) * H_small(i,j)) / p_noise;
        snr_db = 10 * log10(snr_temp);
        
        if snr_db < snr_th(1)
            SE_matrix(i,j) = 0;
        elseif snr_db > snr_th(end)
            SE_matrix(i,j) = se_actual(end);
        else
            idx = find(snr_db >= snr_th, 1, 'last');
            SE_matrix(i,j) = se_actual(idx);
        end
        
        if SE_matrix(i,j) < SE_th
            SE_matrix(i,j) = 0;
        end
    end
end
end

function SE_matrix = calculate_5G_SE(n_devices, n_channels, H_large, H_small, p, p_noise)
% 计算5G系统SE矩阵
SE_matrix = zeros(n_devices, n_channels);
SE_th = 1; % SE阈值

% 5G系统的SNR-CQI映射表

se_actual_RE = [ 0.1523, 0.3770, 0.8770, 1.4766, 1.9141, 2.4063, 2.7305, ...
             3.3223, 3.9023, 4.5234, 5.1152, 5.5547, 6.2266, 6.9141, 7.4063];

se_actual = se_actual_RE * 14 / 15 ;
snr_th = 10*log10((2.^(se_actual*1.25)-1)/0.66);

for i = 1:n_devices
    for j = 1:n_channels
        snr_temp = 10^(p/10) * (H_large(i) * H_small(i,j)) / p_noise;
        snr_db = 10 * log10(snr_temp);
        
        if snr_db < snr_th(1)
            SE_matrix(i,j) = 0;
        elseif snr_db > snr_th(end)
            SE_matrix(i,j) = se_actual(end);
        else
            idx = find(snr_db >= snr_th, 1, 'last');
            SE_matrix(i,j) = se_actual(idx);
        end
        
        if SE_matrix(i,j) < SE_th
            SE_matrix(i,j) = 0;
        end
    end
end
end

function [H_large, H_small] = generate_channels(n_devices, n_channels, radius, shadow_factor)
% 生成大尺度和小尺度信道
H_large = zeros(n_devices, 1);
H_small = zeros(n_devices, n_channels);

% 生成用户位置
position = zeros(n_devices, 2);
radius_dev = radius * sqrt(rand(n_devices, 1));
phase = rand(n_devices, 1) * 2 * pi;
position(:, 1) = radius_dev .* cos(phase);
position(:, 2) = radius_dev .* sin(phase);

% 大尺度衰落
for i = 1:n_devices
    d = sqrt(position(i,1)^2 + position(i,2)^2);
    pl = 128.1 + 37.6 * log10(d/1000);
    H_large(i) = 10^(-(pl + shadow_factor)/10);
end

% 小尺度衰落（瑞利衰落）
h_real = randn(n_devices, n_channels);
h_image = randn(n_devices, n_channels);
H_small = (h_real.^2 + h_image.^2) / 2;
end

function snr = calculate_SNR(h_large, h_small, p, p_noise)
% 计算SNR
snr_linear = 10^(p/10) * h_large * h_small / p_noise;
snr = 10 * log10(snr_linear);
end

function similarity = get_semantic_similarity(k, snr, snr_range, sem_table, f_th)
% 获取语义相似度
if snr < min(snr_range)
    snr = min(snr_range);
elseif snr > max(snr_range)
    snr = max(snr_range);
end

snr_index = round(snr) - min(snr_range) + 1;
if snr_index < 1
    snr_index = 1;
elseif snr_index > size(sem_table, 2)
    snr_index = size(sem_table, 2);
end

similarity = sem_table(k, snr_index);

if similarity < f_th
    similarity = 0;
end
end