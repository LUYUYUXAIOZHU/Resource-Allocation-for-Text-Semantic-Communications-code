clc
clear

tic  % 统一计时，统计总仿真时间

%% -------------------------- 公共参数定义（两模型共用/合并）--------------------------
n_devices = 5;          % 用户数量（两模型一致）
radius = 500;           % 小区半径（两模型一致）
p_noise = 180000 * 10^(-17.4); % 噪声功率（毫瓦），约-121dBm（两模型一致）
shadow_factor = 6;      % 阴影衰落因子（两模型一致）
p = 10;                 % 发射功率（dBm）（两模型一致）
snr_range = -10:1:20;   % DeepSC性能表的SNR范围（两模型一致）
load('sem_table.mat');  % DeepSC性能表：行--n_sym; 列--snr_range（两模型共用，仅加载一次）
mento = 1000;           % 仿真次数（两模型一致）
f_th = 0.9;             % 语义相似度阈值（两模型一致）
SE_th = 1;              % 传统模型频谱效率阈值（1bit/s/Hz，仅传统模型用）
t_factor = 40;          % 转换因子（仅Proposed模型用）
s_th_proposed = 1/t_factor; % Proposed模型语义频谱效率阈值（仅Proposed模型用）
s_th_conv = SE_th/t_factor; % 传统模型语义频谱效率阈值（仅传统模型用）

% 各模型特有参数
n_sym_proposed = 1:1:20;       % Proposed模型：传输语义符号数量范围
sym_list_conv = [3, 5, 7, 9, 11]; % 传统模型：传输语义符号数量列表
colors = ['b', 'r', 'g', 'm', 'k', 'c']; % 绘图颜色（新增，用于区分两模型曲线）


%% -------------------------- 1. Proposed Model 计算 --------------------------
% 结果存储：Proposed模型的平均S-SE（信道数量1-10）
SE_results_proposed = []; 

for n_channels = 1:1:10 % 遍历信道数量（1-10）
    SE_mento_proposed = []; % 单次仿真的最优总S-SE（存储1000次仿真结果）
    
    for sim_times = 1:1:mento % 遍历仿真次数（1000次）
%% 大尺度衰落计算（路径损耗+阴影衰落）
d = zeros(n_devices, 1);          % 用户到基站距离
h_large_scale = zeros(n_devices, 1); % 大尺度信道增益
position = zeros(n_devices, 2);   % 用户位置（基站在(0,0)）

radius_dev = radius * sqrt(rand(n_devices, 1)); 
phase = rand(n_devices, 1) * 2 * pi; 
position(:, 1) = radius_dev.*cos(phase);  % x坐标
position(:, 2) = radius_dev.*sin(phase);  % y坐标

for i = 1:n_devices
    d(i) = sqrt(position(i,1)^2 + position(i,2)^2);
    pl = 128.1 + 37.6 * log10(d(i)/1000); % 路径损耗公式
    h_large_scale(i) = 10^(-(pl + shadow_factor)/10); % 大尺度增益
end

%% 小尺度衰落计算（瑞利衰落）
h_real = randn(n_devices, n_channels);
h_image = randn(n_devices, n_channels);
h_small_scale = (h_real.^2 + h_image.^2)/2; % 小尺度增益

        %% 计算每个用户-信道对的最大S-SE
        SE_proposed = zeros(n_devices, n_channels); % 用户-信道对的S-SE
        sym_results_proposed = zeros(n_devices, n_channels); % 最优语义符号数
        f_n_results_proposed = zeros(n_devices, n_channels); % 对应语义相似度
        snr_temp_proposed = [];

        % 遍历每个用户
        for i = 1:n_devices
            h_large = h_large_scale(i);
            % 遍历每个信道
            for j = 1:n_channels
                h_small = h_small_scale(i, j);
                % 遍历发射功率（仅1个值，保留原循环逻辑）
                for p_index = 1:length(p)
                    % 计算SNR
                    snr = 10*log10(10^(p(p_index)/10) * h_large * h_small / p_noise);
                    snr_temp_proposed = [snr_temp_proposed, snr];
                    
                    % SNR范围限制（避免超出sem_table范围）
                    if snr < min(snr_range)
                        snr = min(snr_range);
                    elseif snr > max(snr_range)
                        snr = max(snr_range);
                    end

                    % 遍历语义符号数量，找最大S-SE
                    SE_temp_proposed = zeros(1, length(n_sym_proposed));
                    f_temp_proposed = zeros(1, length(n_sym_proposed));
                    for sym_index = 1:length(n_sym_proposed)
                        snr_index = round(snr) - min(snr_range) + 1;
                        f_n = sem_table(n_sym_proposed(sym_index), snr_index);
                        
                        % 语义相似度阈值判断
                        if f_n < f_th
                            f_n = 0;
                        end
                        f_temp_proposed(p_index, sym_index) = f_n;
                        
                        % 计算S-SE并判断阈值
                        SE_n = f_n / n_sym_proposed(sym_index);
                        if SE_n < s_th_proposed
                            SE_n = 0;
                        end
                        SE_temp_proposed(p_index, sym_index) = SE_n;
                    end
                end

                % 提取当前用户-信道对的最大S-SE
                [max_se, sym_index_opt] = max(SE_temp_proposed);
                SE_proposed(i, j) = max_se;
                f_n_results_proposed(i, j) = f_temp_proposed(1, sym_index_opt);
                sym_results_proposed(i, j) = n_sym_proposed(sym_index_opt);
            end
        end

        % 匈牙利算法分配信道，求总S-SE
        [alpha_proposed, cost_proposed] = Hungarian(-SE_proposed);
        SE_mento_proposed = [SE_mento_proposed, -cost_proposed];
    end

    % 计算当前信道数量的平均S-SE
    SE_results_proposed = [SE_results_proposed, mean(SE_mento_proposed)];
end


%% -------------------------- 2. Conventional Model 计算 --------------------------
% 结果存储：传统模型（不同sym值的S-SE，cell数组）
SE_results_all_conv = cell(1, length(sym_list_conv));
legend_labels_conv = cell(1, length(sym_list_conv));

% 遍历传统模型的语义符号数量
for sym_idx = 1:length(sym_list_conv)
    sym_conv = sym_list_conv(sym_idx);
    fprintf('Conventional Model: 正在仿真 sym = %d\n', sym_conv);
    
    SE_results_conv = []; % 当前sym值的S-SE（信道数量1-10）

    for n_channels = 1:1:10 % 遍历信道数量（1-10）
        SE_mento_conv = []; % 单次仿真的总S-SE（1000次仿真）
        
        for sim_times = 1:1:mento % 遍历仿真次数（1000次）
            %% 大尺度衰落计算（与Proposed模型一致）
            d_conv = zeros(n_devices, 1);
            h_large_scale_conv = zeros(n_devices, 1);
            position_conv = zeros(n_devices, 2);
            
            radius_dev_conv = radius * sqrt(rand(n_devices, 1));
            phase_conv = rand(n_devices, 1) * 2 * pi;
            position_conv(:, 1) = radius_dev_conv.*cos(phase_conv);
            position_conv(:, 2) = radius_dev_conv.*sin(phase_conv);
            
            for i = 1:n_devices
                d_conv(i) = sqrt(position_conv(i,1)^2 + position_conv(i,2)^2);
                pl_conv = 128.1 + 37.6 * log10(d_conv(i)/1000);
                h_large_scale_conv(i) = 10^(-(pl_conv + shadow_factor)/10);
            end

            %% 小尺度衰落计算（与Proposed模型一致）
            h_real_conv = randn(n_devices, n_channels);
            h_image_conv = randn(n_devices, n_channels);
            h_small_scale_conv = (h_real_conv.^2 + h_image_conv.^2)/2;

            %% 计算传统模型的频谱效率（香农公式）
            SE_conv = zeros(n_devices, n_channels); % 传统频谱效率
            for i = 1:n_devices
                h_large_conv = h_large_scale_conv(i);
                for j = 1:n_channels
                    h_small_conv = h_small_scale_conv(i, j);
                    snr_temp_conv = 10^(p/10) * (h_large_conv * h_small_conv) / p_noise;
                    SE_conv(i, j) = log2(1 + snr_temp_conv); % 香农频谱效率
                    
                    % 频谱效率阈值判断
                    if SE_conv(i, j) < SE_th
                        SE_conv(i, j) = 0;
                    end
                end
            end

            % 匈牙利算法分配信道（基于传统频谱效率）
            [alpha_conv, cost_conv] = Hungarian(-SE_conv);

            %% 计算传统模型的语义频谱效率（S-SE）
            f_n_temp_conv = zeros(n_devices, n_channels);
            SE_n_temp_conv = zeros(n_devices, n_channels);
            
            for i = 1:n_devices
                for j = 1:n_channels
                    if alpha_conv(i, j) == 1 % 仅计算分配到的信道
                        % 计算SNR
                        snr_conv = 10*log10(10^(p/10) * h_large_scale_conv(i) * h_small_scale_conv(i, j) / p_noise);
                        
                        % SNR范围限制
                        if snr_conv < min(snr_range)
                            snr_conv = min(snr_range);
                        elseif snr_conv > max(snr_range)
                            snr_conv = max(snr_range);
                        end

                        % 查表获取语义相似度
                        snr_index_conv = round(snr_conv) - min(snr_range) + 1;
                        f_n_conv = sem_table(sym_conv, snr_index_conv);
                        
                        % 语义相似度阈值判断
                        if f_n_conv < f_th
                            f_n_conv = 0;
                        end
                        f_n_temp_conv(i, j) = f_n_conv;
                        
                        % 计算S-SE并判断阈值
                        SE_n_conv = f_n_conv / sym_conv;
                        if SE_n_conv < s_th_conv
                            SE_n_conv = 0;
                        end
                        SE_n_temp_conv(i, j) = SE_n_conv;
                    end
                end
            end

            % 存储当前仿真的总S-SE
            SE_mento_conv = [SE_mento_conv, sum(sum(SE_n_temp_conv))];
        end

        % 计算当前信道数量的平均S-SE
        SE_results_conv = [SE_results_conv, mean(SE_mento_conv)];
    end

    % 保存当前sym值的结果与图例
    SE_results_all_conv{sym_idx} = SE_results_conv;
    legend_labels_conv{sym_idx} = sprintf('Conventional, k_n = %d', sym_conv);
end


%% -------------------------- 3. 合并绘图（两模型结果同图对比） --------------------------
figure; 
hold on; grid on;
% 绘制Proposed模型曲线
plot(SE_results_proposed, [colors(1), '-o'], 'LineWidth', 1.5, 'DisplayName', 'Proposed Model');
% 绘制传统模型各sym值曲线
for sym_idx = 1:length(sym_list_conv)
    plot(SE_results_all_conv{sym_idx}, [colors(sym_idx+1), '-o'], 'LineWidth', 1.5, 'DisplayName', legend_labels_conv{sym_idx});
end
% 图形设置
xlabel('Number of channels, M', 'FontSize', 10);
ylabel('S-SE,Φ(suts/s/Hz)', 'FontSize', 10);
title('The S-SE of the semantic-aware network with different models', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 7);
hold off;

toc  % 结束总计时