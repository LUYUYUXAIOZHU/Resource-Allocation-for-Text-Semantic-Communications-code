clc
clear

n_devices = 5;          % 用户数
n_channels = 5;         % 信道数
radius = 500;           % 小区半径 (m)
p_noise = 180000 * 10^(-17.4); % 噪声功率 (mW)
shadow_factor = 6;      % 阴影衰落因子
p_range = -40:2:30;     % 发射功率范围 (dBm)
mento = 1000;            % 蒙特卡洛仿真次数

% 语义通信参数
n_sym = 1:20;           % 每词语义符号数
snr_range = -10:1:20;   % SNR 范围
load('sem_table.mat');  % 加载 DeepSC 性能表
f_th = 0.9;             % 语义相似度阈值
t_factor = 40;          % 转换因子
s_th = 1 / t_factor;    % S-SE 阈值

% 4G 系统参数
se_actual_RE = [0.15, 0.23, 0.38, 0.6, 0.88, 1.18, 1.48, 1.91, 2.41, 2.73, 3.32, 3.9, 4.52, 5.12, 5.55];
se_actual = se_actual_RE * 14 / 15; % 转换为 bits/s/Hz
snr_th_4G = 10 * log10((2.^(se_actual * 1.25) - 1) / 0.66);

% 5G 系统参数（假设比 4G 高 20%）
se_actual_5G = se_actual * 1.2;
snr_th_5G = 10 * log10((2.^(se_actual_5G * 1.25) - 1) / 0.66);

% 初始化结果存储
S_SE_semantic = zeros(length(p_range), 1);
S_SE_ideal = zeros(length(p_range), 1);
S_SE_4G = zeros(length(p_range), 1);
S_SE_5G = zeros(length(p_range), 1);

% 主循环：遍历发射功率
for p_idx = 1:length(p_range)
    p_dBm = p_range(p_idx);
    p_mW = 10^(p_dBm / 10);
    
    semantic_SE_mento = [];
    ideal_SE_mento = [];
    g4_SE_mento = [];
    g5_SE_mento = [];
    
    for sim_times = 1:mento
        %% 大尺度衰落：路径损耗和阴影
        d = zeros(n_devices, 1);
        h_large_scale = zeros(n_devices, 1);
        position = zeros(n_devices, 2);
        radius_dev = radius * sqrt(rand(n_devices, 1));
        phase = rand(n_devices, 1) * 2 * pi;
        position(:, 1) = radius_dev .* cos(phase);
        position(:, 2) = radius_dev .* sin(phase);
        for i = 1:n_devices
            d(i) = sqrt(position(i,1)^2 + position(i,2)^2);
            pl = 128.1 + 37.6 * log10(d(i)/1000);
            h_large_scale(i) = 10^(-(pl + shadow_factor)/10);
        end

        %% 小尺度衰落：瑞利衰落
        h_real = randn(n_devices, n_channels);
        h_image = randn(n_devices, n_channels);
        h_small_scale = (h_real.^2 + h_image.^2) / 2;

        %% 语义通信系统 S-SE 计算
        SE_semantic = zeros(n_devices, n_channels);
        for i = 1:n_devices
            h_large = h_large_scale(i);
            for j = 1:n_channels
                h_small = h_small_scale(i, j);
                snr = 10 * log10(p_mW * h_large * h_small / p_noise);
                if snr < min(snr_range), snr = min(snr_range); end
                if snr > max(snr_range), snr = max(snr_range); end
                snr_index = round(snr) - min(snr_range) + 1;
                SE_temp = zeros(1, length(n_sym));
                for sym_index = 1:length(n_sym)
                    f_n = sem_table(n_sym(sym_index), snr_index);
                    if f_n < f_th, f_n = 0; end
                    SE_n = f_n / n_sym(sym_index);
                    if SE_n < s_th, SE_n = 0; end
                    SE_temp(sym_index) = SE_n;
                end
                SE_semantic(i, j) = max(SE_temp);
            end
        end
        [~, cost] = Hungarian(-SE_semantic);
        semantic_SE_mento = [semantic_SE_mento, -cost];

        %% 理想系统 S-SE 计算
        SE_ideal = zeros(n_devices, n_channels);
        for i = 1:n_devices
            h_large = h_large_scale(i);
            for j = 1:n_channels
                h_small = h_small_scale(i, j);
                snr_linear = p_mW * h_large * h_small / p_noise;
                SE_bits = log2(1 + snr_linear);
                SE_ideal(i, j) = SE_bits / t_factor;
            end
        end
        [~, cost] = Hungarian(-SE_ideal);
        ideal_SE_mento = [ideal_SE_mento, -cost];

        %% 4G 系统 S-SE 计算
        SE_4G = zeros(n_devices, n_channels);
        for i = 1:n_devices
            h_large = h_large_scale(i);
            for j = 1:n_channels
                h_small = h_small_scale(i, j);
                snr_temp = 10 * log10(p_mW * h_large * h_small / p_noise);
                if snr_temp < snr_th_4G(1)
                    SE_4G(i, j) = 0;
                elseif snr_temp > snr_th_4G(end)
                    SE_4G(i, j) = se_actual(end) / t_factor;
                else
                    idx = find(snr_temp >= snr_th_4G, 1, 'last');
                    SE_4G(i, j) = se_actual(idx) / t_factor;
                end
            end
        end
        [~, cost] = Hungarian(-SE_4G);
        g4_SE_mento = [g4_SE_mento, -cost];

        %% 5G 系统 S-SE 计算
        SE_5G = zeros(n_devices, n_channels);
        for i = 1:n_devices
            h_large = h_large_scale(i);
            for j = 1:n_channels
                h_small = h_small_scale(i, j);
                snr_temp = 10 * log10(p_mW * h_large * h_small / p_noise);
                if snr_temp < snr_th_5G(1)
                    SE_5G(i, j) = 0;
                elseif snr_temp > snr_th_5G(end)
                    SE_5G(i, j) = se_actual_5G(end) / t_factor;
                else
                    idx = find(snr_temp >= snr_th_5G, 1, 'last');
                    SE_5G(i, j) = se_actual_5G(idx) / t_factor;
                end
            end
        end
        [~, cost] = Hungarian(-SE_5G);
        g5_SE_mento = [g5_SE_mento, -cost];
    end

    % 平均 S-SE
    S_SE_semantic(p_idx) = mean(semantic_SE_mento);
    S_SE_ideal(p_idx) = mean(ideal_SE_mento);
    S_SE_4G(p_idx) = mean(g4_SE_mento);
    S_SE_5G(p_idx) = mean(g5_SE_mento);
end

%% 绘图
figure;
plot(p_range, S_SE_semantic, 'b-o', 'LineWidth', 2, 'DisplayName', 'Semantic');
hold on;
plot(p_range, S_SE_ideal, 'r-o', 'LineWidth', 2, 'DisplayName', 'Ideal');
plot(p_range, S_SE_4G, 'g-o', 'LineWidth', 2, 'DisplayName', '4G');
plot(p_range, S_SE_5G, 'm-o', 'LineWidth', 2, 'DisplayName', '5G');
xlabel('Transmit Power,pn (dBm)');
ylabel('S-SE,Φ(suts/s/Hz)');
grid on;
legend('show');
title('The S-SE versus the transmit power');
