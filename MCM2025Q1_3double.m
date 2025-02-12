clear, clc;

% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3; % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.2; % 台阶初始高度 (m)

% 台阶使用年限和踩踏频率
T = 300 * 365; % 台阶年限为300年
f = 10; % 每天的踩踏频率

% 脚踩区域尺寸（椭圆参数）
L_foot = 0.25; % 脚的长度 (m)
W_foot = 0.12; % 脚的宽度 (m)
a = L_foot / 2; % 椭圆的长轴半径
b = W_foot / 2; % 椭圆的短轴半径

% 材料硬度
H = 5e10; % 硬度 (Pa)

% 摩擦系数范围
k_min = 0.02; % 脚尖处摩擦系数
k_max = 0.08; % 脚跟处摩擦系数

% 压力和滑动距离
P = 30000; % 压力 (N/m^2)
Delta_s = 0.01; % 滑动距离 (m)

% 高斯分布参数
% x方向双高斯分布
mu_x1 = -0.2; % 双高斯路径1中心
mu_x2 = 0.2; % 双高斯路径2中心
sigma_x1 = 0.1; % 路径1的标准差
sigma_x2 = 0.1; % 路径2的标准差
weight1 = 0.5; %  权重

% y方向单高斯分布
mu_y = -Ly / 2 + L_foot; % 高斯分布中心 (y方向)
sigma_y = 0.03; % 高斯分布标准差 (y方向)

% 上下行比例
alpha_values = [1, 0]; % 上行比例数组

% 网格划分
Nx = 100; % x方向网格数
Ny = 100; % y方向网格数
y = linspace(-Ly/2, Ly/2, Ny); % 台阶宽度方向 -> 脚的长度方向
x = linspace(-Lx/2, Lx/2, Nx); % 台阶长度方向 -> 脚的宽度方向
[X, Y] = meshgrid(x, y);

% 初始化台阶高度矩阵
H_step = Lz * ones(Ny, Nx); % 台阶初始高度

% 模拟脚踩行为
num_steps = 1000; % 脚踩次数

figure('Position', [100, 100, 1200, 800]);
for idx = 1:length(alpha_values)
    alpha = alpha_values(idx); % 当前上行比例
    H_step = Lz * ones(Ny, Nx); % 重置台阶初始高度

    for step = 1:num_steps
        % x方向双高斯分布
        if rand < weight1
            foot_center_x = normrnd(mu_x1, sigma_x1); % 路径1位置
        else
            foot_center_x = normrnd(mu_x2, sigma_x2); % 路径2位置
        end

        % y方向单高斯分布
        foot_center_y = normrnd(mu_y, sigma_y); % y方向位置

        % 计算相对位置和椭圆掩码
        X_rel = X - foot_center_x; % x方向相对位置
        Y_rel = Y - foot_center_y; % y方向相对位置
        ellipse_mask = (X_rel.^2 / b^2 + Y_rel.^2 / a^2) <= 1; % 椭圆区域掩码

        % 计算摩擦系数（非线性变化）
        f_y = (Y_rel + a) / (2 * a); % y方向归一化
        k_up = k_min + (k_max - k_min) * f_y; % 上行摩擦系数
        k_down = k_max + (k_min - k_max) * f_y; % 下行摩擦系数

        % 计算磨损量
        wear = k_up .* P * Delta_s * alpha * f * T / H + ...
               k_down .* P * Delta_s * (1 - alpha) * f * T / H;

        % 更新台阶高度
        H_step = H_step - wear .* ellipse_mask;
    end

    % 确保台阶高度非负
    H_step = max(H_step, 0);

    min_val = min(H_step(:));
    max_val = max(H_step(:));
    c_range = [min_val, max_val];

       % 自定义颜色图，将最高值设置为大理石颜色
    custom_colormap = parula(256);
    marble_color = [0.9, 0.9, 0.9]; % 浅灰色
    custom_colormap(end, :) = marble_color;

    % 绘制凹陷效果的3D图
    subplot(2, 2, idx);
    surf(X, Y, H_step);
    colorbar;
    shading interp;
    colormap(custom_colormap);
    title(['Upward Ratio: ', num2str(alpha)], 'FontSize', 14);
    xlabel('Step Length (m)');
    ylabel('Step Width (m)');
    zlabel('Step Height (m)');
    axis tight;
    axis equal;
    zlim([0, Lz]);
    view(135, 30);
    
    subplot(2, 2, idx+2);
    surf(X, Y, H_step);
    colorbar;
    shading interp;
    colormap(custom_colormap);
    title(['view xOy'], 'FontSize', 14);
    xlabel('Step Length (m)');
    ylabel('Step Width (m)');
    zlabel('Step Height (m)');
    axis tight;
    axis equal;
    zlim([0, Lz]);
    view(0, 270);
end
