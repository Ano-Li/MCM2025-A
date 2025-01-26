%%
%新版本采用矩阵运算
% 统一colorbar范围
% 清理工作区和设置参数
clear, clc;

% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3; % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.2; % 台阶初始高度 (m)

% 台阶使用年限和踩踏频率
T = 300 * 365; % 台阶年限为300年
f = 1000; % 每天的踩踏频率为5000次

% 脚踩区域尺寸（椭圆参数）
L_foot = 0.25; % 脚的长度 (m)
W_foot = 0.12; % 脚的宽度 (m)
a = L_foot / 2; % 椭圆的长轴半径
b = W_foot / 2; % 椭圆的短轴半径

% 材料硬度
H = 5e10; % 硬度 (Pa)，如大理石

% 摩擦系数范围
k_min = 0.02; % 脚尖处摩擦系数
k_max = 0.08; % 脚跟处摩擦系数

% 压力和滑动距离
P = 30000; % 压力 (N/m^2)
Delta_s = 0.01; % 滑动距离 (m)

% 非线性控制参数
n = 2; % 非线性程度控制

% 上下行比例
alpha_values = [0.7,0.3]; % 上行比例数组

% 网格划分
Nx = 100; % x方向网格数
Ny = 100; % y方向网格数
y = linspace(-Ly/2, Ly/2, Ny); % 台阶宽度方向 -> 脚的长度方向
x = linspace(-Lx/2, Lx/2, Nx); % 台阶长度方向 -> 脚的宽度方向
[X, Y] = meshgrid(x, y);

% 定义脚踩区域中心（现在脚踩位置改为台阶宽度方向边缘）
foot_center_x = 0; % 脚掌中心的x坐标
foot_center_y = -Ly / 2 + L_foot; % 脚掌中心的y坐标，位于台阶上方宽度边缘

% 初始化台阶高度矩阵
H_step = Lz * ones(Ny, Nx); % 台阶初始高度

% 记录所有子图的高度范围
min_height = 0;
max_height = Lz;

% 绘制每个子图
figure('Position', [100, 100, 1200, 800]);
% 预处理每个子图的高度范围
for idx = 1:length(alpha_values)
    alpha = alpha_values(idx); % 当前上行比例
    H_step = Lz * ones(Ny, Nx); % 重置台阶初始高度

    % 计算椭圆脚掌区域内的摩擦系数和磨损分布
    % 计算每个点是否在椭圆区域内
    X_rel = X - foot_center_x; % 计算x方向相对位置
    Y_rel = Y - foot_center_y; % 计算y方向相对位置
    ellipse_mask = (X_rel.^2 / b^2 + Y_rel.^2 / a^2) <= 1; % 计算椭圆区域的掩码

    % 计算摩擦系数（非线性变化）
    f_y = (Y_rel + a) / (2 * a); % y方向的归一化位置
    k_up = k_min + (k_max - k_min) * f_y; % 上行的摩擦系数
    k_down = k_max + (k_min - k_max) * f_y; % 下行的摩擦系数

     k_up = k_up.*ellipse_mask;
     k_down = k_down.*ellipse_mask;

    % 计算磨损量（矩阵运算）
    wear = k_up .* P * Delta_s * alpha*f * T / H + k_down .* P * Delta_s * (1-alpha)*f * T / H; % 每年累积磨损量

    % 更新磨损分布和台阶高度
    W_total = wear .* ellipse_mask; % 只有在椭圆区域内才有磨损
    H_step = H_step - W_total; % 更新台阶高度
    H_step = max(H_step, 0); % 确保高度非负

    % 更新高度范围
    min_height = min(min_height, min(H_step(:)));
    max_height = max(max_height, max(H_step(:)));

    % 自定义颜色图，将最高值设置为大理石颜色
    custom_colormap = parula(256); % 使用默认 parula 颜色图作为基础
    marble_color = [0.9, 0.9, 0.9]; % 大理石的浅灰色 (RGB 值，可根据实际调整)
    custom_colormap(end, :) = marble_color; % 将最高值对应的颜色改为大理石颜色
    
    % 绘制凹陷效果的3D图
    subplot(2, 2, idx);
    surf(X, Y, H_step); % 使用 surf 绘制3D图
    colorbar;
    caxis([min(H_step(:)), max(H_step(:))]); % 设置统一的颜色范围
    shading interp; % 平滑颜色过渡
    colormap(custom_colormap); % 应用自定义颜色图
    title(['Upward Ratio: ', num2str(alpha)]); % 英文标题
    xlabel('Step Length  (m)'); % 对应脚宽方向
    ylabel('Step Width  (m)'); % 对应脚长方向
    zlabel('Step Height (m)'); % 对应台阶高度
    axis tight; % 缩放坐标轴范围至有效数据区域
    axis equal;
    zlim([0, 0.2]); % 限制 z 轴的范围，突出磨损区域
    view(135, 30); % 设置为3D视角

    subplot(2, 2, idx+2);
    surf(X, Y, H_step); % 使用 surf 绘制3D图
    colorbar;
    caxis([min(H_step(:)), max(H_step(:))]); % 设置统一的颜色范围
    shading interp; % 平滑颜色过渡
    colormap(custom_colormap); % 应用自定义颜色图
    title(['Upward Ratio: ', num2str(alpha)]); % 英文标题
    xlabel('Step Length  (m)'); % 对应脚宽方向
    ylabel('Step Width  (m)'); % 对应脚长方向
    zlabel('Step Height (m)'); % 对应台阶高度
    axis tight; % 缩放坐标轴范围至有效数据区域
    axis equal;
    zlim([0, 0.2]); % 限制 z 轴的范围，突出磨损区域
    view(90, 0); % 设置为3D视角
end
disp("计算完成！");