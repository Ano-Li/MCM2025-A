% 定义参数
H = 10;        % 硬度
k = 0.02;      % 磨损系数
p = 100;       % 法向载荷
s = 0.5;       % 步幅

% 定义磨损量范围
W_values = linspace(1, 20, 100); % 磨损量从 1 到 20

% 计算对应的使用频率
f_values = H .* W_values ./ (k .* p .* s);

% 绘制曲线
figure;
plot(W_values, f_values, 'LineWidth', 2);
xlabel('磨损量 W (立方厘米)');
ylabel('使用频率 f (次/年)');
title('使用频率与磨损量的关系');
grid on;
%%
% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3;   % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.15; % 台阶高度 (m)

% 脚踩区域尺寸
L_foot = 0.25; % 脚的长度 (m)
W_foot = 0.12; % 脚的宽度 (m)

% 材料硬度
H = 2e9; % 硬度 (Pa)，如大理石

% 摩擦系数范围
k_min = 0.2; % 脚尖处摩擦系数
k_max = 0.8; % 脚跟处摩擦系数

% 压力和滑动距离
P = 100; % 压力 (N/m^2)
Delta_s = 0.01; % 滑动距离 (m)

% 非线性控制参数
n = 2; % 多项式模型的非线性程度

% 网格划分
Nx = 50; % x方向网格数
Ny = 50; % y方向网格数
y = linspace(-Ly/2, Ly/2, Ny); % 台阶宽度方向 -> 脚的长度方向
x = linspace(-Lx/2, Lx/2, Nx); % 台阶长度方向 -> 脚的宽度方向
[X, Y] = meshgrid(x, y);

% 定义脚踩区域中心
foot_center_x = 0; % 脚踩区域中心的x坐标
foot_center_y = 0; % 脚踩区域中心的y坐标

% 确定脚踩区域范围
foot_x_min = foot_center_x - W_foot / 2;
foot_x_max = foot_center_x + W_foot / 2;
foot_y_min = foot_center_y - L_foot / 2;
foot_y_max = foot_center_y + L_foot / 2;

% 初始化磨损分布矩阵
W_total = zeros(Nx, Ny);

% 计算脚踩区域内的摩擦系数和磨损分布
for ix = 1:Nx
    for iy = 1:Ny
        % 判断当前位置是否在脚踩区域内
        if x(ix) >= foot_x_min && x(ix) <= foot_x_max && ...
           y(iy) >= foot_y_min && y(iy) <= foot_y_max
            % 非线性摩擦系数分布（基于y方向）
            f_y = (y(iy) - foot_y_min) / (foot_y_max - foot_y_min); % 归一化位置
            k_y = k_min + (k_max - k_min) * f_y; % 多项式模型

            % 计算磨损量
            W_total(iy,ix) = k_y * P * Delta_s / H;
        end
    end
end

% 绘制热图
figure;
imagesc(x, y, W_total); % 注意交换 x 和 y
colorbar;
title('固定脚踩区域的磨损分布');
xlabel('台阶长度方向 (m)'); % 对应脚宽方向
ylabel('台阶宽度方向 (m)'); % 对应脚长方向
axis equal;
set(gca, 'YDir', 'normal'); % 调整坐标轴方向

%%
%最后用这个，不同上下行占比
clear,clc;
% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3;   % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.15; % 台阶高度 (m)

% 脚踩区域尺寸
L_foot = 0.25; % 脚的长度 (m)
W_foot = 0.12; % 脚的宽度 (m)

% 材料硬度
H = 2e9; % 硬度 (Pa)，如大理石

% 摩擦系数范围
k_min = 0.2; % 脚尖处摩擦系数
k_max = 0.8; % 脚跟处摩擦系数

% 压力和滑动距离
P = 100; % 压力 (N/m^2)
Delta_s = 0.01; % 滑动距离 (m)

% 非线性控制参数
n = 2; % 非线性程度控制

% 上下行比例
alpha_values = [1, 0, 0.7]; % 上行比例数组

% 网格划分
Nx = 50; % x方向网格数
Ny = 50; % y方向网格数
y = linspace(-Ly/2, Ly/2, Ny); % 台阶宽度方向 -> 脚的长度方向
x = linspace(-Lx/2, Lx/2, Nx); % 台阶长度方向 -> 脚的宽度方向
[X, Y] = meshgrid(x, y);

% 定义脚踩区域中心
foot_center_x = 0; % 脚踩区域中心的x坐标
foot_center_y = 0; % 脚踩区域中心的y坐标

% 确定脚踩区域范围
foot_x_min = foot_center_x - W_foot / 2;
foot_x_max = foot_center_x + W_foot / 2;
foot_y_min = foot_center_y - L_foot / 2;
foot_y_max = foot_center_y + L_foot / 2;

% 绘制热图
figure;
for idx = 1:length(alpha_values)
    alpha = alpha_values(idx); % 当前上行比例
    W_total = zeros(Ny, Nx); % 初始化磨损分布矩阵

    % 计算脚踩区域内的摩擦系数和磨损分布
    for ix = 1:Nx
        for iy = 1:Ny
            % 判断当前位置是否在脚踩区域内
            if x(ix) >= foot_x_min && x(ix) <= foot_x_max && ...
               y(iy) >= foot_y_min && y(iy) <= foot_y_max
                % 非线性摩擦系数分布（基于y方向）
                f_y = (y(iy) - foot_y_min) / (foot_y_max - foot_y_min); % 归一化位置
                
                % 上行与下行的摩擦系数
                k_up = k_min + (k_max - k_min) * f_y^n; % 上行模式
                k_down = k_max + (k_min - k_max) * f_y^n; % 下行模式
                
                % 综合摩擦系数
                k_total = alpha * k_up + (1 - alpha) * k_down;
                
                % 计算磨损量
                W_total(iy, ix) = k_total * P * Delta_s / H;
            end
        end
    end

    % 绘制热图
    subplot(1, length(alpha_values), idx);
    imagesc(x, y, W_total);
    colorbar;
    title(['上行比例: ', num2str(alpha)]);
    xlabel('台阶长度方向 (m)'); % 对应脚宽方向
    ylabel('台阶宽度方向 (m)'); % 对应脚长方向
    axis equal;
    set(gca, 'YDir', 'normal'); % 调整坐标轴方向
end








