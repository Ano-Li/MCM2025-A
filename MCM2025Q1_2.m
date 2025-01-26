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
%%
% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3;   % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.15; % 台阶高度 (m)

% 脚踩区域尺寸（椭圆参数）
L_foot = 0.25; % 脚的长度 (m)
W_foot = 0.12; % 脚的宽度 (m)
a = L_foot / 2; % 椭圆的长轴半径
b = W_foot / 2; % 椭圆的短轴半径

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
Nx = 100; % x方向网格数
Ny = 100; % y方向网格数
y = linspace(-Ly/2, Ly/2, Ny); % 台阶宽度方向 -> 脚的长度方向
x = linspace(-Lx/2, Lx/2, Nx); % 台阶长度方向 -> 脚的宽度方向
[X, Y] = meshgrid(x, y);

% 定义脚踩区域中心
foot_center_x = 0; % 脚掌中心的x坐标
foot_center_y = 0; % 脚掌中心的y坐标

% 绘制热图
figure;
for idx = 1:length(alpha_values)
    alpha = alpha_values(idx); % 当前上行比例
    W_total = NaN(Ny, Nx); % 初始化磨损分布矩阵，默认为 NaN（背景）

    % 计算椭圆脚掌区域内的摩擦系数和磨损分布
    for ix = 1:Nx
        for iy = 1:Ny
            % 判断当前位置是否在台阶区域内
            if abs(y(iy)) <= Lx/2 && abs(x(ix)) <= Ly/2
                % 判断当前位置是否在椭圆区域内
                if ((x(ix) - foot_center_x)^2 / b^2 + (y(iy) - foot_center_y)^2 / a^2) <= 1
                    % 非线性摩擦系数分布（基于y方向）
                    f_y = (y(iy) - foot_center_y + a) / (2 * a); % 椭圆y方向归一化位置
                    
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
    end

    % Plot the heatmap
    subplot(1, length(alpha_values), idx);
    imagesc(x, y, W_total);
    colorbar;
    caxis([0 max(W_total(:))]); % Standardize color range
    title(['Upward Ratio: ', num2str(alpha)]);
    xlabel('Step Length Direction (m)'); % Corresponding to foot width direction
    ylabel('Step Width Direction (m)'); % Corresponding to foot length direction
    axis tight; % Scale the axis to fit the valid data region
    axis equal;
    set(gca, 'YDir', 'normal'); % Adjust the axis direction
end
%%
%四种比例
clear,clc;
% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3;   % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.15; % 台阶高度 (m)

% 脚踩区域尺寸（椭圆参数）
L_foot = 0.25; % 脚的长度 (m)
W_foot = 0.12; % 脚的宽度 (m)
a = L_foot / 2; % 椭圆的长轴半径
b = W_foot / 2; % 椭圆的短轴半径

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
alpha_values = [1, 0, 0.7, 0.5]; % 上行比例数组，添加了 0.5 的情况

% 网格划分
Nx = 10000; % x方向网格数
Ny = 10000; % y方向网格数
y = linspace(-Ly/2, Ly/2, Ny); % 台阶宽度方向 -> 脚的长度方向
x = linspace(-Lx/2, Lx/2, Nx); % 台阶长度方向 -> 脚的宽度方向
[X, Y] = meshgrid(x, y);

% 定义脚踩区域中心
foot_center_x = 0; % 脚掌中心的x坐标
foot_center_y = 0; % 脚掌中心的y坐标

% 绘制热图
figure('Position', [100, 100, 1200, 800]); % 增加图形的显示尺寸
for idx = 1:length(alpha_values)
    alpha = alpha_values(idx); % 当前上行比例
    W_total = NaN(Ny, Nx); % 初始化磨损分布矩阵，默认为 NaN（背景）

    % 计算椭圆脚掌区域内的摩擦系数和磨损分布
    for ix = 1:Nx
        for iy = 1:Ny
            % 判断当前位置是否在台阶区域内
            if abs(y(iy)) <= Ly/2 && abs(x(ix)) <= Lx/2
                % 判断当前位置是否在椭圆区域内
                if ((x(ix) - foot_center_x)^2 / b^2 + (y(iy) - foot_center_y)^2 / a^2) <= 1
                    % 非线性摩擦系数分布（基于y方向）
                    f_y = (y(iy) - foot_center_y + a) / (2 * a); % 椭圆y方向归一化位置
                    
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
    end

    % 绘制热图
    subplot(2, 2, idx); % 使用 2x2 的布局显示四个子图
    imagesc(x, y, W_total);
    colorbar;
    caxis([0 max(W_total(:))]); % 规范颜色范围
    title(['Upward Ratio: ', num2str(alpha)]); % 英文标题
    xlabel('Step Length Direction (m)'); % 对应脚宽方向
    ylabel('Step Width Direction (m)'); % 对应脚长方向
    axis tight; % 缩放坐标轴范围至有效数据区域
    axis equal;
    set(gca, 'YDir', 'normal'); % 调整坐标轴方向
end

%%
%以台阶高度作为初始值
% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3;   % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.15; % 台阶初始高度 (m)

T=300*365; %台阶年限为300年
f= 5000; %每天的踩踏频率为10000次；


% 脚踩区域尺寸（椭圆参数）
L_foot = 0.25; % 脚的长度 (m)
W_foot = 0.12; % 脚的宽度 (m)
a = L_foot / 2; % 椭圆的长轴半径
b = W_foot / 2; % 椭圆的短轴半径

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
alpha_values = [1, 0, 0.7, 0.5]; % 上行比例数组，添加了 0.5 的情况

% 网格划分
Nx = 100; % x方向网格数
Ny = 100; % y方向网格数
y = linspace(-Ly/2, Ly/2, Ny); % 台阶宽度方向 -> 脚的长度方向
x = linspace(-Lx/2, Lx/2, Nx); % 台阶长度方向 -> 脚的宽度方向
[X, Y] = meshgrid(x, y);

% 定义脚踩区域中心
foot_center_x = 0; % 脚掌中心的x坐标
foot_center_y = 0; % 脚掌中心的y坐标

% 初始化台阶高度矩阵
H_step = Lz * ones(Ny, Nx); % 台阶初始高度

% 绘制热图
figure('Position', [100, 100, 1200, 800]);
for idx = 1:length(alpha_values)
    alpha = alpha_values(idx); % 当前上行比例
    H_step = Lz * ones(Ny, Nx); % 重置台阶初始高度
        % 计算椭圆脚掌区域内的摩擦系数和磨损分布
        W_total = zeros(Ny, Nx); % 初始化磨损分布矩阵
        for ix = 1:Nx
            for iy = 1:Ny
                % 判断当前位置是否在台阶区域内
                if abs(y(iy)) <= Ly/2 && abs(x(ix)) <= Lx/2
                    % 判断当前位置是否在椭圆区域内
                    if ((x(ix) - foot_center_x)^2 / b^2 + (y(iy) - foot_center_y)^2 / a^2) <= 1
                        % 非线性摩擦系数分布（基于y方向）
                        f_y = (y(iy) - foot_center_y + a) / (2 * a); % 椭圆y方向归一化位置

                        % 上行与下行的摩擦系数
                        k_up = k_min + (k_max - k_min) * f_y; % 上行模式
                        k_down = k_max + (k_min - k_max) * f_y; % 下行模式

                        % 综合摩擦系数
                        k_total = alpha * k_up + (1 - alpha) * k_down;

                        % 计算当前磨损量
                        wear = k_total * P * Delta_s* f *T/ H;

                        % 更新磨损分布
                        W_total(iy, ix) = wear;

                        % 更新台阶高度（减去磨损量）
                        H_step(iy, ix) = H_step(iy, ix) - wear;
                        H_step(iy, ix) = max(H_step(iy, ix), 0); % 确保高度非负
                    end
                end
            end
    end

    % 绘制凹陷效果的3D图
    subplot(2, 2, idx);
    surf(X, Y, H_step); % 使用 surf 绘制3D图
    colorbar;
    caxis([0 Lz]); % 颜色范围与初始高度匹配
    title(['Upward Ratio: ', num2str(alpha)]); % 英文标题
    xlabel('Step Length  (m)'); % 对应脚宽方向
    ylabel('Step Width  (m)'); % 对应脚长方向
    zlabel('Step Height (m)'); % 对应台阶高度
    axis tight; % 缩放坐标轴范围至有效数据区域
    axis equal;
    view(3); % 设置为3D视图
end

%%
%改了踩踏位置靠近台阶边缘
% 清理工作区和设置参数
clear, clc;
% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3; % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.2; % 台阶初始高度 (m)

% 台阶使用年限和踩踏频率
T = 300 * 365; % 台阶年限为300年
f = 3000; % 每天的踩踏频率为5000次

% 脚踩区域尺寸（椭圆参数）
L_foot = 0.25; % 脚的长度 (m)
W_foot = 0.12; % 脚的宽度 (m)
a = L_foot / 2; % 椭圆的长轴半径
b = W_foot / 2; % 椭圆的短轴半径

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
alpha_values = [1, 0, 0.7, 0.5]; % 上行比例数组

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

% 绘制热图
figure('Position', [100, 100, 1200, 800]);
for idx = 1:length(alpha_values)
    alpha = alpha_values(idx); % 当前上行比例
    H_step = Lz * ones(Ny, Nx); % 重置台阶初始高度

    % 计算椭圆脚掌区域内的摩擦系数和磨损分布
    W_total = zeros(Ny, Nx); % 初始化磨损分布矩阵
    for ix = 1:Nx
        for iy = 1:Ny
            % 判断当前位置是否在台阶区域内
            if abs(y(iy)) <= Ly/2 && abs(x(ix)) <= Lx/2
                % 判断当前位置是否在椭圆区域内
                if ((x(ix) - foot_center_x)^2 / b^2 + (y(iy) - foot_center_y)^2 / a^2) <= 1
                    % 非线性摩擦系数分布（基于y方向）
                    f_y = (y(iy) - foot_center_y + a) / (2 * a); % 椭圆y方向归一化位置

                    % 上行与下行的摩擦系数
                    k_up = k_min + (k_max - k_min) * f_y; % 上行模式
                    k_down = k_max + (k_min - k_max) * f_y; % 下行模式

                    % 综合摩擦系数
                    k_total = alpha * k_up + (1 - alpha) * k_down;

                    % 计算当前磨损量
                    wear = k_total * P * Delta_s * f * T / H; % 每年累积磨损量

                    % 更新磨损分布
                    W_total(iy, ix) = wear;

                    % 更新台阶高度（减去磨损量）
                    H_step(iy, ix) = H_step(iy, ix) - wear;
                    H_step(iy, ix) = max(H_step(iy, ix), 0); % 确保高度非负
                end
            end
        end
    end

    % 绘制凹陷效果的3D图
    subplot(2, 2, idx);
    surf(X, Y, H_step); % 使用 surf 绘制3D图
    colorbar;
    caxis([0 Lz]); % 颜色范围与初始高度匹配
    title(['Upward Ratio: ', num2str(alpha)]); % 英文标题
    xlabel('Step Length  (m)'); % 对应脚宽方向
    ylabel('Step Width  (m)'); % 对应脚长方向
    zlabel('Step Height (m)'); % 对应台阶高度
    axis tight; % 缩放坐标轴范围至有效数据区域
    axis equal;
    view(3); % 设置为3D视图
end










