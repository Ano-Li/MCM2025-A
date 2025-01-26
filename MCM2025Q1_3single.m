clear, clc;
% 台阶尺寸
Lx = 1; % 台阶宽度 (m) -> 对应脚的长度（y方向）
Ly = 0.3; % 台阶长度 (m) -> 对应脚的宽度（x方向）
Lz = 0.5; % 台阶初始高度 (m)

% 台阶使用年限和踩踏频率
T = 300 * 365; % 台阶年限为300年
f = 10; % 每天的踩踏频率为20次

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
alpha_values = [1, 0, 0.7, 0.5]; % 上行比例数组

% 网格划分
Nx = 100; % x方向网格数
Ny = 100; % y方向网格数
y = linspace(-Ly/2, Ly/2, Ny); % 台阶宽度方向 -> 脚的长度方向
x = linspace(-Lx/2, Lx/2, Nx); % 台阶长度方向 -> 脚的宽度方向
[X, Y] = meshgrid(x, y);

% 初始化台阶高度矩阵
H_step = Lz * ones(Ny, Nx); % 台阶初始高度

% 单高斯分布参数（x方向）
mu_x = 0;        % 高斯分布中心位置（脚踩区域中心）
sigma_x = 0.2;   % 高斯分布标准差，控制脚踩位置的分散程度

% \(y\) 方向随机扰动参数
mu_y = -Ly / 2 + L_foot; % 脚掌中心在 \(y\) 方向的初始位置
sigma_y = 0.05;          % \(y\) 方向的随机扰动幅度

% 绘制热图
figure('Position', [100, 100, 1200, 800]);
for idx = 1:length(alpha_values)
    alpha = alpha_values(idx); % 当前上行比例
    H_step = Lz * ones(Ny, Nx); % 重置台阶初始高度

    % 计算多次磨损叠加效果
    W_total = zeros(Ny, Nx); % 初始化磨损分布矩阵
    num_steps = 2000; % 模拟多次脚踩叠加
    for step = 1:num_steps
        % 当前脚踩位置的随机扰动
        foot_center_x = normrnd(mu_x, sigma_x); % \(x\) 方向服从单高斯分布
        foot_center_y = normrnd(mu_y, sigma_y); % \(y\) 方向加入随机性

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
                        wear = k_total * P * Delta_s * f * T / H; % 每次脚踩累积磨损量

                        % 更新磨损分布
                        W_total(iy, ix) = W_total(iy, ix) + wear;

                        % 更新台阶高度（减去磨损量）
                        H_step(iy, ix) = H_step(iy, ix) - wear;
                        H_step(iy, ix) = max(H_step(iy, ix), 0); % 确保高度非负
                    end
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
    zlim([0,0.5]);
    view(135,30); % 设置为3D视图
%     shading interp; % 平滑颜色
end
