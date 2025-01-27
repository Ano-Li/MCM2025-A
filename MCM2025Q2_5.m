clear, clc, close all;
fit_all=[];
for idx = 1:6
% file_name = ['old-stairs-1_stairs-0',num2str(idx),'.mat'];
file_name = ['stairs-of-the-17th-century_stairs-0', num2str(idx), '.mat'];
load(file_name); 
Wear = Wear - min(Wear, [], 'all');
% 矩阵维度
[Ny, Nx] = size(Wear);
x = linspace(0, 1, Nx); % x方向坐标
y = linspace(0, 0.3, Ny); % y方向坐标

% -------- 消除边缘异常值 --------
% 设置边缘区域大小
edge_cut_x = round(0.05 * Nx); % 切除x方向5%的边缘
edge_cut_y = round(0.05 * Ny); % 切除y方向5%的边缘

% 去除矩阵边缘
cleaned_matrix = Wear(edge_cut_y+1:end-edge_cut_y, edge_cut_x+1:end-edge_cut_x);

% 更新坐标范围
cleaned_x = x(edge_cut_x+1:end-edge_cut_x);
cleaned_y = y(edge_cut_y+1:end-edge_cut_y);

% 对x方向的数据进行高斯拟合
wear_x = sum(cleaned_matrix, 1); % 按列求和（x方向累积磨损量）

% 高斯拟合函数
gauss_eqn = 'a * exp(-((x-b)^2)/(2*c^2))';
x_fit = fit(cleaned_x', wear_x', gauss_eqn, 'StartPoint', [max(wear_x), 0.5, 0.1]);

% 提取拟合参数
mu_x = x_fit.b; % 高斯分布中心（均值）
sigma_x = x_fit.c; % 高斯分布标准差

% 对y方向的数据进行高斯拟合
wear_y = sum(cleaned_matrix, 2); % 按行求和（y方向累积磨损量）

% 高斯拟合函数
y_fit = fit(cleaned_y', wear_y, gauss_eqn, 'StartPoint', [max(wear_y), 0.15, 0.05]);

% 提取拟合参数
mu_y = y_fit.b; % 高斯分布中心（均值）
sigma_y = y_fit.c; % 高斯分布标准差

fit_params = [mu_x, mu_y; sigma_x, sigma_y];
fit_all = [fit_all ; fit_params];
% 显示拟合结果
disp('Gaussian Fit Results for X Direction:');
disp(['Mu_x: ', num2str(mu_x)]);
disp(['Sigma_x: ', num2str(sigma_x)]);

disp('Gaussian Fit Results for Y Direction:');
disp(['Mu_y: ', num2str(mu_y)]);
disp(['Sigma_y: ', num2str(sigma_y)]);

% 绘制拟合结果
figure('Position', [100, 100, 1200, 600]);

% X方向高斯拟合
subplot(2, 1, 1);
plot(cleaned_x, wear_x, 'b-', 'LineWidth', 1.5); hold on;
plot(cleaned_x, x_fit(cleaned_x), 'r--', 'LineWidth', 1.5);
title('Gaussian Fit for X Direction');
xlabel('X Coordinate (m)');
ylabel('Accumulated Wear (mm^3)');
legend('Original Data', 'Gaussian Fit' , 'Location', 'northeast');
grid on;

% Y方向高斯拟合
subplot(2, 1, 2);
plot(cleaned_y, wear_y, 'b-', 'LineWidth', 1.5); hold on;
plot(cleaned_y, y_fit(cleaned_y), 'r--', 'LineWidth', 1.5);
title('Gaussian Fit for Y Direction');
xlabel('Y Coordinate (m)');
ylabel('Accumulated Wear');
legend('Original Data', 'Gaussian Fit','Location', 'northeast');
grid on;

% saveas(gcf, [file_name(1:end-4),'.fig']); % 保存为 MATLAB 可编辑格式
end
