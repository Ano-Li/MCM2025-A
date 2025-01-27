clear, clc;
for idx = 1:6
% file_name = ['old-stairs-1_stairs-0',num2str(idx),'.mat'];
file_name = ['stairs-of-the-17th-century_stairs-0', num2str(idx), '.mat'];
load(file_name); 
Wear = Wear - min(Wear, [], 'all');
% 矩阵维度
[Ny, Nx] = size(Wear);
x = linspace(-1, 1, Nx); % x方向坐标
y = linspace(-0.15, 0.15, Ny); % y方向坐标

% 对x方向的数据进行高斯拟合
wear_x = sum(Wear, 1); % 按列求和（x方向累积磨损量）

% 高斯拟合函数
gauss_eqn = 'a * exp(-((x-b)^2)/(2*c^2))';
x_fit = fit(x', wear_x', gauss_eqn, 'StartPoint', [max(wear_x), 0, 0.1]);

% 提取拟合参数
mu_x = x_fit.b; % 高斯分布中心（均值）
sigma_x = x_fit.c; % 高斯分布标准差

% 对y方向的数据进行高斯拟合
wear_y = sum(Wear, 2); % 按行求和（y方向累积磨损量）

% 高斯拟合函数
y_fit = fit(y', wear_y, gauss_eqn, 'StartPoint', [max(wear_y), 0, 0.01]);

% 提取拟合参数
mu_y = y_fit.b; % 高斯分布中心（均值）
sigma_y = y_fit.c; % 高斯分布标准差

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
plot(x, wear_x, 'b-', 'LineWidth', 1.5); hold on;
plot(x, x_fit(x), 'r--', 'LineWidth', 1.5);
title('Gaussian Fit for X Direction');
xlabel('X Coordinate');
ylabel('Accumulated Wear');
legend('Original Data', 'Gaussian Fit');
grid on;

% Y方向高斯拟合
subplot(2, 1, 2);
plot(y, wear_y, 'b-', 'LineWidth', 1.5); hold on;
plot(y, y_fit(y), 'r--', 'LineWidth', 1.5);
title('Gaussian Fit for Y Direction');
xlabel('Y Coordinate');
ylabel('Accumulated Wear');
legend('Original Data', 'Gaussian Fit');
grid on;

saveas(gcf, [file_name(1:end-4),'.fig']); % 保存为 MATLAB 可编辑格式
end
