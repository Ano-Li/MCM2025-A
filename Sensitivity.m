% 生成 T 和 Fn 的取值范围
T_values = 450:10:550; % 年份范围
Fn_values_N = linspace(450, 800, 6); % 体重（作用力）范围
Fn_values = round(Fn_values_N / 9.81); % 转换为 kg

k_vals = 0.05;  % 磨损系数
H      = 5e10;  % 硬度 [Pa]
V_measured = [0.00344070319919769, 0.00319007091714000, 0.00248152525640727, ...
    0.00173072545425254, 0.00146184367861409, 0.00200781611432773, ...
    0.00102799793499907, 0.000789410771073454, 0.000787570697464598, ...
    0.00142666029777259, 0.000914022931391297]; % 磨损体积 [m^3]
V_current = V_measured(1); % 当前选择的磨损体积

% 创建网格
[T_grid, Fn_grid] = meshgrid(T_values, Fn_values);

% 初始化每年使用频率的网格值
f_per_year_grid = zeros(size(T_grid));

% 根据每个 T 和 Fn 计算每年使用频率
for i = 1:length(T_values)
    for j = 1:length(Fn_values)
        T = T_values(i);
        Fn = Fn_values(j);

        % 根据磨损体积，计算总踩踏次数
        f_steps_current = (V_current * H) / (k_vals * Fn); 
        f_per_year_grid(j, i) = round(f_steps_current / T); % 每年使用频率并取整
    end
end

% 绘制三维柱状图
figure;

% 绘制柱状图
h = bar3(f_per_year_grid);

% 为每个柱设置颜色
for k = 1:length(h)
    z_data = h(k).ZData;
    h(k).CData = z_data; % 映射颜色
    h(k).FaceColor = 'interp'; % 设置为插值颜色
end

% 设置颜色映射
colormap('jet'); % 使用 'jet' 色图来映射值
colorbar; % 添加颜色条

% 设置图形标签和标题
xlabel('Years (T)');
ylabel('Weight Force (kg)');
zlabel('Footfalls per Year');
title(' Footfalls per Year');

% 设置横坐标具体的值
xticks(1:length(T_values));
xticklabels(arrayfun(@num2str, T_values, 'UniformOutput', false)); % 显示具体 T 值

% 设置纵坐标具体的值
yticks(1:length(Fn_values));
yticklabels(arrayfun(@num2str, Fn_values, 'UniformOutput', false)); % 显示具体 Fn 值

% % 设置颜色条范围并取整刻度
% c_min = ceil(min(f_per_year_grid(:))); % 取颜色条最小值向上取整
% c_max = floor(max(f_per_year_grid(:))); % 取颜色条最大值向下取整
% tick_values = c_min:round((c_max - c_min) / 5):c_max; % 生成均匀的整数刻度
% 
% colorbar('Ticks', tick_values, ...
%     'TickLabels', arrayfun(@num2str, tick_values, 'UniformOutput', false)); % 设置整数刻度标签

%%
% 生成 T 和 Fn 的取值范围
T_values = linspace(450, 550, 10); % 年份范围
Fn_values = linspace(600, 900, 6); % 体重（作用力）范围

k_vals = 0.05;  % 磨损系数
H      = 5e10;  % 硬度 [Pa]
V_measured = [0.00344070319919769, 0.00319007091714000, 0.00248152525640727, ...
    0.00173072545425254, 0.00146184367861409, 0.00200781611432773, ...
    0.00102799793499907, 0.000789410771073454, 0.000787570697464598, ...
    0.00142666029777259, 0.000914022931391297]; % 磨损体积 [m^3]
V_current = V_measured(1); % 当前选择的磨损体积

% 创建网格
[T_grid, Fn_grid] = meshgrid(T_values, Fn_values);

% 初始化每年使用频率的网格值
f_per_year_grid = zeros(size(T_grid));

% 根据每个 T 和 Fn 计算每年使用频率
for i = 1:length(T_values)
    for j = 1:length(Fn_values)
        T = T_values(i);
        Fn = Fn_values(j);

        % 根据磨损体积，计算总踩踏次数
        f_steps_current = (V_current * H) / (k_vals * Fn); 
        f_per_year_grid(j, i) = round(f_steps_current / T); % 每年使用频率并取整
    end
end

% 绘制三维柱状图
figure;

% 绘制柱状图
h = bar3(f_per_year_grid);

% 为每根柱子设置颜色
for k = 1:length(h)
    z_data = h(k).ZData;
    % 将每个柱子的颜色与对应的 Fn 值相关联
    % 使用 Fn 值来映射颜色
    % 获取柱子的 x 和 y 坐标，从而对应到 Fn 值
    Fn_idx = ceil(h(k).YData(1)); % 获取当前柱子对应的 Fn 值的索引
    color_map = jet(length(Fn_values)); % 使用 jet 色图
    h(k).CData = repmat(color_map(Fn_idx, :), size(z_data, 1), 1); % 根据 Fn 值映射颜色
    h(k).FaceColor = 'flat'; % 使用平面颜色
end

% 设置颜色映射
colormap('jet'); % 使用 'jet' 色图来映射值
colorbar; % 添加颜色条

% 设置图形标签和标题
xlabel('Years (T)');
ylabel('Weight Force (Fn)');
zlabel('Footfalls per Year');
title('Sensitivity Analysis: Footfalls per Year vs T and Fn');

% 设置横坐标具体的值
xticks(1:length(T_values));
xticklabels(arrayfun(@num2str, T_values, 'UniformOutput', false)); % 显示具体 T 值

% 设置纵坐标具体的值
yticks(1:length(Fn_values));
yticklabels(arrayfun(@num2str, Fn_values, 'UniformOutput', false)); % 显示具体 Fn 值

% 设置颜色条范围并取整刻度
c_min = ceil(min(f_per_year_grid(:))); % 取颜色条最小值向上取整
c_max = floor(max(f_per_year_grid(:))); % 取颜色条最大值向下取整
tick_values = c_min:round((c_max - c_min) / 5):c_max; % 生成均匀的整数刻度

colorbar('Ticks', tick_values, ...
    'TickLabels', arrayfun(@num2str, tick_values, 'UniformOutput', false)); % 设置整数刻度标签



