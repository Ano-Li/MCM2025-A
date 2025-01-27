clear,clc

% 参数定义
k_vals = 0.05;  % 磨损系数
H = 5e10;       % 硬度 [Pa]
Fn = 686;       % 力 [N]
s0 = 0.01;      % 每次脚步的滑动距离 [m]
f0 = 2000 * 365;      % 每年的脚步次数
r = 0.02;       % 增长率参数
M = 0.5;         %容纳量

V_measured = [0.00344070319919769, 0.00319007091714000, 0.00248152525640727, ...
    0.00173072545425254, 0.00146184367861409, 0.00200781611432773, ...
    0.00102799793499907, 0.000789410771073454, 0.000787570697464598, ...
    0.00142666029777259, 0.000914022931391297];   % 磨损体积 [m^3]

% 结果存储
T_values = zeros(size(V_measured)); % 存储每个体积对应的T值

% 遍历每个磨损体积
for i = 1:length(V_measured)
    V_current = V_measured(i);
    
    % Archard's law => 总滑动距离 s_total
    s_total = (V_current * H) / (k_vals * Fn);
    
    % 定义方程 s_total = s0 * T * f0 / (1 + exp(-r * T))
    wear_equation = @(T) (M * s0 * T * f0 / (1 + (1-M) * exp(-r * T))) - s_total;
    
    % 使用数值方法求解 T 值
    T_initial_guess = 100; % 初始猜测值
    T_values(i) = fzero(wear_equation, T_initial_guess);
end

% 输出结果
disp('T values corresponding to each measured wear volume:');
disp(T_values);


