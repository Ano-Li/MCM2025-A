%体积
clear

k_vals = 0.05;  % Range of possible wear coefficients
H      = 5e10;  % Hardness [Pa]
Fn     = 686;    % wight [N]
s0     = 0.01;  % Slip distance per footstep [m]
years_used = 500;

V_measured = [0.00344070319919769	0.00319007091714000	0.00248152525640727...
    0.00173072545425254	0.00146184367861409	0.00200781611432773...
    0.00102799793499907	0.000789410771073454	0.000787570697464598	0.00142666029777259	0.000914022931391297];   % Volume of wear [m^3]

% RESULTS STORAGE
f_steps_array = zeros(size(V_measured));

for i = 1:length(V_measured)
    V_current = V_measured(i);
    % Archardns law => s = V * H / (k * Fn)
    s_total = (V_current * H) / (k_vals * Fn);
    
    % Number of footfalls
    f_steps = s_total / s0;
    f_steps_array(i) = f_steps;
end

% Convert to usage per year (assuming continuous usage over 300 years)
f_steps_per_year = f_steps_array / years_used;

% Display some results
disp('Possible usage rates (footfalls/year) for various wear coefficients:');
disp(table(V_measured', f_steps_per_year', 'VariableNames', ...
           {'V', 'f'}));

%%
%深度
% clear
% 
% k_vals = 0.05;  % Range of possible wear coefficients
% H      = 5e10;  % Hardness [Pa]
% Pn     = 30000;    % wight [Pa]
% s0     = 0.01;  % Slip distance per footstep [m]
% years_used = 500;
% 
% W_measured = [2.75494166405037	2.61586063889118	1.75092293274812...
%     3.53323352844741	1.25984857145743	2.95858067693070];   % depth of wear [cm]
% 
% % RESULTS STORAGE
% f_steps_array = zeros(size(W_measured));
% 
% for i = 1:length(W_measured)
%     W_current = W_measured(i);
%     % Archardns law => s = W * H / (k * Pn)
%     s_total = (W_current * H) / (k_vals * Pn);
%     
%     % Number of footfalls
%     f_steps = s_total / s0;
%     f_steps_array(i) = f_steps;
% end
% 
% % Convert to usage per year (assuming continuous usage over 300 years)
% f_steps_per_year = f_steps_array / years_used;
% 
% % Display some results
% disp('Possible usage rates (footfalls/year) for various wear coefficients:');
% disp(table(W_measured', f_steps_per_year', 'VariableNames', ...
%            {'W', 'f'}));
%%

plotmatrix = [V_measured;f_steps_per_year];
% 按第一行的值排序
[sorted_row, sort_index] = sort(plotmatrix(1, :)); % 获取排序后的第一行及其索引
pm_sorted = plotmatrix(:, sort_index); % 按索引对整个矩阵的列进行重排

% 创建x轴索引
x = 1:size(pm_sorted, 2);

% 绘图
figure;
yyaxis right; % 右侧y轴用于柱状图
bar(x, pm_sorted(2, :), 0.5, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none'); % 柱状图颜色设为蓝色
ylabel('steps per year');

hold on; % 保持绘图状态
yyaxis left; % 左侧y轴用于折线图
plot(x, pm_sorted(1, :), '-o', 'LineWidth', 2, 'Color', [0.9 0.3 0.3], 'MarkerFaceColor', [1 0 0]); % 折线图颜色为红色
ylabel('Volume of wear (m^3)');
ylim([0,0.0038]);

% 设置标签、标题和图例
xlabel('different stairs');
legend('steps per year', 'Volume of wear ', 'Location', 'northwest');
% title('Combined Plot: V as Line and f as Bar');
grid on;

