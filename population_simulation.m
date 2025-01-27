% 定义参数
K = 20;       % 环境承载能力（最终稳定值）
r = 0.02;       % 增长速率
t = 0:0.1:400; % 时间，从0到100，间隔0.1

% 逻辑斯蒂增长函数
f = K ./ (1 + (K - 1) * exp(-r * t));

% 绘制增长曲线
figure;
plot(t, f, 'LineWidth', 2);
grid on;
xlabel('Time (t)', 'FontSize', 12);
ylabel('Population (f(t))', 'FontSize', 12);
title('Logistic Growth Curve', 'FontSize', 14);


%% 
K = 20;        % 环境承载能力
r = 0.02;      % 增长速率

% 定义函数句柄
f_handle = @(t) K ./ (1 + (K - 1) * exp(-r * t));

% Romberg积分
I_romberg = romberg(f_handle, 0, 400, 1e-6); % 精度1e-6
disp(['龙贝格积分法结果：', num2str(I_romberg)]);

% 辅助函数
function I = romberg(f, a, b, tol)
    % Romberg积分
    R = zeros(1, 1);
    h = b - a;
    R(1,1) = h * (f(a) + f(b)) / 2;
    k = 1;
    while true
        h = h / 2;
        x = a + h:2*h:b-h;
        R(k+1,1) = R(k,1) / 2 + h * sum(f(x));
        for j = 2:k+1
            R(k+1,j) = (4^(j-1) * R(k+1,j-1) - R(k,j-1)) / (4^(j-1) - 1);
        end
        if abs(R(k+1,k+1) - R(k,k)) < tol
            I = R(k+1,k+1);
            return;
        end
        k = k + 1;
    end
end


%% 
K = 20;       % 环境承载能力
r = 0.02;     % 增长速率
t = 0:0.1:400; % 时间从0到400，步长为0.1

% 定义逻辑斯蒂增长函数
f = K ./ (1 + (K - 1) * exp(-r * t));

% 使用trapz（梯形积分法）
I_trapz = trapz(t, f);
disp(['梯形法积分结果：', num2str(I_trapz)]);

% 使用integral（自适应积分法）
f_handle = @(t) K ./ (1 + (K - 1) * exp(-r * t));
I_integral = integral(f_handle, 0, 400);
disp(['自适应积分法结果：', num2str(I_integral)]);
