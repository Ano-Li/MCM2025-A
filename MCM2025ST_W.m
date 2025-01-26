W_all = []; % 初始化磨损体积数组

for idx = 1:6
    % 拼接文件名并加载
    file_name = ['old-stairs-1_stairs-0', num2str(idx), '.mat'];
    load(file_name); % 加载文件，确保文件中包含变量 Wear

    % 归一化 Wear 数据
    Wear = Wear - min(Wear, [], 'all'); % 将 Wear 的最小值调整为 0

    % 计算磨损体积
    W_measured = sum(Wear, 'all') * 4^3 * 1e-9; % [m^3]

    % 记录磨损体积
    W_all = [W_all W_measured]; % 将计算结果存入数组
end
