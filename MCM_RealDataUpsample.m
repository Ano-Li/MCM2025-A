clear

for idx = 1:5
    % 拼接文件名并加载
    file_name = ['stairs-of-the-17th-century_stairs-0', num2str(idx), '.mat'];
    load(file_name); % 加载文件，确保文件中包含变量 Wear

    % 归一化 Wear 数据
    Wear = Wear - min(Wear, [], 'all'); % 将 Wear 的最小值调整为 0
    Wear_upsample = imresize(Wear, 4*3/2);
    Wmax = max(max(Wear));
    % 计算磨损体积
    W_measured = sum(Wear, 'all') * 1e-9; % [m^3]
    W_measured_upsample = sum(Wear_upsample, 'all') * 1e-9; % [m^3]

    % 记录磨损体积
    W_maxdepth = [W_maxdepth Wmax];
    W_all = [W_all W_measured]; % 将计算结果存入数组
    W_upsample_all = [W_upsample_all W_measured_upsample];
end
