clear,clc;
% 定义给定的五个数
data = [410.614	315.558	314.827	569.778	365.162];

% 生成所有三个数的组合
combinations = nchoosek(data, 3);

% 初始化结果存储
num_combinations = size(combinations, 1);
std_results = zeros(num_combinations, 1);

% 计算每个组合的标准差
for i = 1:num_combinations
    std_results(i) = std(combinations(i, :));
end

% 将组合及其标准差存入数组，并按标准差排序
results = [combinations, std_results];
results = sortrows(results, 4);