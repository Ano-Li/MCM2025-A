% 第一组数据
data1 = [1374.13762498490 1274.04086311502 991.064043877598 691.212233751347 801.875563378475];

% 第二组数据
data2 = [410.614394569643 315.558734826235 314.827392459691 569.778473305859 365.162609506248];

figure('Color', 'white');  % 新建图窗，背景色为白色

% 使用 boxplot 直接传入矩阵，每列作为一组
boxplot([data1', data2'], 'Labels', {'The First Stairwell','The Second Stairwell'});
hold on;

% 绘制散点（不加随机抖动）
scatter(ones(size(data1)), data1, 70, 'filled', 'MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(2*ones(size(data2)), data2, 70, 'filled', 'MarkerFaceColor','k','MarkerEdgeColor','k');

% 或者如果想让散点稍微分散，避免完全重叠，可以添加一点随机抖动：
% scatter(1 + 0.1*(rand(size(data1))-0.5), data1, 70, 'filled', 'MarkerFaceColor','k','MarkerEdgeColor','k');
% scatter(2 + 0.1*(rand(size(data2))-0.5), data2, 70, 'filled', 'MarkerFaceColor','k','MarkerEdgeColor','k');

hold off;

% 设置坐标轴标题
xlabel('Dataset Name', 'FontSize', 12);
ylabel('Age (Year)', 'FontSize', 12);
title('Age Prediction on Two Dataset', 'FontSize', 14);

% 可选：调整坐标轴刻度、图表范围等
set(gca, 'FontSize', 12);  % 让刻度字体更大一些
