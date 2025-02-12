figure;
imagesc(Wear);      % 将矩阵数值映射为颜色
colormap('jet');        % 设置颜色映射方案，可换成 'parula', 'hot' 等
colorbar;               % 显示颜色刻度条
axis equal tight;       % 坐标轴设置：等比例、紧凑布局
title('Wear Data Heatmap');
xlabel('X-axis (column index)');
ylabel('Y-axis (row index)');


%% 
% 假设 wearData 仍为 700×2700
Wear = Wear - min(Wear, [], 'all');
[Ny, Nx] = size(Wear);
edge_cut_x = round(0.02 * Nx); % 切除x方向5%的边缘
edge_cut_y = round(0.05 * Ny); % 切除y方向5%的边缘
cleaned_matrix = Wear(edge_cut_y+1:end-edge_cut_y*0.2, edge_cut_x*4+1:end-edge_cut_x);
cleaned_matrix(cleaned_matrix < 1.2) = NaN;
figure;

surf(wear_1, 'EdgeColor', 'none');  % 绘制表面图，隐藏网格线
daspect([1 2 0.005]); 
colormap('jet');
colorbar;
title('3D Surface of Wear Data');
xlabel('X-axis (column index)');
ylabel('Y-axis (row index)');
zlabel('Wear Value');
axis off
grid off

wear_1 = cleaned_matrix;


% 调整视角
view(45, 30);  % 可根据需要改变俯仰角度

