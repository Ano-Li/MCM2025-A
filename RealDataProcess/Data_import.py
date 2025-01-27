# %% [markdown]
# # txt数据读入

# %%
import numpy as np
import open3d as o3d
import os

vertices_path = "/Users/ano/Library/CloudStorage/OneDrive-个人/文件/事件/硕士研究生文档/数学建模/MCM-2025/真实数据/old-stairs-1/stairs-01/Vertices.txt"
dirpath = os.path.dirname(vertices_path)
parent_dir = os.path.basename(dirpath)
grandparent_dir = os.path.basename(os.path.dirname(dirpath))

mat_filename = f"{grandparent_dir}_{parent_dir}.mat"
mat_path = os.path.join(dirpath, mat_filename)

# 1. 读入顶点数据 (Nx3)
#    假设 txt 文件是空格或 tab 分隔
vertices = np.loadtxt(vertices_path)  # 每行: x y z

# 2. 读入三角面数据 (Mx3)
#    如果来自 Blender, 通常是 1-based 索引，需要减 1 转为 0-based
faces_path = os.path.join(dirpath, "Faces.txt")
faces = np.loadtxt(faces_path, dtype=np.int32)
# faces -= 1  # 将 1-based 转为 0-based

# 检查下数据形状
print("Vertices shape:", vertices.shape)  # 期望 (N, 3)
print("Faces shape:", faces.shape)        # 期望 (M, 3)


# %% [markdown]
# # Open3D中构造三角网格并可视化

# %%
# 3. 创建 TriangleMesh 对象
mesh = o3d.geometry.TriangleMesh()
mesh.vertices = o3d.utility.Vector3dVector(vertices)
mesh.triangles = o3d.utility.Vector3iVector(faces)

# 4. 如果需要法线，计算一下顶点/面法线
mesh.compute_vertex_normals()

# 5. 可视化检查
print(mesh)  # 打印网格信息 (顶点数、面数等)
# o3d.visualization.draw_geometries([mesh])

# (B) Loop 细分（增大面数、保留大体形状）
mesh_subdiv = mesh.subdivide_loop(number_of_iterations=1)
mesh_subdiv.compute_vertex_normals()

# o3d.visualization.draw_geometries([mesh_subdiv])

# %% [markdown]
# # 网格采样点云

# %%
pcd = mesh.sample_points_poisson_disk(number_of_points=50000)

# %% [markdown]
# # 网格滤波

# %%
cl, ind = pcd.remove_statistical_outlier(nb_neighbors=50, std_ratio=2.0)
# cl 是过滤后的点云, ind 是保留点索引
pcd_filtered = cl
# o3d.visualization.draw_geometries([pcd_filtered], window_name="Filtered Cloud")

# %% [markdown]
# # 找到滤波后的最高点Z0

# %%
pts_filtered = np.asarray(pcd_filtered.points)  # shape: (Nf, 3)
z0 = np.max(pts_filtered[:, 2])
z1 = np.min(pts_filtered[:, 2])
print("Max Z after filtering =", z0)
print("Min Z after filtering =", z1)

# %% [markdown]
# # 构造规则网络

# %%
# 6. 构造规则网格 (xq, yq)
x_vals = pts_filtered[:, 0]
y_vals = pts_filtered[:, 1]
xmin, xmax = x_vals.min(), x_vals.max()
ymin, ymax = y_vals.min(), y_vals.max()

spacing = 0.01  # 1 mm

# 6.2 根据 spacing 计算所需的采样点个数
width_x = xmax - xmin
width_y = ymax - ymin

# 确保至少有 2 个点, 以免出现 width=0 或极小导致报错
Nx = max(2, int(np.ceil(width_x / spacing)) + 1)
Ny = max(2, int(np.ceil(width_y / spacing)) + 1)

print(f"Sampling Nx={Nx}, Ny={Ny}, spacing={spacing} m")

Xq, Yq = np.meshgrid(
    np.linspace(xmin, xmax, Nx),
    np.linspace(ymin, ymax, Ny)
)

# %% [markdown]
# # 用griddata差值磨损面

# %%
from scipy.interpolate import griddata
Z_worn = griddata(
    (x_vals, y_vals),         # 已过滤点云的 (x, y)
    pts_filtered[:, 2]*10,       # 已过滤点云的 z
    (Xq, Yq),
    method='linear'
)


# %% [markdown]
# # 构建初始平面

# %%
# 8. 构建“理想初始面” Z_plane (Z = z0 + 小噪声)
noise_amplitude = 0.001  # 1 mm
random_noise = np.random.normal(loc=0.0, scale=noise_amplitude, size=Xq.shape)
Z_plane = z0 + random_noise

# %% [markdown]
# # 计算磨损矩阵 wear matrix

# %%
Wear = None
if Z_worn is not None:
    Z_worn_filled = np.nan_to_num(Z_worn, nan=z0)
    Wear = Z_plane - Z_worn_filled
else:
    print("Warning: Z_worn is None, griddata failed. Check your data coverage.")


# %%
import matplotlib.pyplot as plt
import scipy.io
if Wear is not None:
    # plt.figure(figsize=(8,6))
    # plt.imshow(
    #     Wear, origin='lower', cmap='jet',
    #     extent=[xmin, xmax, ymin, ymax]
    # )
    # plt.colorbar(label="Wear Depth (m)")
    # plt.title("Wear Matrix = Z_plane - Z_worn")
    # plt.xlabel("X")
    # plt.ylabel("Y")
    # plt.show()

    scipy.io.savemat(mat_path, {"Wear": Wear})


