import numpy as np
import matplotlib.pyplot as plt
# 计算频率范围

def D_function(k, phi, U_10):
    Omega = 0.84
    cp = U_10 / Omega
    km = 370  # 定义波数 km，单位为 rad/m
    g = 9.8  # 重力加速度 m/s^2
    c = np.sqrt(g * (1 + k ** 2 / km ** 2) / k)
    u = U_10 * np.sqrt(0.001 * (0.81 + 0.65 * U_10))  # 海浪风速
    c_m = 0.23  # cm 参数
    delta = np.tanh(np.log(2) / 4 + 4 * (c / cp) ** 2.5 + 0.13 * u / c_m * (c_m / c) ** 2.5)
    result = np.empty((n, m))
    for i in range(n):
        for j in range(m):
            result[i, j] = (1 + delta[ i]) * np.cos(2 * phi[ j]) / (2 * np.pi)
    return result
def Sk(k,U_10):
    # 物理常数
    g = 9.8  # 重力加速度 m/s^2
    km = 370  # 定义波数 km，单位为 rad/m

    # 定义参数和变量
    omega = 0.84  # 对于充分发展的海域
    k_p = g / (U_10 / omega)**2  # 波数 kp，单位为 rad/m
    
    # 计算参数
    u = U_10 * np.sqrt(0.001 * (0.81 + 0.65 * U_10))  # 海浪风速
    c_m = 0.23  # cm 参数
    c = np.sqrt(g * (1 + k**2 / km**2) / k)
    a_m = 1.4 * np.exp(-2) * u / c_m  # am 参数

    # 计算Fm
    F_m = np.exp(-0.25 * ((k / k_p) - 1)**2)

    # 计算BH
    B_H = 0.5 * a_m * c_m * F_m / c

    # 计算BL
    B_L = 0.003 * np.sqrt(omega * k / k_p) * np.exp(-omega * (np.sqrt(k / k_p) - 1) / np.sqrt(10))

    # 计算L
    L = np.exp(-1.25 * (k_p / k)**2)

    # 计算Gamma和gamma
    Gamma = 1.7
    gamma = np.exp(-(np.sqrt(k) - np.sqrt(k_p))**2 / (2 * (0.08 * (1 + 4 / omega**3))**2 * k_p))

    # 计算SE
    SE = k**-3 * (B_L + B_H) * L * Gamma**gamma
    return(SE)
'''
SE_1 = Sk(3)
SE_2 = Sk(5)
SE_3 = Sk(7)
SE_4 = Sk(9)

# 绘图
plt.figure(figsize=(10, 6))
plt.plot(np.log10(k), np.log10(SE_1), label=f'U_10 = {3} m/s')
plt.plot(np.log10(k), np.log10(SE_2), label=f'U_10 = {5} m/s')
plt.plot(np.log10(k), np.log10(SE_3), label=f'U_10 = {7} m/s')
plt.plot(np.log10(k), np.log10(SE_4), label=f'U_10 = {9} m/s')
plt.xlabel('lg Wave Number (rad/m)')
plt.ylabel('lg Spectral Density (m^2/rad)')
plt.title('Elfouhaily Wave Spectrum')
plt.grid(True)
plt.legend()
plt.xlim([-3, 5])
plt.ylim([-20, 5])
plt.show()'''
import numpy as np
import math
# 定义函数计算振幅 aij

# 定义函数计算海面模型
def generate_sea_surface_model(x, y, t, n, m, A, omega, k, phi, epsilon):
    z = np.zeros((200,200))
    for q in range(200):
        for l in range(200):
            for i in range(n):
                for j in range(m):
                    z[q][l] += A[i][j] *math.cos(omega[i] * t - k[i] * (x[q] * math.cos(phi[j]) + y[l] * math.sin(phi[j])) + epsilon[i][j])

    return z

# 测试代码
# 定义参数
n = 30  # 波的数量
m = 10  # 方向的数量
x = np.linspace(0, 50, 200)  # 水平方位坐标
y = np.linspace(0, 50, 200)  # 水平方位坐标
t = 0  # 时刻
omega = np.random.rand(n)  # 随机角频率
k=np.linspace(0.01,10, n)
dk=k[1]-k[0]
phi =np.linspace(0.01, 2*np.pi, m)
dphi=phi[1]-phi[0]
e_x=np.linspace(0.01, 2*np.pi, n)
e_y=np.linspace(0.01, 2*np.pi, m)
epsilon = np.random.rand(n,m)*np.pi*2 # 随机初始相位
se=Sk(k,4)

D=D_function(k,phi,4)

A = np.zeros((n,m))
for i in range(n):
    for j in range(m):
        A[i][j]=np.sqrt(math.fabs(2*se[i]*D[i][j]/k[i]*dphi*dk))
# 生成海面模型
z = generate_sea_surface_model(x, y, t, n, m, A, omega, k, phi, epsilon)
# 绘图
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

X, Y = np.meshgrid(x, y)
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, z, cmap='viridis')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
plt.title('Sea Surface Model')
plt.colorbar(surf, label='Z (m)')
plt.show()
