import numpy as np
import matplotlib.pyplot as plt
# 计算频率范围
k = np.linspace(0.001, 100000, 10000000)  # 波数范围，单位为 rad/m
def Sk(U_10):
    # 物理常数
    g = 9.8  # 重力加速度 m/s^2
    km = 370  # 定义波数 km，单位为 rad/m

    # 定义参数和变量  # 10米高的速度 m/s
    omega = 0.84  # 对于充分发展的海域
    k_p = g / (U_10 / omega)**2  # 波数 kp，单位为 rad/m
    global k
    
    # 计算参数
    u = U_10 * np.sqrt(0.001 * (0.81 + 0.65 * U_10))  # 海浪风速
    c_m = 0.23  # cm 参数
    c = np.sqrt(g * (1 + k**2 / km**2) / k)
    a_m = 1.4 * np.exp(-2) * u / c_m  # am 参数

    # 计算Fm
    F_m = np.exp(-0.25 * ((k / k_p) - 1)**2)

    # 计算BH
    B_H = 0.5*a_m * c_m * F_m / c

    # 计算BL
    B_L = 0.003 * np.sqrt(omega * k / k_p) * np.exp(-omega * (np.sqrt(k / k_p) - 1) / np.sqrt(10))

    # 计算L
    L = np.exp(-1.25 * (k_p / k)**2)

    # 计算Gamma和gamma
    Gamma = np.where((0.83 < omega) & (omega < 1), 1.7, 1.7 + 6 * np.log10(omega))
    gamma = np.exp(-(np.sqrt(k) - np.sqrt(k_p))**2 / (2 * (0.08 * (1 + 4 / omega**3))**2 * k_p))

    # 计算SE
    SE = k**-3 * (B_L + B_H) * L * Gamma**gamma
    return(SE)
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
plt.ylabel('lg Spectral Density (m^3/rad)')
plt.title('Elfouhaily Wave Spectrum')
plt.grid(True)
plt.legend()
plt.xlim([-3, 5])
plt.ylim([-20, 5])
plt.show()
