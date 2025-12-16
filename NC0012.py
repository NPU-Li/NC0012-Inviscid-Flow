import os
from math import *
import numpy as np
import matplotlib.pyplot as plt

'''迭代函数'''


def grid_gen(Locx, Locy, xi, yj, Δη, Δξ):
    if yj + 2 > len(Locx[xi]):
        yjz1 = 0
    else:
        yjz1 = yj + 1
    α = ((Locx[xi + 1][yj] - Locx[xi - 1][yj]) / 2 / Δη) ** 2 + ((Locy[xi + 1][yj] - Locy[xi - 1][yj]) / 2 / Δη) ** 2
    β = (Locx[xi][yjz1] - Locx[xi][yj - 1]) / 2 / Δξ * (Locx[xi + 1][yj] - Locx[xi - 1][yj]) / 2 / Δη + (
                Locy[xi][yjz1] - Locy[xi][yj - 1]) / 2 / Δξ * (Locy[xi + 1][yj] - Locy[xi - 1][yj]) / 2 / Δη
    γ = ((Locx[xi][yjz1] - Locx[xi][yj - 1]) / 2 / Δξ) ** 2 + ((Locy[xi][yjz1] - Locy[xi][yj - 1]) / 2 / Δξ) ** 2
    xij = 0.5 * ((α * (Locx[xi][yjz1] + Locx[xi][yj - 1]) + γ * (Locx[xi + 1][yj] + Locx[xi - 1][yj])) - 0.5 * β * (
                Locx[xi + 1][yjz1] + Locx[xi - 1][yj - 1] - Locx[xi - 1][yjz1] - Locx[xi + 1][yj - 1])) / (α + γ)
    yij = 0.5 * ((α * (Locy[xi][yjz1] + Locy[xi][yj - 1]) + γ * (Locy[xi + 1][yj] + Locy[xi - 1][yj])) - 0.5 * β * (
                Locy[xi + 1][yjz1] + Locy[xi - 1][yj - 1] - Locy[xi - 1][yjz1] - Locy[xi + 1][yj - 1])) / (α + γ)
    return [xij, yij]


'''翼型方程'''
fy_ul_up = lambda x: 0.6 * (-0.1015 * x ** 4 + 0.2843 * x ** 3 - 0.3576 * x ** 2 - 0.1221 * x + 0.2969 * x ** 0.5)
fy_ul_down = lambda x: -0.6 * (-0.1015 * x ** 4 + 0.2843 * x ** 3 - 0.3576 * x ** 2 - 0.1221 * x + 0.2969 * x ** 0.5)

'''网格参数'''
c = 1  # 弦长
n_boundary = 201  # 边界点数，原则上为偶数
n_lay = 101  # 网格层数，不得小于3
r_far = 50 * c  # 远场半径
Δη = 1;
Δξ = 1
'''远场边界曲线'''
fy_far_up = lambda x: (r_far ** 2 - (x - c / 2) ** 2) ** 0.5
fy_far_down = lambda x: -(r_far ** 2 - (x - c / 2) ** 2) ** 0.5
'''确定壁面及远场点坐标'''
f_xwall = lambda x: c / 2 * (1 + cos(x))  # 壁面横坐标分布
X_wall = []  # 壁面横坐标
Y_wall = []  # 壁面纵坐标
f_xfar = lambda x: c / 2 + r_far * cos(x)
X_far = []  # 远场横坐标
Y_far = []  # 远场纵坐标
for i in range(n_boundary):  # 确定壁面点坐标
    X_wall.append(f_xwall(2 * pi * i / n_boundary))
    if i < n_boundary / 2:
        Y_wall.append(fy_ul_down(X_wall[i]))
    else:
        Y_wall.append(fy_ul_up(X_wall[i]))

    X_far.append(f_xfar(2 * pi * i / n_boundary))
    if i < n_boundary / 2:
        Y_far.append(fy_far_down(X_far[i]))
    else:
        Y_far.append(fy_far_up(X_far[i]))
# 绘制翼型壁面点
# plt.plot(X_wall,Y_wall);plt.gca().set_aspect(1)
# plt.plot(X_far,Y_far);plt.gca().set_aspect(1)

'''插值，获得初始条件'''
Location_x = [X_wall];
Location_y = [Y_wall]
for i in range(n_lay - 2):
    Location_x.append([]);
    Location_y.append([]);
    for j in range(n_boundary):
        Location_x[i + 1].append((X_far[j] - X_wall[j]) / (n_lay - 1) * (i + 1) + X_wall[j])
        Location_y[i + 1].append((Y_far[j] - Y_wall[j]) / (n_lay - 1) * (i + 1) + Y_wall[j])

    # plt.plot(Location_x[i+1],Location_y[i+1]);plt.gca().set_aspect(1)
Location_x.append(X_far);
Location_y.append(Y_far)
'''求解网格坐标'''
maxesp = 10
while maxesp > 1e-3:
    Esp = []
    for i in range(1, n_lay - 1):
        for j in range(n_boundary):
            x = grid_gen(Location_x, Location_y, i, j, Δη, Δξ)[0]
            y = grid_gen(Location_x, Location_y, i, j, Δη, Δξ)[1]
            Esp.append(abs(Location_x[i][j] - x))
            Esp.append(abs(Location_y[i][j] - y))
            Location_x[i][j] = x
            Location_y[i][j] = y
    maxesp = max(Esp)
    print(maxesp)
for i in range(n_lay):
    Location_x[i].append(Location_x[i][0])
    Location_y[i].append(Location_y[i][0])

'''流场参数'''
Δη = 1;
Δξ = 1
V_far = 25  # 自由来流速度，m/s
Atk = 0  # 翼型攻角,°
Φ = []
'''构建远场边界条件'''
φ_boundary = []
for i in range(n_lay):
    Φ.append([])
    for j in range(n_boundary + 1):
        Φ[i].append(V_far * cos(Atk / 180 * pi) * Location_x[i][j] + V_far * sin(Atk / 180 * pi) * Location_y[i][j])

'''壁面边界条件计算函数'''  # 二阶精度，对壁面βγ值的求解中，对壁面点的二阶精度差分，二阶精度多项式处理


def fφ_wall(Locx, Locy, Φ_wall, Φ_wall1, Φ_wall2, ξi):
    if ξi != 0:
        βξi = (Locx[0][ξi + 1] - Locx[0][ξi - 1]) / 2 / Δξ * (
                    -3 * Locx[0][ξi] + 4 * Locx[1][ξi] - Locx[2][ξi]) / 2 / Δη + (
                          Locy[0][ξi + 1] - Locy[0][ξi - 1]) / 2 / Δξ * (
                          -3 * Locy[0][ξi] + 4 * Locy[1][ξi] - Locy[2][ξi]) / 2 / Δη
        γξi = ((Locx[0][ξi + 1] - Locx[0][ξi - 1]) / 2 / Δξ) ** 2 + ((Locy[0][ξi + 1] - Locy[0][ξi - 1]) / 2 / Δξ) ** 2
        φi0 = 1 / 3 * (4 * Φ_wall1[ξi] - Φ_wall2[ξi] - 2 * βξi / γξi * Δη / Δξ / 2 * (Φ_wall[ξi + 1] - Φ_wall[ξi - 1]))
    else:

        # φi0=(4*Φ_wall1[0]-Φ_wall2[0])/3
        βξi = (Locx[0][ξi + 1] - Locx[0][-2]) / 2 / Δξ * (-3 * Locx[0][ξi] + 4 * Locx[1][ξi] - Locx[2][ξi]) / 2 / Δη + (
                    Locy[0][ξi + 1] - Locy[0][-2]) / 2 / Δξ * (
                          -3 * Locy[0][ξi] + 4 * Locy[1][ξi] - Locy[2][ξi]) / 2 / Δη
        γξi = ((Locx[0][ξi + 1] - Locx[0][-2]) / 2 / Δξ) ** 2 + ((Locy[0][ξi + 1] - Locy[0][-2]) / 2 / Δξ) ** 2
        φi0 = 1 / 3 * (4 * Φ_wall1[ξi] - Φ_wall2[ξi] - 2 * βξi / γξi * Δη / Δξ * (Φ_wall[ξi + 1] - Φ_wall[ξi]))

    return φi0


'''初始化'''  # 令流场各点速度等于自由来流速度

'''迭代求解'''


def fΦ(Φ_list, Locx, Locy, xi, yj, Δη, Δξ):
    if yj != 0:
        α = ((Locx[xi + 1][yj] - Locx[xi - 1][yj]) / 2 / Δη) ** 2 + (
                    (Locy[xi + 1][yj] - Locy[xi - 1][yj]) / 2 / Δη) ** 2
        β = (Locx[xi][yj + 1] - Locx[xi][yj - 1]) / 2 / Δξ * (Locx[xi + 1][yj] - Locx[xi - 1][yj]) / 2 / Δη + (
                    Locy[xi][yj + 1] - Locy[xi][yj - 1]) / 2 / Δξ * (Locy[xi + 1][yj] - Locy[xi - 1][yj]) / 2 / Δη
        γ = ((Locx[xi][yj + 1] - Locx[xi][yj - 1]) / 2 / Δξ) ** 2 + (
                    (Locy[xi][yj + 1] - Locy[xi][yj - 1]) / 2 / Δξ) ** 2
        Φij = 0.5 * ((α * (Φ_list[xi][yj + 1] + Φ_list[xi][yj - 1]) + γ * (
                    Φ_list[xi + 1][yj] + Φ_list[xi - 1][yj])) - 0.5 * β * (
                                 Φ_list[xi + 1][yj + 1] + Φ_list[xi - 1][yj - 1] - Φ_list[xi - 1][yj + 1] -
                                 Φ_list[xi + 1][yj - 1])) / (α + γ)
    else:

        # Φij3=(Φ_list[xi][1]+Φ_list[xi][2])/2;Φij=Φij3

        # bc-6边界条件
        α = ((Locx[xi + 1][yj] - Locx[xi - 1][yj]) / 2 / Δη) ** 2 + (
                    (Locy[xi + 1][yj] - Locy[xi - 1][yj]) / 2 / Δη) ** 2
        β = (Locx[xi][yj + 1] - Locx[xi][-2]) / 2 / Δξ * (Locx[xi + 1][yj] - Locx[xi - 1][yj]) / 2 / Δη + (
                    Locy[xi][yj + 1] - Locy[xi][-2]) / 2 / Δξ * (Locy[xi + 1][yj] - Locy[xi - 1][yj]) / 2 / Δη
        γ = ((Locx[xi][yj + 1] - Locx[xi][-2]) / 2 / Δξ) ** 2 + ((Locy[xi][yj + 1] - Locy[xi][-2]) / 2 / Δξ) ** 2
        Φij = 0.5 * ((α * (Φ_list[xi][yj + 1] + Φ_list[xi][-2] - Γ) + γ * (
                    Φ_list[xi + 1][yj] + Φ_list[xi - 1][yj])) - 0.5 * β * (
                                 Φ_list[xi + 1][yj + 1] + Φ_list[xi - 1][-2] - Γ - Φ_list[xi - 1][yj + 1] -
                                 Φ_list[xi + 1][-2] + Γ)) / (α + γ)

    return Φij


Γ = 100
Γ1 = 0
while abs(Γ - Γ1) > 1e-5:
    # for t in range(1000):
    Γ1 = Γ
    Γ = Φ[0][-2] - Φ[0][1]
    # Γ=0
    print(Γ, Γ - Γ1)

    for j in range(n_lay):
        Φ[j][-1] == Φ[j][n_boundary]
        Φ[j][-1] = Φ[j][0] + Γ
    for j in range(2, n_lay):
        for i in range(n_boundary):
            Φ[n_lay - j][i] = fΦ(Φ, Location_x, Location_y, n_lay - j, i, Δη, Δξ)
    for i in range(n_boundary):
        Φ[0][i] = fφ_wall(Location_x, Location_y, Φ[0], Φ[1], Φ[2], i)

'''由速度势求解速度'''
U = [[0]];
V = [[0]];
Cp = [[-1]]
for i in range(1, n_boundary):
    xξ = (Location_x[0][i + 1] - Location_x[0][i - 1]) / 2 / Δξ
    yξ = (Location_y[0][i + 1] - Location_y[0][i - 1]) / 2 / Δξ
    xη = (-3 * Location_x[0][i] + 4 * Location_x[1][i] - Location_x[2][i]) / 2 / Δη
    yη = (-3 * Location_y[0][i] + 4 * Location_y[1][i] - Location_y[2][i]) / 2 / Δη
    φξ = (Φ[0][i + 1] - Φ[0][i - 1]) / 2 / Δξ
    φη = (-3 * Φ[0][i] + 4 * Φ[1][i] - Φ[2][i]) / 2 / Δη
    α = xη ** 2 + yη ** 2
    β = xξ * xη + yξ * yη
    γ = xξ ** 2 + yξ ** 2
    J = xξ * yη - xη * yξ
    U[0].append(1 / J * (φξ * yη - φη * yξ))
    V[0].append(1 / J * (φη * xξ - φξ * xη))
    Cp[0].append((U[0][i] ** 2 + V[0][i] ** 2) / V_far ** 2 - 1)
U[0].append(0)
V[0].append(0)
Cp[0].append(-1)
for j in range(1, n_lay - 1):
    U.append([]);
    V.append([]);
    Cp.append([])
    for i in range(n_boundary + 1):
        if i == n_boundary or i == 0:
            xξ = (Location_x[j][1] - Location_x[j][-2]) / 2 / Δξ
            yξ = (Location_y[j][1] - Location_y[j][-2]) / 2 / Δξ
        else:
            xξ = (Location_x[j][i + 1] - Location_x[j][i - 1]) / 2 / Δξ
            yξ = (Location_y[j][i + 1] - Location_y[j][i - 1]) / 2 / Δξ
        xη = (Location_x[j + 1][i] - Location_x[j - 1][i]) / 2 / Δη
        yη = (Location_y[j + 1][i] - Location_y[j - 1][i]) / 2 / Δη
        φη = (Φ[j + 1][i] - Φ[j - 1][i]) / 2 / Δη
        if i == 0:
            φξ = (-3 * Φ[j][i] + 4 * Φ[j][i + 1] - Φ[j][i + 2]) / 2 / Δξ
        else:
            if i == n_boundary:
                φξ = (-3 * Φ[j][i] + 4 * Φ[j][i - 1] - Φ[j][i - 2]) / (-2) / Δξ
            else:
                φξ = (Φ[j][i + 1] - Φ[j][i - 1]) / 2 / Δξ
        α = xη ** 2 + yη ** 2
        β = xξ * xη + yξ * yη
        γ = xξ ** 2 + yξ ** 2
        J = xξ * yη - xη * yξ
        U[j].append(1 / J * (φξ * yη - φη * yξ))
        V[j].append(1 / J * (φη * xξ - φξ * xη))
        Cp[j].append((U[j][i] ** 2 + V[j][i] ** 2) / V_far ** 2 - 1)

U.append([]);
V.append([]);
Cp.append([])
for i in range(n_boundary + 1):
    U[-1].append(V_far * cos(Atk / 180 * pi))
    V[-1].append(V_far * sin(Atk / 180 * pi))
    Cp[-1].append(0)

plt.scatter(Location_x[0], Cp[0])

text = open('naca0012-%s%s网格-%s攻角.dat' % (n_boundary, n_lay, Atk), 'w')
print('Title=\"NACA0012\"', file=text)
print('VARIABLES=X,Y,U,V', file=text)
print('ZONE T="field", I=%s  J=%s' % (n_boundary + 1, n_lay), file=text)
for j in range(n_lay):
    for i in range(n_boundary + 1):
        print('%s     %s     %s     %s' % (Location_x[j][i], Location_y[j][i], U[j][i], V[j][i]), file=text)
text.close()

text = open('naca0012-%s%s网格-%s攻角-cp.dat' % (n_boundary, n_lay, Atk), 'w')
print('Title=\"NACA0012\"', file=text)
print('VARIABLES=x/c,Cp', file=text)
print('ZONE T="cal"', file=text)

for i in range(n_boundary + 1):
    print('%s     %s' % (Location_x[0][i], Cp[0][i]), file=text)
text.close()