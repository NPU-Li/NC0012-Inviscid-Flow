clear all
clc

% 物理域绘制
% 加密左侧和上侧外界区域
x = -30:0.1:50; % 减小步长以整体加密
for i = 1:length(x)
    if i <= 301 % 加密半圆部分
        y(i) = (900 - x(i)^2)^0.5; % 半圆区域
    elseif i <= length(x)
        y(i) = 30; % 矩形上边界
    end
end
x0 = x; 
y0 = y;

% 翼型生成 (翼型单侧一共 320 个点)
x22 = 0:0.003125:1;  
y22 = 0.1781*x22.^0.5 - 0.0756*x22 - 0.2122*x22.^2 + 0.1705*x22.^3 - 0.0609*x22.^4;

% 翼型数据储存
x_wall = [x22];
y_wall = [y22];
y_wall(:,end) = 0; % 尾缘点置为 0
l_wall = length(y_wall);

% 下边界网格生成
x_low = 3.33:0.5825:50; % 增加分布点数
y_low = zeros(size(x_low));
l_low = length(y_low);

% 左侧边界数据储存
x_left = [x0(1:l_wall)];
y_left = [y0(1:l_wall)];
l_left = length(y_left);

% 上侧矩形部分边界
x_up = [x0(l_wall+1:l_wall+l_low-1)];
y_up = [y0(l_wall+1:l_wall+l_low-1)];
l_up = length(y_up);

% 网格插值层数
n = 99; % 增加层数使整体网格更密集

% 储存位置的大数组
X = zeros(n+1, l_wall + l_low); 
Y = zeros(n+1, l_wall + l_low);

% 给 X 和 Y 的第一层赋值
X(1,1:l_wall) = x_wall; 
Y(1,1:l_wall) = y_wall;
X(1,l_wall+1:l_wall+l_low) = x_low; 
Y(1,l_wall+1:l_wall+l_low) = y_low;
X(n+1,1:l_wall) = x_left; 
Y(n+1,1:l_wall) = y_left;
X(n+1,l_wall+1:l_wall+l_up) = x_up; 
Y(n+1,l_wall+1:l_wall+l_up) = y_up;

% 插值计算网格
for j = 1:l_wall + l_up
    a1 = Y(n+1,j);
    a2 = X(n+1,j);
    a = (Y(n+1,j) - Y(1,j)) / (X(n+1,j) - X(1,j)); % 斜率
    b = Y(1,j);
    h = (X(n+1,j) - X(1,j)) / n; % 插值步长
    for i = 1:n
        % 增加尾部的加密权重
        weight = 1.0; % 默认权重
        if X(1,j) > 0.8 % 对尾部区域增加权重
            weight = 1.5;
        end
        xi = X(1,j) + i*h / weight; 
        yi = a*(xi - X(1,j)) + b;
        X(i+1,j) = xi; 
        Y(i+1,j) = yi;
    end
end

% 删除重复列
X(:,l_wall+l_up) = 50 * ones(n+1,1);
X(:,end) = []; 
Y(:,end) = [];

% 物理域所有位置的坐标矩阵
X0 = [fliplr(X) X]; 
Y0 = [fliplr(Y) -Y];
Z0 = ones(size(X0));

% 网格可视化
figure
mesh(X0, Y0, Z0); 
view(2)

% 获取网格尺寸
[M, N] = size(X0);

% 初始化
X2 = X0; 
X2(:,N/2+1) = []; 
Y2 = Y0; 
Y2(:,N/2+1) = [];
X3 = X2; 
Y3 = Y2;

[M, N] = size(X2);

% 设置迭代误差
e0 = 0.001;

% 网格优化迭代
for k = 1:3000
    for i = 2:M-1
        for j = 2:N-1
            a = ((X2(i,j+1) - X2(i,j-1))^2 + (Y2(i,j+1) - Y2(i,j-1))^2) / 4;
            b = -(((X2(i+1,j) - X2(i-1,j)) * (X2(i,j+1) - X2(i,j-1))) / 4 + ...
                 ((Y2(i+1,j) - Y2(i-1,j)) * (Y2(i,j+1) - Y2(i,j-1))) / 4);
            g = ((X2(i+1,j) - X2(i-1,j))^2 + (Y2(i+1,j) - Y2(i-1,j))^2) / 4;
            X2(i,j) = (a*(X2(i+1,j) + X2(i-1,j)) + 0.5*b*(X2(i-1,j+1) + X2(i+1,j-1) - X2(i+1,j+1) - X2(i-1,j-1)) + g*(X2(i,j+1) + X2(i,j-1))) ...
                     / (2*a + 2*g);
            Y2(i,j) = (a*(Y2(i+1,j) + Y2(i-1,j)) + 0.5*b*(Y2(i-1,j+1) + Y2(i+1,j-1) - Y2(i+1,j+1) - Y2(i-1,j-1)) + g*(Y2(i,j+1) + Y2(i,j-1))) ...
                     / (2*a + 2*g);
        end
    end
    % 对称边界更新
    X2(1,N-l_up+1:N) = fliplr(X2(1,1:l_up));
    EX = abs(X3 - X2); 
    EY = abs(Y3 - Y2);
    e = max(max(EX), max(EY));
    if e < e0
        break
    else
        X3 = X2; 
        Y3 = Y2;
    end
end

% 优化后网格可视化
figure
Z2 = ones(size(X2));
mesh(X2, Y2, Z2);
view(2)
