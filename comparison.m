% 参数设置
gamma = 1.76e11;  % 磁旋比 (rad/s/T)
alpha_Al = 0.005;  % 铝的阻尼系数
alpha_Au = 0.01;   % 黄金的阻尼系数
alpha_GaAs = 0.002;  % 砷化镓的阻尼系数
alpha_Ga = 0.01;     % 镓的阻尼系数
Ms_Al = 500e3;     % 铝的饱和磁化强度 (A/m)
Ms_Au = 800e3;     % 黄金的饱和磁化强度 (A/m)
Ms_GaAs = 600e3;   % 砷化镓的饱和磁化强度 (A/m)
Ms_Ga = 700e3;     % 镓的饱和磁化强度 (A/m)
H_ext = [0; 0; 1e3];  % 外部磁场 (A/m)

% 初始条件
M0_Al = [Ms_Al; 0; 0];  % 铝的初始磁化强度矢量
M0_Au = [Ms_Au; 0; 0];  % 黄金的初始磁化强度矢量
M0_GaAs = [Ms_GaAs; 0; 0];  % 砷化镓的初始磁化强度矢量
M0_Ga = [Ms_Ga; 0; 0];      % 镓的初始磁化强度矢量

% 时间设置
tspan = [0 1e-9];  % 模拟时间范围

% 使用ODE求解器进行数值模拟（铝）
[t_Al, M_Al] = ode45(@(t, M) lLG_equation(t, M, gamma, alpha_Al, Ms_Al, H_ext), tspan, M0_Al);

% 使用ODE求解器进行数值模拟（黄金）
[t_Au, M_Au] = ode45(@(t, M) lLG_equation(t, M, gamma, alpha_Au, Ms_Au, H_ext), tspan, M0_Au);

% 使用ODE求解器进行数值模拟（砷化镓）
[t_GaAs, M_GaAs] = ode45(@(t, M) lLG_equation(t, M, gamma, alpha_GaAs, Ms_GaAs, H_ext), tspan, M0_GaAs);

% 使用ODE求解器进行数值模拟（镓）
[t_Ga, M_Ga] = ode45(@(t, M) lLG_equation(t, M, gamma, alpha_Ga, Ms_Ga, H_ext), tspan, M0_Ga);

% 绘制 M_x 分量的磁化动态图
figure;
plot(t_Au, M_Au(:,1), 'r--', 'LineWidth', 1.5); % 黄金的线条
hold on;
plot(t_Ga, M_Ga(:,1), 'm:', 'LineWidth', 1.5);  % 镓的线条
plot(t_GaAs, M_GaAs(:,1), 'g-.', 'LineWidth', 1.5); % 砷化镓的线条
plot(t_Al, M_Al(:,1), 'b-', 'LineWidth', 2);    % 铝的线条
hold off;
xlabel('Time (s)');
ylabel('M_x (A/m)');
legend('Gold', 'Gallium', 'Gallium Arsenide', 'Aluminum');
title('Magnetization Dynamics M_x');

% 绘制 M_y 分量的磁化动态图
figure;
plot(t_Au, M_Au(:,2), 'r--', 'LineWidth', 1.5); % 黄金的线条
hold on;
plot(t_Ga, M_Ga(:,2), 'm:', 'LineWidth', 1.5);  % 镓的线条
plot(t_GaAs, M_GaAs(:,2), 'g-.', 'LineWidth', 1.5); % 砷化镓的线条
plot(t_Al, M_Al(:,2), 'b-', 'LineWidth', 2);    % 铝的线条
hold off;
xlabel('Time (s)');
ylabel('M_y (A/m)');
legend('Gold', 'Gallium', 'Gallium Arsenide', 'Aluminum');
title('Magnetization Dynamics M_y');

% 绘制 M_z 分量的磁化动态图
figure;
plot(t_Au, M_Au(:,3), 'r--', 'LineWidth', 1.5); % 黄金的线条
hold on;
plot(t_Ga, M_Ga(:,3), 'm:', 'LineWidth', 1.5);  % 镓的线条
plot(t_GaAs, M_GaAs(:,3), 'g-.', 'LineWidth', 1.5); % 砷化镓的线条
plot(t_Al, M_Al(:,3), 'b-', 'LineWidth', 2);    % 铝的线条
hold off;
xlabel('Time (s)');
ylabel('M_z (A/m)');
legend('Gold', 'Gallium', 'Gallium Arsenide', 'Aluminum');
title('Magnetization Dynamics M_z');

% 定义Landau-Lifshitz-Gilbert方程
function dMdt = lLG_equation(~, M, gamma, alpha, Ms, H_ext)
    H_eff = H_ext;  % 假设有效磁场为常数
    dMdt = -gamma * cross(M, H_eff) + (alpha / Ms) * cross(M, cross(M, H_eff));
end
