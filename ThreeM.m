clear;
clc;

% 物理参数
gamma = 1.76e11; % 磁旋比 (rad/Ts)
alpha = 0.1;     % 增大阻尼系数
Ms = 1e5;        % 磁化强度 (A/m)
H0 = [0; 0; 1];  % 初始外加磁场 (T)
tspan = [0 5e-9]; % 延长仿真时间范围 (s)
dt = 1e-11;      % 增大时间步长 (s)
m0 = [1; 0; 0];  % 初始磁化方向 (单位向量)

% 定义时间依赖的磁场（随时间变化）
H = @(t) H0 + [0.1*sin(2*pi*1e9*t); 0; 0.1*cos(2*pi*1e9*t)];

% Runge-Kutta法求解LLG方程
[t_RK, m_RK] = RungeKuttaLLG(m0, H, alpha, gamma, Ms, tspan, dt);

% Adams-Bashforth法求解LLG方程
[t_AB, m_AB] = AdamsBashforthLLG(m0, H, alpha, gamma, Ms, tspan, dt);

% 自适应Runge-Kutta-Fehlberg法（RK45）求解LLG方程
options = odeset('RelTol',1e-6,'AbsTol',1e-8); % 设置求解精度
[t_RK45, m_RK45] = ode45(@(t,m) LLG(m, H(t), alpha, gamma, Ms), tspan, m0, options);

% 绘制结果：所有分量的变化
figure;

% m_x 分量
subplot(3, 1, 1);
plot(t_RK, m_RK(1,:), 'r', 'DisplayName', 'Runge-Kutta m_x');
hold on;
plot(t_AB, m_AB(1,:), 'b--', 'DisplayName', 'Adams-Bashforth m_x');
plot(t_RK45, m_RK45(:, 1), 'g-.', 'DisplayName', 'RK45 m_x');
xlabel('Time (s)');
ylabel('m_x');
legend;
title('m_x(t) using Runge-Kutta, Adams-Bashforth, and RK45 methods');
grid on;

% m_y 分量
subplot(3, 1, 2);
plot(t_RK, m_RK(2,:), 'r', 'DisplayName', 'Runge-Kutta m_y');
hold on;
plot(t_AB, m_AB(2,:), 'b--', 'DisplayName', 'Adams-Bashforth m_y');
plot(t_RK45, m_RK45(:, 2), 'g-.', 'DisplayName', 'RK45 m_y');
xlabel('Time (s)');
ylabel('m_y');
legend;
title('m_y(t) using Runge-Kutta, Adams-Bashforth, and RK45 methods');
grid on;

% m_z 分量
subplot(3, 1, 3);
plot(t_RK, m_RK(3,:), 'r', 'DisplayName', 'Runge-Kutta m_z');
hold on;
plot(t_AB, m_AB(3,:), 'b--', 'DisplayName', 'Adams-Bashforth m_z');
plot(t_RK45, m_RK45(:, 3), 'g-.', 'DisplayName', 'RK45 m_z');
xlabel('Time (s)');
ylabel('m_z');
legend;
title('m_z(t) using Runge-Kutta, Adams-Bashforth, and RK45 methods');
grid on;

% LLG方程的Runge-Kutta求解函数
function [t, m] = RungeKuttaLLG(m0, H, alpha, gamma, Ms, tspan, dt)
    t = tspan(1):dt:tspan(2);
    m = zeros(3, length(t));
    m(:, 1) = m0;

    for i = 1:length(t)-1
        k1 = LLG(m(:,i), H(t(i)), alpha, gamma, Ms);
        k2 = LLG(m(:,i) + 0.5*dt*k1, H(t(i)+0.5*dt), alpha, gamma, Ms);
        k3 = LLG(m(:,i) + 0.5*dt*k2, H(t(i)+0.5*dt), alpha, gamma, Ms);
        k4 = LLG(m(:,i) + dt*k3, H(t(i)+dt), alpha, gamma, Ms);
        m(:, i+1) = m(:, i) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        m(:, i+1) = m(:, i+1) / norm(m(:, i+1)); % 保持单位向量
    end
end

% LLG方程的Adams-Bashforth求解函数
function [t, m] = AdamsBashforthLLG(m0, H, alpha, gamma, Ms, tspan, dt)
    t = tspan(1):dt:tspan(2);
    m = zeros(3, length(t));
    m(:, 1) = m0;
    
    % 使用Runge-Kutta法计算前几步的值
    [~, m_temp] = RungeKuttaLLG(m(:,1), H, alpha, gamma, Ms, [t(1) t(2)], dt);
    m(:, 2) = m_temp(:,end);  % 取第二个时间点的结果
    
    [~, m_temp] = RungeKuttaLLG(m(:,2), H, alpha, gamma, Ms, [t(2) t(3)], dt);
    m(:, 3) = m_temp(:,end);  % 取第三个时间点的结果

    for i = 3:length(t)-1
        f1 = LLG(m(:, i), H(t(i)), alpha, gamma, Ms);
        f2 = LLG(m(:, i-1), H(t(i-1)), alpha, gamma, Ms);
        f3 = LLG(m(:, i-2), H(t(i-2)), alpha, gamma, Ms);
        m(:, i+1) = m(:, i) + dt * (23*f1 - 16*f2 + 5*f3) / 12;
        m(:, i+1) = m(:, i+1) / norm(m(:, i+1)); % 保持单位向量
    end
end

% LLG方程
function dm = LLG(m, H, alpha, gamma, Ms)
    Heff = H; % 时间依赖的外加磁场
    dm = -gamma/(1 + alpha^2) * (cross(m, Heff) + alpha*cross(m, cross(m, Heff)));
end
