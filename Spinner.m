clear;
clc;

% 物理参数设置
gamma = 1.76e11; % 磁旋比 (rad/Ts)
Ms_sv = 8e5;     % 金属薄膜自旋阀的磁化强度 (A/m)
Ms_mtj = 1.5e6;  % 磁隧道结的磁化强度 (A/m)
Ms_stno = 5e5;   % 自旋注入自旋扭矩纳米振荡器的磁化强度 (A/m)
alpha_sv = 0.2;  % 增加自旋阀的阻尼系数
alpha_mtj = 0.05; % 增加磁隧道结的阻尼系数
alpha_stno = 0.01; % 减小自旋扭矩纳米振荡器的阻尼系数
H0 = [0; 0; 5];  % 大幅增加外加静磁场强度 (T)
f_thz = 1e12;    % 增加太赫兹频率 (Hz)
omega_thz = 2 * pi * f_thz; % 太赫兹频段角频率 (rad/s)
tspan = [0 1e-11]; % 延长时间范围 (s)
dt = 1e-14;      % 时间步长 (s)
m0 = [0.707; 0.707; 0];  % 改变初始磁化方向 (单位向量)

% 定义时间依赖的外加磁场 (包括强烈的太赫兹频段振荡)
H = @(t) H0 + [0.3*cos(omega_thz*t); 0.3*sin(omega_thz*t); 0.1*cos(2*omega_thz*t)];

% 自旋器1: 金属薄膜自旋阀 (Spin Valve)
[t_sv, m_sv] = RungeKuttaLLG(m0, H, alpha_sv, gamma, Ms_sv, tspan, dt);

% 自旋器2: 磁隧道结 (Magnetic Tunnel Junction)
[t_mtj, m_mtj] = RungeKuttaLLG(m0, H, alpha_mtj, gamma, Ms_mtj, tspan, dt);

% 自旋器3: 自旋扭矩纳米振荡器 (Spin Torque Nano-Oscillator)
[t_stno, m_stno] = RungeKuttaLLG(m0, H, alpha_stno, gamma, Ms_stno, tspan, dt);

% 绘制结果：比较不同自旋器的所有分量变化
figure;

subplot(3, 1, 1);
plot(t_sv, m_sv(1,:), 'r', 'DisplayName', 'Spin Valve m_x');
hold on;
plot(t_mtj, m_mtj(1,:), 'b--', 'DisplayName', 'MTJ m_x');
plot(t_stno, m_stno(1,:), 'g-.', 'DisplayName', 'STNO m_x');
xlabel('Time (s)');
ylabel('m_x');
legend;
title('Comparison of m_x(t) between Spintronic Devices in Terahertz Range');
grid on;

subplot(3, 1, 2);
plot(t_sv, m_sv(2,:), 'r', 'DisplayName', 'Spin Valve m_y');
hold on;
plot(t_mtj, m_mtj(2,:), 'b--', 'DisplayName', 'MTJ m_y');
plot(t_stno, m_stno(2,:), 'g-.', 'DisplayName', 'STNO m_y');
xlabel('Time (s)');
ylabel('m_y');
legend;
title('Comparison of m_y(t) between Spintronic Devices in Terahertz Range');
grid on;

subplot(3, 1, 3);
plot(t_sv, m_sv(3,:), 'r', 'DisplayName', 'Spin Valve m_z');
hold on;
plot(t_mtj, m_mtj(3,:), 'b--', 'DisplayName', 'MTJ m_z');
plot(t_stno, m_stno(3,:), 'g-.', 'DisplayName', 'STNO m_z');
xlabel('Time (s)');
ylabel('m_z');
legend;
title('Comparison of m_z(t) between Spintronic Devices in Terahertz Range');
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

% LLG方程
function dm = LLG(m, H, alpha, gamma, Ms)
    Heff = H; % 时间依赖的外加磁场
    dm = -gamma/(1 + alpha^2) * (cross(m, Heff) + alpha*cross(m, cross(m, Heff)));
end
