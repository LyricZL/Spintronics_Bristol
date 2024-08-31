% CoFeB/MgO/Heavy Metal Structure Simulations

% 参数设置
gamma = 1.76e11;  % 磁旋比 (rad/s/T)
tspan = [0 1e-9];  % 时间范围

% CoFeB/MgO/重金属结构参数
alpha_CoFeB = 0.01;  % CoFeB的阻尼系数
Ms_CoFeB = 1.4e6;  % CoFeB的饱和磁化强度 (A/m)

% 初始条件
M0_CoFeB = [Ms_CoFeB; 0; 0];

% 自旋波生成模拟
H_ext_spin_wave = [0; 0; 1e3];  % 原外部磁场
[t_spin_wave, M_spin_wave] = ode45(@(t, M) lLG_equation(t, M, gamma, alpha_CoFeB, Ms_CoFeB, H_ext_spin_wave), tspan, M0_CoFeB);

% 信号放大 (STT) - 增强自旋极化电流
J_stt = 5e11;  % 增强的自旋极化电流密度 (A/m^2)
[t_stt, M_stt] = ode45(@(t, M) stt_equation_custom(t, M, gamma, alpha_CoFeB, Ms_CoFeB, H_ext_spin_wave, J_stt), tspan, M0_CoFeB);

% 信号探测 (Spin Pumping) - 改变层厚度或外部磁场
t_spin_pump = 2e-9;  % 改变的层厚度 (m)
H_ext_spin_pump = [0; 1e3; 1e3];  % 改变外部磁场方向
[t_spin_pump, M_spin_pump] = ode45(@(t, M) spin_pumping_equation(t, M, gamma, alpha_CoFeB, Ms_CoFeB, H_ext_spin_pump), tspan, M0_CoFeB);

% 绘图 - 自旋波生成
figure;
subplot(3,1,1);
plot(t_spin_wave, M_spin_wave(:,1), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_x (A/m)');
title('Magnetization Dynamics M_x (Spin Wave Generation)');

subplot(3,1,2);
plot(t_spin_wave, M_spin_wave(:,2), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_y (A/m)');
title('Magnetization Dynamics M_y (Spin Wave Generation)');

subplot(3,1,3);
plot(t_spin_wave, M_spin_wave(:,3), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_z (A/m)');
title('Magnetization Dynamics M_z (Spin Wave Generation)');

% 绘图 - 信号放大
figure;
subplot(3,1,1);
plot(t_stt, M_stt(:,1), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_x (A/m)');
title('Magnetization Dynamics M_x (STT Amplification)');

subplot(3,1,2);
plot(t_stt, M_stt(:,2), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_y (A/m)');
title('Magnetization Dynamics M_y (STT Amplification)');

subplot(3,1,3);
plot(t_stt, M_stt(:,3), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_z (A/m)');
title('Magnetization Dynamics M_z (STT Amplification)');

% 绘图 - 信号探测
figure;
subplot(3,1,1);
plot(t_spin_pump, M_spin_pump(:,1), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_x (A/m)');
title('Magnetization Dynamics M_x (Spin Pumping)');

subplot(3,1,2);
plot(t_spin_pump, M_spin_pump(:,2), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_y (A/m)');
title('Magnetization Dynamics M_y (Spin Pumping)');

subplot(3,1,3);
plot(t_spin_pump, M_spin_pump(:,3), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_z (A/m)');
title('Magnetization Dynamics M_z (Spin Pumping)');

% LLG方程的定义
function dMdt = lLG_equation(~, M, gamma, alpha, Ms, H_ext)
    H_eff = H_ext;  % 假设有效磁场为常数
    dMdt = -gamma * cross(M, H_eff) + (alpha / Ms) * cross(M, cross(M, H_eff));
end

% 自定义STT方程，加入增强自旋极化电流
function dMdt = stt_equation_custom(~, M, gamma, alpha, Ms, H_ext, J)
    H_eff = H_ext;  
    P = 0.7;  
    e = 1.6e-19;  
    hbar = 1.054e-34;  
    t = 1e-9;  
    
    tau_stt = (hbar * J * P / (2 * e * Ms * t)) * cross(M, [0; 0; 1]);
    dMdt = -gamma * cross(M, H_eff) + (alpha / Ms) * cross(M, cross(M, H_eff)) + tau_stt;
end

% 自旋泵浦方程的定义
function dMdt = spin_pumping_equation(~, M, gamma, alpha, Ms, H_ext)
    H_eff = H_ext;  % 假设有效磁场为常数
    dMdt = -gamma * cross(M, H_eff) + (alpha / Ms) * cross(M, cross(M, H_eff));
end
