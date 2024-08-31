function stt_amplification_simulation
    % 参数设置
    gamma = 1.76e11;  % 磁旋比 (rad/s/T)
    H_ext = [0; 0; 1e3];  % 外部磁场 (A/m)
    tspan = [0 1e-9];  % 时间范围

    % CoFeB/MgO/重金属结构参数
    alpha_CoFeB = 0.01;  % CoFeB的阻尼系数
    Ms_CoFeB = 1.4e6;  % CoFeB的饱和磁化强度 (A/m)

    % 初始条件
    M0_CoFeB = [Ms_CoFeB; 0; 0];

    % 使用ODE求解器进行数值模拟
    [t, M] = ode45(@(t, M) stt_equation(t, M, gamma, alpha_CoFeB, Ms_CoFeB, H_ext), tspan, M0_CoFeB);
    
    % 绘图结果
    plot_results(t, M, 'CoFeB/MgO/Heavy Metal Structure');
end

function dMdt = stt_equation(~, M, gamma, alpha, Ms, H_ext)
    H_eff = H_ext;  % 假设有效磁场为常数
    J = 1e11;  % 自旋极化电流密度 (A/m^2)
    P = 0.7;  % 自旋极化率
    e = 1.6e-19;  % 电子电荷 (C)
    hbar = 1.054e-34;  % 约化普朗克常数 (J·s)
    t = 1e-9;  % 层厚度 (m)
    
    % 自旋转移力矩项
    tau_stt = (hbar * J * P / (2 * e * Ms * t)) * cross(M, [0; 0; 1]);
    dMdt = -gamma * cross(M, H_eff) + (alpha / Ms) * cross(M, cross(M, H_eff)) + tau_stt;
end
