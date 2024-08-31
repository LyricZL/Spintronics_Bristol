function spin_wave_simulation
    % 参数设置
    gamma = 1.76e11;  % 磁旋比 (rad/s/T)
    H_ext = [0; 0; 1e3];  % 外部磁场 (A/m)
    tspan = [0 1e-9];  % 时间范围

    % 材料参数
    materials = {'Aluminum', 'Gold', 'Gallium Arsenide', 'Gallium'};
    alpha = [0.005, 0.01, 0.002, 0.01];  % 阻尼系数
    Ms = [500e3, 800e3, 600e3, 700e3];  % 饱和磁化强度 (A/m)

    % 初始条件
    M0 = cellfun(@(x) [x; 0; 0], num2cell(Ms), 'UniformOutput', false);

    % 模拟每种材料
    for i = 1:length(materials)
        [t, M] = ode45(@(t, M) lLG_equation(t, M, gamma, alpha(i), Ms(i), H_ext), tspan, M0{i});
        plot_results(t, M, materials{i});
    end
end

function dMdt = lLG_equation(~, M, gamma, alpha, Ms, H_ext)
    H_eff = H_ext;  % 假设有效磁场为常数
    dMdt = -gamma * cross(M, H_eff) + (alpha / Ms) * cross(M, cross(M, H_eff));
end

function plot_results(t, M, material)
    figure;
    subplot(3,1,1);
    plot(t, M(:,1), 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('M_x (A/m)');
    title(['Magnetization Dynamics M_x for ', material]);

    subplot(3,1,2);
    plot(t, M(:,2), 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('M_y (A/m)');
    title(['Magnetization Dynamics M_y for ', material]);

    subplot(3,1,3);
    plot(t, M(:,3), 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('M_z (A/m)');
    title(['Magnetization Dynamics M_z for ', material]);
end

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

function spin_pumping_simulation
    % 参数设置
    gamma = 1.76e11;  % 磁旋比 (rad/s/T)
    H_ext = [0; 0; 1e3];  % 外部磁场 (A/m)
    tspan = [0 1e-9];  % 时间范围

    % YIG/Pt双层结构参数
    alpha_YIG = 0.001;  % YIG的阻尼系数
    Ms_YIG = 140e3;  % YIG的饱和磁化强度 (A/m)

    % 初始条件
    M0_YIG = [Ms_YIG; 0; 0];

    % 使用ODE求解器进行数值模拟
    [t, M] = ode45(@(t, M) spin_pumping_equation(t, M, gamma, alpha_YIG, Ms_YIG, H_ext), tspan, M0_YIG);
    
    % 绘图结果
    plot_results(t, M, 'YIG/Pt Bilayer Structure');
end

function dMdt = spin_pumping_equation(~, M, gamma, alpha, Ms, H_ext)
    H_eff = H_ext;  % 假设有效磁场为常数
    dMdt = -gamma * cross(M, H_eff) + (alpha / Ms) * cross(M, cross(M, H_eff));
end
