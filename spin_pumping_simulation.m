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
