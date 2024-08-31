function integrated_spintronic_thz_system
    % 参数设置
    gamma = 1.76e11;  % 磁旋比 (rad/s/T)
    H_ext = [0; 0; 1e3];  % 外部磁场 (A/m)
    tspan = [0 1e-9];  % 时间范围

    % 材料参数设置
    materials = {'Aluminum', 'Gold', 'Gallium Arsenide', 'Gallium'};
    alpha = [0.005, 0.01, 0.002, 0.01];  % 阻尼系数
    Ms = [500e3, 800e3, 600e3, 700e3];  % 饱和磁化强度 (A/m)

    % 初始条件
    M0 = cellfun(@(x) [x; 0; 0], num2cell(Ms), 'UniformOutput', false);

    % 信号生成与放大模拟
    for i = 1:length(materials)
        % 自旋波生成
        [t_gen, M_gen] = ode45(@(t, M) lLG_equation(t, M, gamma, alpha(i), Ms(i), H_ext), tspan, M0{i});
        
        % 信号放大
        [t_amp, M_amp] = ode45(@(t, M) stt_equation(t, M, gamma, alpha(i), Ms(i), H_ext), tspan, M_gen(end,:));
        
        % 信号探测
        [t_det, M_det] = ode45(@(t, M) spin_pumping_equation(t, M, gamma, alpha(i), Ms(i), H_ext), tspan, M_amp(end,:));

        % 绘图结果
        figure;
        subplot(3,1,1);
        plot(t_gen, M_gen(:,1), 'LineWidth', 1.5);
        hold on;
        plot(t_amp, M_amp(:,1), '--', 'LineWidth', 1.5);
        plot(t_det, M_det(:,1), ':', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('M_x (A/m)');
        title(['Magnetization Dynamics M_x for ', materials{i}]);

        subplot(3,1,2);
        plot(t_gen, M_gen(:,2), 'LineWidth', 1.5);
        hold on;
        plot(t_amp, M_amp(:,2), '--', 'LineWidth', 1.5);
        plot(t_det, M_det(:,2), ':', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('M_y (A/m)');
        title(['Magnetization Dynamics M_y for ', materials{i}]);

        subplot(3,1,3);
        plot(t_gen, M_gen(:,3), 'LineWidth', 1.5);
        hold on;
        plot(t_amp, M_amp(:,3), '--', 'LineWidth', 1.5);
        plot(t_det, M_det(:,3), ':', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('M_z (A/m)');
        title(['Magnetization Dynamics M_z for ', materials{i}]);
        legend('Generation', 'Amplification', 'Detection');
    end
end
