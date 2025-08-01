function [xk_plus, Pk_plus, n_feature] = UKF_velocity(r1, v1, P1, et_1, et_2, F, V, N, fov1, fov2, sigma_vel, meas_vel, meas_vel_index, C20, C22)

    alpha = 0.01;
    beta = 2;
    k = 0;
    n = 6;
    lambda = alpha^2 * (n + k) - n;

    omega_body = 2*pi/18972.92*[0; 0; 1];
    mass_eros = 6.687e15;

    x1 = [r1; v1];
    options = odeset('reltol', 1e-12, 'abstol', [1e-8*ones(3,1); 1e-11*ones(3,1)]);

    % Generazione sigma points
    [sp] = sigma_points(x1, P1, alpha);
    chi = zeros(n, size(sp, 2));

    gamma = [];
    
    meas_vel = [meas_vel(:, 1); meas_vel(:, 2)];
    meas_vel_index = [meas_vel_index; meas_vel_index];

    for i = 1:size(sp, 2)
        [~, propag] = ode78(@(t,x) model_rotating_dynamics(t, x, mass_eros, omega_body, C20, C22), [et_1, et_2], sp(:, i), options);
        chi(:, i) = propag(end, :)';
    end

    for i = 1:size(sp, 2)
        [v_cam, v_cam_index] = ground_velocity_measurement(sp(1:3, i)', sp(4:6, i)', chi(1:3, i)', chi(4:6,i)', F, V, N, fov1, fov2, et_1, et_2);
            
            v_cam = [v_cam(:, 1); v_cam(:, 2)];
            v_cam_index = [v_cam_index; v_cam_index];

            [C, ia, ib] = intersect(v_cam_index, meas_vel_index, 'stable');
            
            v_cam_index = v_cam_index(ia);
            v_cam = v_cam(ia);
            meas_vel_index = meas_vel_index(ib);
            meas_vel = meas_vel(ib);
 
               if i > 1
                   gamma = gamma(ib, :);
               end
            %end
        
        gamma = [gamma, v_cam];
        
    end

    if isempty(meas_vel)
        chi0 = chi(:,1);
        chii = chi(:,2:end);
        mean_UT_f = zeros(6, 1);
        P_UT_f = zeros(6,6);
        
        %Compute mean's weights
        for i = 1:size(chii, 2)
            Wi = 1/(2*(n+k));
            mean_UT_f = mean_UT_f + Wi*chii(:, i);
        end
        
        W0m = 1 - n/(alpha^2*(n+k));
        %compute mean
        xk_plus = mean_UT_f + W0m*chi0;
        
        %compute P's weights
        for i = 1:size(chii, 2)
            Wi = 1/(2*(n+k));
            P_UT_f = P_UT_f + Wi*(chii(:, i) - xk_plus)*(chii(:, i) - xk_plus)';
        end
        W0c = (2-alpha^2+beta) - n/(alpha^2*(n+k));
        
        %compute P
        Pk_plus = P_UT_f + W0c*(chi0 - xk_plus)*(chi0 - xk_plus)';

        n_feature = 0;
        return
    end

    %FIND MINUS VALUES FOR POSITIONS
    chi0 = chi(:, 1);
    chi_i = chi(:, 2:end);
    xk_minus = 0;

    for j = 1:size(chi_i, 2)
    Wi = 1/(2*(n+lambda));
    xk_minus = xk_minus + Wi*chi_i(:, j);
    end    
    W0m = (lambda)/(lambda+n);
    xk_minus = xk_minus + W0m*chi0;
    
    %FIND MINUS VALUES FOR MEASUREMENTS 
    gamma0 = gamma(:,1);
    gamma_i = gamma(:,2:end);
    yk_minus = 0;

    for j = 1:size(gamma_i, 2)
    Wi = 1/(2*(n+lambda));
    yk_minus = yk_minus + Wi*gamma_i(:, j);
    end
    W0m = (lambda)/(lambda+n);
    yk_minus = yk_minus + W0m*gamma0;
    
    %FIND COVARIANCE MATRICES VALUES
    Pk_minus = zeros(length(xk_minus));
    Peek = zeros(length(yk_minus));
    Pxyk = zeros(length(xk_minus), length(yk_minus));
    
    yk = meas_vel;
    variances = [repmat(sigma_vel, size(meas_vel, 1), 1)];
    Rk = diag(variances);
    
    for j = 1:size(gamma_i, 2)
    Wic = 1/(2*(n+lambda));
    Pk_minus = Pk_minus + Wic*(chi_i(:, j) - xk_minus)*(chi_i(:, j) - xk_minus)';
    Peek = Peek + Wic*(gamma_i(:, j) - yk_minus)*(gamma_i(:, j) - yk_minus)' ;
    Pxyk = Pxyk + Wic*(chi_i(:, j) - xk_minus)*(gamma_i(:, j) - yk_minus)';
    end

    Peek = Peek + Rk;

    W0c = (lambda)/(lambda+n) + (1-alpha^2+beta);
    Pk_minus = Pk_minus + W0c*(chi0 - xk_minus)*(chi0 - xk_minus)';
    Peek = Peek + W0c*(gamma0 - yk_minus)*(gamma0 - yk_minus)';
    Pxyk = Pxyk + W0c*(chi0 - xk_minus)*(gamma0 - yk_minus)';

    Kk = Pxyk/Peek;
    
    %COMPUTE "PLUS" VALUES
    xk_plus = xk_minus + Kk*(yk - yk_minus);
    Pk_plus = Pk_minus - Kk*Peek*Kk';

    n_feature = length(yk)/2;
end
