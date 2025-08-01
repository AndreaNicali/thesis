function [xk_plus, Pk_plus, n_feature] = UKF(r0, v0, P0, et_i, et_f, F, V, N, features, fov1, fov2, sigma_angular, sigma_range, measurements_angular, measurement_range, C20, C22)

    alpha = 1;
    beta = 2;
    k = 0;
    n = 6;
    lambda = alpha^2 * (n + k) - n;

    omega_body = 2*pi/18972.92*[0; 0; 1];
    mass_eros = 6.687e15;

    x0 = [r0; v0];
    options = odeset('reltol', 1e-12, 'abstol', [1e-8*ones(3,1); 1e-11*ones(3,1)]);

    % Generazione sigma points
    [sp] = sigma_points(x0, P0, alpha);
    chi = zeros(n, size(sp, 2));
    
    % Separate angular measurements
    measurements_coalt = measurements_angular.val(1:(length(measurements_angular.val)/2));
    measurements_azim = measurements_angular.val((length(measurements_angular.val)/2 + 1):end);

    measurements_coalt_index = measurements_angular.index;
    measurements_azim_index = measurements_angular.index;
    
    % Initialize variables
    gamma_azim = [];
    gamma_coalt = [];
    gamma_range = [];
    
    %Propagate sigma points
    for i = 1:size(sp, 2)
        [~, propag] = ode78(@(t,x) model_rotating_dynamics(t, x, mass_eros, omega_body, C20, C22), [et_i, et_f], sp(:, i), options);
        chi(:, i) = propag(end, :)';
    end
    
    % Calculating measurements from sigma points
    for i = 1:size(sp, 2)
        [obs_coalt, obs_azim, obs_range] = feature_position_measurement(chi(1:3,i)', chi(4:6,i)', F, V, N, features, et_f, fov1, fov2);
            
            %To avoid mismatching measurements, for each sigma points i
            %check if the measurement is coeherent with the real
            %measurements and, if not, i don't consider the measurement in
            %the filter
            [C, ia, ib] = intersect(obs_coalt.index, measurements_coalt_index, 'stable');
            
            obs_coalt.index = obs_coalt.index(ia);
            obs_coalt.val = obs_coalt.val(ia);
            measurements_coalt_index = measurements_coalt_index(ib);
            measurements_coalt = measurements_coalt(ib);

            obs_azim.index = obs_azim.index(ia);
            obs_azim.val = obs_azim.val(ia);
            measurements_azim_index = measurements_azim_index(ib);
            measurements_azim = measurements_azim(ib);
            
            %Also eliminate the eventual mismatching measurement from the
            %previous sigma points
               if i > 1
                   gamma_coalt = gamma_coalt(ib, :);
                   gamma_azim = gamma_azim(ib, :);
               end
        
        gamma_azim = [gamma_azim, obs_azim.val];
        gamma_coalt = [gamma_coalt, obs_coalt.val];
        gamma_range = [gamma_range, obs_range];
    end

    gamma = [gamma_coalt; gamma_azim; gamma_range];
    
    %If no measurement exists or all the measurements mismatch, perform
    %unfiltered propagation

    if isempty(measurements_coalt)
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

    %UKF
    
    %State variables
    chi0 = chi(:, 1);
    chi_i = chi(:, 2:end);
    xk_minus = 0;

    for j = 1:size(chi_i, 2)
    Wi = 1/(2*(n+lambda));
    xk_minus = xk_minus + Wi*chi_i(:, j);
    end    
    W0m = (lambda)/(lambda+n);
    xk_minus = xk_minus + W0m*chi0;
    
    %Measurements variables
    gamma0 = gamma(:,1);
    gamma_i = gamma(:,2:end);
    yk_minus = 0;

    for j = 1:size(gamma_i, 2)
    Wi = 1/(2*(n+lambda));
    yk_minus = yk_minus + Wi*gamma_i(:, j);
    end
    W0m = (lambda)/(lambda+n);
    yk_minus = yk_minus + W0m*gamma0;
    
    Pk_minus = zeros(length(xk_minus));
    Peek = zeros(length(yk_minus));
    Pxyk = zeros(length(xk_minus), length(yk_minus));

    %Add measurements errors
    yk = [measurements_coalt; measurements_azim; measurement_range];
    variances = [repmat(sigma_angular, size([measurements_coalt; measurements_azim], 1), 1); sigma_range];
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
    
    %Compute final values
    xk_plus = xk_minus + Kk*(yk - yk_minus);
    Pk_plus = Pk_minus - Kk*Peek*Kk';
    
    %Compute the number of features that actually gave a measurements
    n_feature = length(yk)/2-1;
end
