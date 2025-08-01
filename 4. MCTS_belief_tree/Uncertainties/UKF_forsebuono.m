function [xx_filtered, P_filtered] = UKF_forsebuono(x0, tt, measurements, P0, spacecraft_data, sigma_meas, sigma_acc)

C20 = spacecraft_data.data_asteroids.C20;
C22 = spacecraft_data.data_asteroids.C22;

mass_eros = spacecraft_data.data_asteroids.mass;
omega_body = spacecraft_data.data_asteroids.omega;

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

alpha = 0.01;
beta = 2;
n = length(x0);
k = 0;
lambda = alpha^2*(n+k) - n;

xx_filtered = [];
P_filtered = P0;

for i = 1:(length(tt)-1)
    
    tao = tt(i+1)-tt(i);
    q11 = (1/4)*tao^4;
    q12 = (1/2)*tao^3;
    q22 = tao^2;

    Qk = sigma_acc^2 * [ ...
        q11, 0,   0,   q12, 0,   0;
        0,   q11, 0,   0,   q12, 0;
        0,   0,   q11, 0,   0,   q12;
        q12, 0,   0,   q22, 0,   0;
        0,   q12, 0,   0,   q22, 0;
        0,   0,   q12, 0,   0,   q22];

    meas_azim = measurements.azimuth.val(measurements.azimuth.time == tt(i+1));
    meas_coal = measurements.coaltitude.val(measurements.coaltitude.time == tt(i+1));

    feat_azim = measurements.azimuth.feature(measurements.azimuth.time == tt(i+1));
    feat_coal = measurements.coaltitude.feature(measurements.coaltitude.time == tt(i+1));
    
    nmax = 10;
    if length(meas_azim) > nmax
        meas_azim = meas_azim(1:nmax);
        meas_coal = meas_coal(1:nmax);
        feat_azim = feat_azim(1:nmax);
        feat_coal = feat_coal(1:nmax);
    end

    %GENERATE SIGMA POINTS
    [sp] = sigma_points(x0, P0+Qk, alpha);

    chi = zeros(size(sp));
    gamma = [];

    %GENERATE PROPAGATED POSITIONS AND MEASURMENTS
    for j = 1:size(sp, 2)

        [~ , prop] = ode78(@(t,x) model_rotating_dynamics(t, x, mass_eros, omega_body, C20, C22), [tt(i), tt(i+1)], sp(:, j), options);
        chi(:, j) = prop(end, :);
            
        sim_measurements = feature_measurements(chi(:, j)', tt(i+1), spacecraft_data);
        
        sim_meas_azim = sim_measurements.azimuth.val;
        sim_meas_coal = sim_measurements.coaltitude.val;

        sim_feat_azim = sim_measurements.azimuth.feature;
        sim_feat_coal = sim_measurements.coaltitude.feature;        
        
        %MEASUREMENT ALLIGNMENT
        if j == 1
            [common, ia, ib] = intersect(sim_feat_coal, feat_coal);
            
            sim_feat_azim = sim_feat_azim(ia);
            sim_meas_azim = sim_meas_azim(ia);
            sim_feat_coal = sim_feat_coal(ia);
            sim_meas_coal = sim_meas_coal(ia);

            feat_azim = feat_azim(ib);
            meas_azim = meas_azim(ib);
            meas_coal = meas_coal(ib);
            feat_coal= feat_coal(ib);

            skip_this_timestep = false;

            if any(~isreal(sp))
                a = 1;
            end
            
            if isempty(sim_meas_azim)
            warning("Nessuna feature visibile comune al tempo %f", tt(i+1));
            [x, P] = UT(x0(1:3), x0(4:6), P0, tt(i), tt(i+1), C20, C22);
            xx_filtered = [xx_filtered; x'];
            P_filtered(:, :, i+1) = P;

            skip_this_timestep = true;
            break;
            end

            gamma(:, j) = [ sim_meas_coal'; sim_meas_azim'];
            valid_ind = sim_feat_coal';

        else
            [common, ia, ib] = intersect(sim_feat_coal, feat_coal);
            
            sim_feat_azim = sim_feat_azim(ia);
            sim_meas_azim = sim_meas_azim(ia);
            sim_feat_coal = sim_feat_coal(ia);
            sim_meas_coal = sim_meas_coal(ia);

            feat_azim = feat_azim(ib);
            meas_azim = meas_azim(ib);
            meas_coal = meas_coal(ib);
            feat_coal = feat_coal(ib);

            skip_this_timestep = false;
            if isempty(sim_meas_azim)
            warning("Nessuna feature visibile comune al tempo %f", tt(i+1));
            [x, P] = UT(x0(1:3), x0(4:6), P0, tt(i), tt(i+1), C20, C22);
            xx_filtered = [xx_filtered; x'];
            P_filtered(:, :, i+1) = P;

            skip_this_timestep = true;
            break;
            end
            
            [common, ia, ib] = intersect(valid_ind, sim_feat_coal');
            
            gamma1 = gamma(1:(size(gamma, 1)/2), :);
            gamma2 = gamma((size(gamma, 1)/2 + 1):end, :);
            
            gamma1 = gamma1(ia, :);
            gamma2 = gamma2(ia, :);

            gamma = [gamma1; gamma2];

            sim_feat_azim = sim_feat_azim(ib);
            sim_meas_azim = sim_meas_azim(ib);
            sim_feat_coal = sim_feat_coal(ib);
            sim_meas_coal = sim_meas_coal(ib);
            feat_azim = feat_azim(ib);
            meas_azim = meas_azim(ib);
            meas_coal = meas_coal(ib);
            feat_coal = feat_coal(ib);

            valid_ind = sim_feat_coal';
            gamma(:, j) = [ sim_meas_coal'; sim_meas_azim'];
        end
    end
    
    if skip_this_timestep
        x0 = x;
        P0 = P;
        continue;
    end
    
    if isempty(gamma)
    warning("Nessuna feature visibile comune al tempo %f", tt(i+1));
    [x, P] = UT(x0(1:3), x0(4:6), P0, tt(i), tt(i+1), C20, C22);
    xx_filtered = [xx_filtered; x'];
    P_filtered(:, :, i+1) = P;

    x0 = x;
    P0 = P;
    continue;
    end
    
    meas_azim = wrapTo2Pi(meas_azim);
    meas_coal = wrapTo2Pi(meas_coal);
    gamma(1:end/2, :) = wrapTo2Pi(gamma(1:end/2, :));
    gamma(end/2+1:end, :) = wrapTo2Pi(gamma(end/2+1:end, :)); 

    Rk = eye(size(gamma, 1))*sigma_meas;
    meas_vec = [meas_coal'; meas_azim'];

    %FIND MINUS VALUES FOR POSITIONS
    chi0 = chi(:, 1);
    chi_i = chi(:, 2:end);
    x_minus = 0;

    for j = 1:size(chi_i, 2)
    Wi = 1/(2*(n+lambda));
    x_minus = x_minus + Wi*chi_i(:, j);
    end    
    W0m = (lambda)/(lambda+n);
    x_minus = x_minus + W0m*chi0;
    
    %FIND MINUS VALUES FOR MEASUREMENTS 
    gamma0 = gamma(:,1);
    gamma_i = gamma(:,2:end);
    y_minus = 0;

    for j = 1:size(gamma_i, 2)
    Wi = 1/(2*(n+lambda));
    y_minus = y_minus + Wi*gamma_i(:, j);
    end
    W0m = (lambda)/(lambda+n);
    y_minus = y_minus + W0m*gamma0;
    
    %FIND COVARIANCE MATRICES VALUES
    Pk_minus = zeros(6,6);
    Peek = zeros(size(gamma, 1));
    Pxyk = zeros(6, size(gamma, 1));

    for j = 1:size(gamma_i, 2)
    Wic = 1/(2*(n+lambda));
    Pk_minus = Pk_minus + Wic*(chi_i(:, j) - x_minus)*(chi_i(:, j) - x_minus)';
    Peek = Peek + Wic*(gamma_i(:, j) - y_minus)*(gamma_i(:, j) - y_minus)';
    Pxyk = Pxyk + Wic*(chi_i(:, j) - x_minus)*(gamma_i(:, j) - y_minus)';

    end
    %Pk_minus = Pk_minus + Qk;
    Peek = Peek + Rk;
    W0c = (lambda)/(lambda+n) + (1-alpha^2+beta);
    Pk_minus = Pk_minus + W0c*(chi0 - x_minus)*(chi0 - x_minus)';
    Peek = Peek + W0c*(gamma0 - y_minus)*(gamma0 - y_minus)';
    Pxyk = Pxyk + W0c*(chi0 - x_minus)*(gamma0 - y_minus)';

    Kk = Pxyk/Peek;
    
    %COMPUTE "PLUS" VALUES
    x = x_minus + Kk*wrapToPi(meas_vec - y_minus);
    P = Pk_minus - Kk*Peek*Kk';


    x0 = x;
    P0 = P;

    xx_filtered = [xx_filtered; x'];
    P_filtered(:, :, i+1) = P;
end

end



