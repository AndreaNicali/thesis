function [xx_filtered, P_filtered] = EKF(x0, tt, measurements, P0, spacecraft_data, sigma_meas, sigma_acc)

C20 = spacecraft_data.data_asteroids.C20;
C22 = spacecraft_data.data_asteroids.C22;

mass_eros = spacecraft_data.data_asteroids.mass;
omega_body = spacecraft_data.data_asteroids.omega;

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
xx_filtered = zeros(length(tt), 6);
P_filtered = zeros(6, 6, length(tt));
xx_filtered(1, :) = x0;
P_filtered(:, :, 1) = P0;

x_start = x0;
P_start = P0;

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

    meas = measurements.val(measurements.time == tt(i+1));
    meas = meas(:);

    [~, xx, PHI]  = propagateEllipsoid(x0, [tt(i), tt(i+1)], mass_eros, omega_body, C20, C22);

    sim_meas = measurementsFun(xx(end,:), tt(i+1), spacecraft_data);
    sim_meas = sim_meas.val(sim_meas.time == (tti+1));
    sim_meas = sim_meas(:);
    [H, R] = measurementsMatrix(sim_meas, sigma_meas);

    xxkplus1 = xx(end, :);
    Pkplus1 = PHI*Pt*PHI' + Q;
    bkplus1 = meas - H*xxkplus1';

    Kkplus1 = Pkplus1*H'/(H*Pkplus1*H'+R);

    P_upd = Pkplus1-Kkplus1*H*Pkplus1;
    xx_upd = xxkplus1 + Kkplus1*bkplus1;
    
    xx_filtered(i+1, :) = xx_upd;
    P_filtered(:, :, i+1) = P_upd;

    x_start = xx_upd;
    P_start = P_upd;
end

   
end

function [H, R] = measurementsMatrix(sim_meas, sigma_meas)
H = zeros(length(sim_meas), 6);
for i = 1:(length(sim_meas)/3)
    p = sim_meas( (3*(i-1)+1):(3*(i+1)+3 ));
    H((3*(i-1)+1):(3*(i+1)+3 ), 1:3) = (- eye(3)/norm(p)-p*p'/norm(p)^3);

    az = atan2(p(2), p(1));
    el = asin(p(3));
    J = [-cos(el)*sin(az), -sin(el)*cos(az);
        cos(el)*cos(az), -sin(el)*sin(az);
            0        ,       cos(el)];

    R = sigma_meas*(J*J');


end
end
