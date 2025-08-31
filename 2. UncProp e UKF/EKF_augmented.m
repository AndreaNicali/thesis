function [xx_filtered, P_filtered, flag] = EKF_augmented(x0, eta0, tt, measurements, P0, spacecraft_data, sigma_meas, sigma_acc, xx_truth)

C20 = spacecraft_data.data_asteroids.C20;
C22 = spacecraft_data.data_asteroids.C22;

mass_eros = spacecraft_data.data_asteroids.mass;
omega_body = spacecraft_data.data_asteroids.omega;

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
xx_filtered = zeros(length(tt), 9);
P_filtered = zeros(9, 9, length(tt));
xx_filtered(1, :) = [x0, eta0];
P_filtered(:, :, 1) = P0;

x_start = [x0, eta0];
P_start = P0;
flag = [];

for i = 1:(length(tt)-1)

    tao = 24*3600;
    
    meas_t = measurements.time == tt(i+1);
    meas = measurements.val(:, meas_t);

    [~, xx, PHI]  = propagateEllipsoidEKF(x_start', [tt(i), tt(i+1)], mass_eros, omega_body, C20, C22, tao);
    
    % delta = [x_start(1)*10^-8, x_start(2)*10^-8, x_start(3)*10^-8, 10^-8, 10^-8, 10^-8];
    % eps_pert = diag(delta);
    % perturbed_states = [x_start', x_start', x_start', x_start', x_start', x_start']+eps_pert;
    % 
    % perturbed_prop = zeros(size(perturbed_states));
    % PHI_num = zeros(6,6);
    % 
    % for j = 1:6
    %     [~, xx_p, PHI]  = propagateEllipsoid(perturbed_states(:, j), [tt(i), tt(i+1)], mass_eros, omega_body, C20, C22);
    %     perturbed_prop(:, j) = xx_p(end, :)';
    % 
    %     PHI_num(:, j) = ( perturbed_prop(:, j)-xx(end, :)' )/delta(j);
    % end
    % PHI = PHI_num;

    sim_measurements = measurementsFun(xx(end,:), tt(i+1), spacecraft_data);
    sim_meas_t = sim_measurements.time == tt(i+1);
    sim_meas = sim_measurements.val(:, sim_meas_t);

    if any(size(sim_meas)~=size(meas))
        sim_features = sim_measurements.feature(sim_meas_t);
        real_features = measurements.feature(meas_t);
        [comm, ia, ib] = intersect(sim_features, real_features);

        sim_meas = sim_meas(:, ia);
        meas = meas(:, ib);
        sim_meas_t = false(size(sim_meas_t));
        sim_meas_t(ia) = true;

    end
    dt  = tt(i+1)-tt(i);
    tau = 24*3600;                      % come avevi
    q = 2*sigma_acc^2 / tau;            % PSD per DMC
    Qk = Qd_DMC(dt, tau, q, eye(3));    % oppure M_RTN2I(t) se eta in RTN

    sim_meas = sim_meas(:);
    meas = meas(:);

    los_meas = sim_measurements.los(:, sim_meas_t);

    [H, R] = measurementsMatrix(los_meas, sigma_meas);

    xxkplus1 = xx(end, :);
    Pkplus1 = PHI*P_start*PHI' + Qk;

   if isempty(meas)
       xx_filtered(i+1,:) = xxkplus1;
       P_filtered(:, :, i+1) = Pkplus1;
       x_start = xxkplus1;
       P_start = Pkplus1;
       flag = [flag, 0];
       continue
   end

    bkplus1 = wrapToPi(meas - sim_meas);

    Kkplus1 = Pkplus1*H'/(H*Pkplus1*H'+ R);

    I = eye(size(Pkplus1));
    P_upd = (I - Kkplus1*H)*Pkplus1*(I - Kkplus1*H)' + Kkplus1*R*Kkplus1';
    xx_upd = xxkplus1' + Kkplus1*bkplus1;
    
    xx_filtered(i+1, :) = xx_upd';
    P_filtered(:, :, i+1) = P_upd;
    
    flag = [flag, 1];
    x_start = xx_upd';
    P_start = P_upd;

end

   
end


function [H, R] = measurementsMatrix(los_meas, sigma_meas)

N = size(los_meas, 2);   % numero di misure
H = zeros(2*N, 9);

for i = 1:N
    rho = norm(los_meas(:, i));
    rhox = los_meas(1,i)/rho;
    rhoy = los_meas(2,i)/rho;
    rhoz = los_meas(3,i)/rho;
    
    % Jacobiano 2x3 rispetto alla posizione
    J = [ 1/sqrt(1-rhoz^2)*rhox*rhoz/rho,   1/sqrt(1-rhoz^2)*rhoy*rhoz/rho,  1/sqrt(1-rhoz^2)*(-1+rhoz^2)/rho; 
        1/(1+(rhoy/rhox)^2)*rhoy/(rhox^2*rho),  -1/(1+(rhoy/rhox)^2)*1/(rhox*rho),  0 ];
    
    % Inserisco nel blocco 2x6 (velocità non influiscono)
    H(2*(i-1)+1 : 2*i, 1:3) = J;
end

% Covarianza misure (assunte indipendenti e identiche)
R = sigma_meas^2 * eye(2*N);

end

function Qk = Qd_DMC(dt, tau, q_vec, M)
% q_vec: [qx qy qz] intensità PSD su eta (km^2/s^5)
% M: rotazione (es. RTN->I). Usa eye(3) se già in inerziale.
if isscalar(q_vec), q_vec = [q_vec q_vec q_vec]; end
Qtil = M*diag(q_vec)*M.';   % 3x3

gpp = tau^5/2*((1-exp(-2*dt/tau)) + 2*dt/tau*(1-2*exp(-dt/tau)) - 2*(dt/tau)^2 + 2/3*(dt/tau)^3);
gpv = tau^4/2*((exp(-2*dt/tau)-1) - 2*(exp(-dt/tau)-1) + 2*dt/tau*(exp(-dt/tau)-1) + (dt/tau)^2);
gpe = tau^3/2*((1-exp(-2*dt/tau)) - 2*dt/tau*exp(-dt/tau));
gvv = tau^3/2*((1-exp(-2*dt/tau)) - 4*(1-exp(-dt/tau)) + 2*dt/tau);
gve = tau^2/2*(1-exp(-dt/tau))^2;
gee = tau^2*(1-exp(-2*dt/tau));

S = [gpp*Qtil, gpv*Qtil, gpe*Qtil;
     gpv*Qtil, gvv*Qtil, gve*Qtil;
     gpe*Qtil, gve*Qtil, gee*Qtil];

Qk = zeros(9);
Qk([1:3,4:6,7:9],[1:3,4:6,7:9]) = S;
end

