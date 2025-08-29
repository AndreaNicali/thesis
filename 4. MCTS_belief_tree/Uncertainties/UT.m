function [mean_UT_f, P_UT_f] = UT(r0, v0, P0, et_i, et_f, C20, C22)

%Set Unscented Transform parameters
alpha = 1;
beta = 2;
k = 0;
n = 6;
mean_UT_f = zeros(6, 1);
P_UT_f = zeros(6,6);

%Create sigma points
[sp] = sigma_points([r0; v0], P0, alpha);
y = zeros(size(sp));

omega_body = 2*pi/18972.92*[0; 0; 1];
mass_eros = 6.687e15;

%Start Unscented Transform
options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

for i = 1:size(sp, 2)

    %Propagate each sigma points
    [~, propag] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), [et_i, et_f], sp(:, i), options);
    y(:, i) = propag(end, :)';

end

y0 = y(:,1);
yi = y(:,2:end);

%Compute weights for the mean value
for i = 1:size(yi, 2)
    Wi = 1/(2*(n+k));
    mean_UT_f = mean_UT_f + Wi*yi(:, i);
end

W0m = 1 - n/(alpha^2*(n+k));

%Compute mean using weighted propagations
mean_UT_f = mean_UT_f + W0m*y0;

%Compute weights for covariance matrices
for i = 1:size(yi, 2)
    Wi = 1/(2*(n+k));
    P_UT_f = P_UT_f + Wi*(yi(:, i) - mean_UT_f)*(yi(:, i) - mean_UT_f)';
end
W0c = (2-alpha^2+beta) - n/(alpha^2*(n+k));

%Compute final covariance matrix
P_UT_f = P_UT_f + W0c*(y0 - mean_UT_f)*(y0 - mean_UT_f)';

end