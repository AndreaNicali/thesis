function [P_filtered, filt_time] = navigation_for_MCTS(y, P0, t, spacecraft_data)

sigma_acc = 5e-9;

sigma_meas = ( (90/3600*pi/180)^2 + (90/3600*pi/180)^2 ) ;

[measurements] = feature_measurements(y, t, spacecraft_data);

%Perturb measurements 

for i = 1:length(measurements.coaltitude.val)
    measurements.coaltitude.val(i) = mvnrnd(measurements.coaltitude.val(i), sigma_meas);
    measurements.azimuth.val(i) = mvnrnd(measurements.azimuth.val(i), sigma_meas);
end

%Choose only a certain number of measurements
[C,ia,ib] = intersect(measurements.coaltitude.time, t);

nmax = length(t)/2;

if length(C)>nmax
    n = length(C);
    idx = round(linspace(1, n, nmax));
    ib_sel = ib(idx);
    filt_time = [t(1); t(ib_sel); t(end)];
    filt_time = unique(filt_time, 'stable');   
else
    filt_time = [t(1); t(ib); t(end)];
    filt_time = unique(filt_time, 'stable');   
end

y = reshape(y, [length(t), 6]);

[xx_filtered, P_filtered] = UKF_forsebuono(y(1, :)', filt_time, measurements, P0, spacecraft_data, sigma_meas, sigma_acc);



end
