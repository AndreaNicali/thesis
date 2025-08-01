function [nav_score, filt_time] = navigation_score(y, P0, t, spacecraft_data)

sigma_acc = 5e-12;

sigma_meas = ( (100/3600*pi/180)^2 + (100/3600*pi/180)^2 ) ;

[measurements] = feature_measurements(y, t, spacecraft_data);

%Perturb measurements 

for i = 1:length(measurements.coaltitude.val)
    measurements.coaltitude.val(i) = mvnrnd(measurements.coaltitude.val(i), sigma_meas);
    measurements.azimuth.val(i) = mvnrnd(measurements.azimuth.val(i), sigma_meas);
end

%Choose only a certain number of measurements
[C,ia,ib] = intersect(measurements.coaltitude.time, t);

if length(C)>10
    n = length(C);
    idx = round(linspace(1, n, 10));
    ib_sel = ib(idx);
    filt_time = [t(1); t(ib_sel); t(end)];
    filt_time = unique(filt_time, 'stable');   
else
    filt_time = [t(1); t(ib); t(end)];
    filt_time = unique(filt_time, 'stable');   
end

y = reshape(y, [length(t), 6]);

[xx_filtered, P_filtered] = UKF_forsebuono(y(1, :)', filt_time, measurements, P0, spacecraft_data, sigma_meas, sigma_acc);

% for i = 1:length(t)
%     if any(filt_time == t(i))
%         ind = find(filt_time == t(i));
%         nav_score(i) = 1/det(P_filtered(:, :, ind));
%     else
%         ind_min = find(filt_time == max(filt_time( filt_time < t(i) ) ));
%         ind_max = find(filt_time == min(filt_time( filt_time > t(i) ) ));
% 
%         det_min = det(P_filtered(:, :, ind_min));
%         det_max = det(P_filtered(:, :, ind_max));
% 
% end

nav_score = zeros(size(t));
nav_score(end) = log10( 1/det(P_filtered(1:3, 1:3, end)) );

end
