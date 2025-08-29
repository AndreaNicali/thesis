function [P_filtered, xx_filtered, flag] = navigationFilter(y0, y_truth, P0, tt, spacecraft_data)

sigma_acc = 5e-9;

sigma_meas = sqrt( (100/3600*pi/180)^2 + (100/3600*pi/180)^2 ) ;

[measurements] = measurementsFun(y_truth, tt, spacecraft_data);

%Perturb measurements 
pert_meas = zeros(size(measurements.val));
    for i = 1:size(measurements.val, 2)
        elevation = mvnrnd(measurements.val(1,i), sigma_meas^2);
        azimuth = mvnrnd(measurements.val(2,i), sigma_meas^2);
        pert_meas(:, i) = [elevation;azimuth];
    end
measurements.val = pert_meas;

%reshape trajectory for consistency
y_truth = reshape(y_truth, [length(tt), 6]);

[xx_filtered, P_filtered, flag] = EKF(y0, tt, measurements, P0, spacecraft_data, sigma_meas, sigma_acc, y_truth);
flag = [0, flag];

end
