function [P_filtered, xx_filtered, eta_f, flag] = navigationFilter(y0, y_truth, P0, tt, spacecraft_data)
%Flag is 0 if the navigation hasn't been performed at that instant, 1 if it
%has

sigma_acc = spacecraft_data.data_guidance.process_noise;
sigma_meas = spacecraft_data.data_guidance.measurement_noise;

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

eta_start = spacecraft_data.data_guidance.eta0;
[xx, Pp, flag] = EKF_augmented(y0, eta_start, tt, measurements, P0, spacecraft_data, sigma_meas, sigma_acc);
xx_filtered = xx(:, 1:6);
eta_f = xx(end, 7:9);
P_filtered = Pp;
flag = [0, flag]; 

end
