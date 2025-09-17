function [dvCorr1, dvCorr2] = refTargeting(x0, tt, refTraj, spacecraft_data)

C20 = spacecraft_data.data_asteroids.C20;
C22 = spacecraft_data.data_asteroids.C22;
mass = spacecraft_data.data_asteroids.mass;
omega = spacecraft_data.data_asteroids.omega;


x0Targ = refTraj(1, :)';
xfTarg = refTraj(end, :)';
maxIter = 100;
tol = 10^-6;
diff = 1;
iter = 0;
xCorr = x0;
Corr = zeros(3, 1);

while norm(diff) > tol && iter < maxIter

    xCorr = xCorr + [zeros(3,1); Corr];

    [~, xx, PHIf] = propagateModel(xCorr, tt, mass, omega, C20, C22);
    xf = xx(end, :)';
    
    PHIrv = PHIf(1:3,4:6);
    PHIrr = PHIf(1:3,1:3);
    
    diff = xfTarg(1:3)-xf(1:3);
    
    %Differential Correction
    Corr = PHIrv\diff;
    iter = iter+1;
end

dvCorr1 = xCorr(4:6) - x0(4:6);
dvCorr2 = xfTarg(4:6) - xf(4:6);

end




