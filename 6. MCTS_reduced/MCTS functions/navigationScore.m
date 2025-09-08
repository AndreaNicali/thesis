function [nav_score] = navigationScore(y, P0, t, spacecraft_data)

omega = spacecraft_data.data_asteroids.omega(3);
minLandToScore = spacecraft_data.data_guidance.minLandmToScore;
DU = spacecraft_data.data_guidance.DU;
detPRef = spacecraft_data.data_guidance.detPRef;

[measurements] = measurementsFun(y, t, spacecraft_data);

nav_bool = zeros(size(t));
for i = 1:length(t)
    list = find(measurements.time == t(i));
    if length(list)>= minLandToScore
        nav_bool(i) = 1/length(t);
    end
end
  
TU = sqrt( DU^3/(astroConstants(1)*spacecraft_data.data_asteroids.mass) );
VU = DU/TU + omega*DU;

P_adim = zeros(size(P0(1:6, 1:6)));
P_adim(1:3, 1:3) = P0(1:3, 1:3)/(DU*DU);
P_adim(4:6, 1:3) = P0(4:6, 1:3)/(DU*VU);
P_adim(1:3, 4:6) = P0(1:3, 4:6)/(DU*VU);
P_adim(4:6, 4:6) = P0(4:6, 4:6)/(VU*VU);

if log10(det(P_adim)) < detPRef
    nav_score = nav_bool*0;
else
    nav_score = 2*(1 - ( log10(det(P_adim) ) / detPRef) ) * nav_bool;
end

end