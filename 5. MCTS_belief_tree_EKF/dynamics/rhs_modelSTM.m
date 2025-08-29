function [acc, G] = rhs_modelSTM(y, mass, C20, C22)

% Extract position and velocity from state
r = y(1:3);

% Set Constants and useful quantities
Ggrav = astroConstants(1);
R0 = 16; %da Determination of Shape, Gravity, and Rotational State of Asteroid 433 Eros
rnorm = norm(r);
lambda = atan2(y(2), y(1));
phi = asin(y(3)/rnorm);
eta = sqrt(y(1)^2+y(2)^2);

%to avoid eta = 0
if eta < 1e-10
    eta = 1e-10; 
end
mu = Ggrav*mass;

%Initialize derivatives that then will compose the carthesian accelerations
dUdr = 0;
dUdlam = 0;
dUdphi = 0;

dUdr2 = 0;
dUdrdphi = 0;
dUdrdlam = 0;
dUdphi2 = 0;
dUdphidlam = 0;
dUdlam2 = 0;

%Start the sum: n goes from 1 ( n=0 terms is accounted as GM/rnorm^3 *r ) to
%length C minus 1 (remember coefficient order has been increased by one). m
%goes from 0 to n.

order = 2;

leg_pol_data = leg_pol_rec(sin(phi), order);    

for n = 1:order

    for m = 0:(n)

        if n == 2 && m == 0
        Cbar = C20;
        Sbar = 0;
        else 
            if n == 2 && m == 2
            Cbar = C22;
            Sbar = 0;
            else
            Cbar = 0;
            Sbar = 0;
            end
        end

            %Evaluate normalized legendre polynomial in sin(phi)
            ind = find(leg_pol_data(:, 1) == n & leg_pol_data(:, 2) == m);
            leg_pol = leg_pol_data(ind, 3);
            
            %Update current value in the summation
            dUdr_plus = (R0/rnorm)^n * (n+1) * leg_pol * (Cbar*cos(m*lambda) + Sbar*sin(m*lambda));
            dUdr = dUdr + dUdr_plus;
            dUdlam_plus = (R0/rnorm)^n  * leg_pol * m * (-Cbar*sin(m*lambda) + Sbar*cos(m*lambda));
            dUdlam = dUdlam + dUdlam_plus;

            
            
            %Compute the derivative term of dUdphi
            if m > 0
                Knm = sqrt((n-m)*(n+m+1));
            else
                Knm = sqrt(n*(n+1)/2);
            end
            Knmplus1 = sqrt( (n-m-1)*(n+m+2) );
            %Impose that legendre polynomial with P(n,m+1) = 0 if m == n
            if m < n
                ind = find(leg_pol_data(:, 1) == n & leg_pol_data(:, 2) == m+1);
                leg_pol_plus1 = leg_pol_data(ind, 3);
            else
                leg_pol_plus1 = 0;
            end

            if m+1 < n
                ind2 = find(leg_pol_data(:, 1) == n & leg_pol_data(:, 2) == m+2);
                leg_pol_plus2 = leg_pol_data(ind2, 3);
            else
                leg_pol_plus2 = 0;
            end

            dPnmdphi = ( -m*tan(phi)*leg_pol + Knm*leg_pol_plus1 );
            %dPnmdphi2 = -m*sec(phi)*leg_pol - tan(phi)*dPnmdphi - Knm*(m+1)*tan(phi)*leg_pol_plus1 + Knm*Knmplus1*leg_pol_plus2;
            dPnmdphi2 = Knmplus1 * leg_pol_plus2 - (2*m+1)*tan(phi)*Knm*leg_pol_plus1 + m*(m*tan(phi)^2 - sec(phi)^2)*leg_pol;
                        
            %Update last derivative value
            dUdphi_plus = (R0/rnorm)^n * dPnmdphi * (Cbar*cos(m*lambda) + Sbar*sin(m*lambda));
            dUdphi = dUdphi + dUdphi_plus;
            
            dUdr2 = dUdr2 + dUdr_plus*(n+2);
            dUdrdlam = dUdrdlam + (n+1)*(R0/rnorm)^n*leg_pol*(-Cbar*sin(m*lambda) + Sbar*cos(m*lambda));
            dUdrdphi = dUdrdphi + dUdphi_plus*(n+1);
            dUdphi2 = dUdphi2 + dUdphi_plus / dPnmdphi * dPnmdphi2;
            dUdphidlam = dUdphidlam + dUdlam_plus/leg_pol*dPnmdphi;
            dUdlam2 = dUdlam2 + dUdr_plus/(n+1)*m^2;
            
    end
end


    %Final moltiplication
    dUdr = -mu/rnorm^2*(1+dUdr);
    dUdlam = mu/rnorm*dUdlam;
    dUdphi = mu/rnorm*dUdphi;

    dUdr2 = mu/rnorm^3*(2+dUdr2);
    dUdrdphi = -mu/rnorm^2*dUdrdphi;
    dUdrdlam = -mu/rnorm^2*dUdrdlam;
    dUdphi2 = mu/rnorm*dUdphi2;
    dUdphidlam = mu/rnorm*dUdphidlam;
    dUdlam2 = -mu/rnorm*dUdlam2;
    
    
    %Building accelerations in carthesian frame
    acc_x = ( (1/rnorm * dUdr - y(3)/(rnorm^2*eta) * dUdphi)*y(1) - (1/eta^2 * dUdlam)*y(2) );
    acc_y = ( (1/rnorm * dUdr - y(3)/(rnorm^2*eta) * dUdphi)*y(2) + (1/eta^2 * dUdlam)*y(1) );
    acc_z = ( 1/rnorm * dUdr * y(3) + eta/rnorm^2 * dUdphi );
    
    acc = [acc_x; acc_y; acc_z];

    rx = y(1);
    ry = y(2);
    rz = y(3);
    vx = y(4);
    vy = y(5);
    vz = y(6);

    drdx = rx/rnorm;
    drdy = ry/rnorm;
    drdz = rz/rnorm;
    dphidx = -rx*rz/(rnorm^2*eta);
    dphidy = -ry*rz/(rnorm^2*eta);
    dphidz = (1-rz^2/rnorm^2)/eta;
    dlamdx = -ry/eta^2;
    dlamdy = rx/eta^2;
    dlamdz = 0;

    drdx2 = (ry^2+rz^2) / rnorm^3;
    drdxdy = -ry*rx/ rnorm^3;
    drdxdz = -rz*rx/rnorm^3;
    drdy2 = (rx^2+rz^2)/rnorm^3;
    drdydz = -ry*rz/ rnorm^3;
    drdz2 = (rx^2 + ry^2)/rnorm^3;


    dphidx2 = (-rnorm^2*eta*rz + rx*rz*(2*rx*eta+rnorm^2*rx/eta)) / (rnorm^4*eta^2);
    dphidxdy = (rx*rz*(2*ry*eta+rnorm^2*ry/eta))/(rnorm^4*eta^2);
    dphidxdz = (-rx*rnorm^2*eta + rx*rz*2*rz*eta)/(rnorm^4*eta^2);
    dphidy2 = (-rnorm^2*eta*rz + ry*rz*(2*ry*eta+rnorm^2*ry/eta)) / (rnorm^4*eta^2);
    dphidydz = (-ry*rnorm^2*eta + ry*rz*2*rz*eta)/(rnorm^4*eta^2);
    dphidz2 = - (2*rz*rnorm^2 - 2*rz^3)/(rnorm^4*eta);


    dlamdx2 = (2*rx*ry) / eta^4;
    dlamdxdy = (-eta^2 + 2*ry^2)/eta^4;
    dlamdxdz = 0;
    dlamdy2 = (-2*rx*ry) / eta^4;
    dlamdydz = 0;
    dlamdz2 = 0;

    
    G11 = drdx^2*dUdr2 + dphidx^2*dUdphi2 + dlamdx^2*dUdlam2 + ...
        2*drdx*dphidx*dUdrdphi + 2*drdx*dlamdx*dUdrdlam + 2*dphidx*dlamdx*dUdphidlam + ...
        dUdr*drdx2 + dUdphi*dphidx2 + dUdlam*dlamdx2;

    G12 = drdy*drdx*dUdr2 + drdy*dphidx*dUdrdphi + drdy*dlamdx*dUdrdlam + ...
        dphidy*drdx*dUdrdphi + dphidy*dphidx*dUdphi2 + dphidy*dlamdx*dUdphidlam + ...
        dlamdy*drdx*dUdrdlam + dlamdy*dphidx*dUdphidlam + dlamdy*dlamdx*dUdlam2 + ...
        dUdr*drdxdy + dUdphi*dphidxdy + dUdlam*dlamdxdy;

    G13 = drdz*drdx*dUdr2 + drdz*dphidx*dUdrdphi + drdz*dlamdx*dUdrdlam + ...
        dphidz*drdx*dUdrdphi + dphidz*dphidx*dUdphi2 + dphidz*dlamdx*dUdphidlam + ...
        dlamdz*drdx*dUdrdlam + dlamdz*dphidx*dUdphidlam + dlamdz*dlamdx*dUdlam2 + ...
        dUdr*drdxdz + dUdphi*dphidxdz + dUdlam*dlamdxdz;

    % G21 = drdy*drdx*dUdr2 + drdx*dphidy*dUdrdphi + drdx*dlamdy*dUdrdlam + ...
    %     dphidx*drdy*dUdrdphi + dphidy*dphidx*dUdphi2 + dphidx*dlamdy*dUdphidlam + ...
    %     dlamdx*drdy*dUdrdlam + dlamdx*dphidy*dUdphidlam + dlamdy*dlamdx*dUdlam2 + ...
    %     dUdr*drdxdy + dUdphi*dphidxdy + dUdlam*dlamdxdy;

    G22 = drdy^2*dUdr2 + dphidy^2*dUdphi2 + dlamdy^2*dUdlam2 + ...
        2*drdy*dphidy*dUdrdphi + 2*drdy*dlamdy*dUdrdlam + 2*dphidy*dlamdy*dUdphidlam + ...
        dUdr*drdy2 + dUdphi*dphidy2 + dUdlam*dlamdy2;
    
    G23 = drdz*drdy*dUdr2 + drdz*dphidy*dUdrdphi + drdz*dlamdy*dUdrdlam + ...
        dphidz*drdy*dUdrdphi + dphidz*dphidy*dUdphi2 + dphidz*dlamdy*dUdphidlam + ...
        dlamdz*drdy*dUdrdlam + dlamdz*dphidy*dUdphidlam + dlamdz*dlamdy*dUdlam2 + ...
        dUdr*drdydz + dUdphi*dphidydz + dUdlam*dlamdydz;

    % G31 = drdz*drdx*dUdr2 + drdx*dphidz*dUdrdphi + drdx*dlamdz*dUdrdlam + ...
    %     dphidx*drdz*dUdrdphi + dphidx*dphidz*dUdphi2 + dphidx*dlamdz*dUdphidlam + ...
    %     dlamdx*drdz*dUdrdlam + dlamdx*dphidz*dUdphidlam + dlamdx*dlamdz*dUdlam2 + ...
    %     dUdr*drdxdz + dUdphi*dphidxdz + dUdlam*dlamdxdz;

    % G32 = drdz*drdy*dUdr2 + drdz*dphidy*dUdrdphi + drdz*dlamdy*dUdrdlam + ...
    %     dphidz*drdy*dUdrdphi + dphidz*dphidy*dUdphi2 + dphidz*dlamdy*dUdphidlam + ...
    %     dlamdz*drdy*dUdrdlam + dlamdz*dphidy*dUdphidlam + dlamdz*dlamdy*dUdlam2 + ...
    %     dUdr*drdydz + dUdphi*dphidydz + dUdlam*dlamdydz;

    G33 = drdz^2*dUdr2 + dphidz^2*dUdphi2 + dlamdz^2*dUdlam2 + ...
        2*drdz*dphidz*dUdrdphi + 2*drdz*dlamdz*dUdrdlam + 2*dphidz*dlamdz*dUdphidlam + ...
        dUdr*drdz2 + dUdphi*dphidz2 + dUdlam*dlamdz2;

    G = [G11, G12, G13;
        G12, G22, G23;
        G13, G23, G33];


end