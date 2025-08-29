function acc = Recursive_Spherical_armonics_perturbation_ellipsoid(y, mass, C20, C22)

% Extract position and velocity from state
r = y(1:3);

% Set Constants and useful quantities
G = astroConstants(1);
R0 = 16; %da Determination of Shape, Gravity, and Rotational State of Asteroid 433 Eros
rnorm = norm(r);
if ~isreal(y)
    ve = 2;
end
lambda = atan2(y(2), y(1));
phi = asin(y(3)/rnorm);
eta = sqrt(y(1)^2+y(2)^2);

%to avoid eta = 0
if eta < 1e-10
    eta = 1e-10; 
end

%Initialize derivatives that then will compose the carthesian accelerations
dUdr = 0;
dUdlam = 0;
dUdphi = 0;

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
            dUdr = dUdr + (R0/rnorm)^n * (n+1) * leg_pol * (Cbar*cos(m*lambda) + Sbar*sin(m*lambda));
            dUdlam = dUdlam + (R0/rnorm)^n  * leg_pol * m * (-Cbar*sin(m*lambda) + Sbar*cos(m*lambda));
            
            %Compute the derivative term of dUdphi
            if m > 0
                Knm = sqrt((n-m)*(n+m+1));
            else
                Knm = sqrt(n*(n+1)/2);
            end
            
            %Impose that legendre polynomial with P(n,m+1) = 0 if m == n
            if m < n
                ind = find(leg_pol_data(:, 1) == n & leg_pol_data(:, 2) == m+1);
                leg_pol_plus1 = leg_pol_data(ind, 3);
            else
                leg_pol_plus1 = 0;
            end
            
            %Update last derivative value
            dUdphi = dUdphi + (R0/rnorm)^n * ( -m*tan(phi)*leg_pol + Knm*leg_pol_plus1 ) * (Cbar*cos(m*lambda) + Sbar*sin(m*lambda));
    end
end


    %Final moltiplication
    dUdr = -G*mass/rnorm^2*dUdr;
    dUdlam = G*mass/rnorm*dUdlam;
    dUdphi = G*mass/rnorm*dUdphi;
    
    %Building accelerations in carthesian frame
    acc_x = ( (1/rnorm * dUdr - y(3)/(rnorm^2*eta) * dUdphi)*y(1) - (1/eta^2 * dUdlam)*y(2) );
    acc_y = ( (1/rnorm * dUdr - y(3)/(rnorm^2*eta) * dUdphi)*y(2) + (1/eta^2 * dUdlam)*y(1) );
    acc_z = ( 1/rnorm * dUdr * y(3) + eta/rnorm^2 * dUdphi );
    
    acc = [acc_x; acc_y; acc_z];

end