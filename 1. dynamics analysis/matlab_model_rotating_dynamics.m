function dy = matlab_model_rotating_dynamics( t, y, mass, omega, degree)

%Position and velocity
r = y(1:3);
v = y(4:6);

C = zeros(degree+1, degree+1);
S = zeros(degree+1, degree+1);

C(1,1)= 1;
S(1,1)= 0;
for i = 1:degree
    for j = 0:i
        C(i+1, j+1) = eros_armonic(i,j, 'c');
        S(i+1, j+1) = eros_armonic(i,j, 's');
    end
end

if degree == 2
    C(3,2) = 0;
    S = zeros(3,3);
end

Re = 16000;

GM = astroConstants(1)*10^(9)*mass;

save('myCustomGravityModel.mat', 'GM', 'Re', 'degree', 'C', 'S');

% Specifica il nome del tuo file dati
my_datafile = 'myCustomGravityModel.mat';

% Ottieni un handle alla tua funzione dfreader
my_dfreader = @myMatFileReader;

% Definisci l'azione per gli input fuori intervallo (opzionale)
action = 'Error'; % O 'Warning', 'None'

% Chiama la funzione gravitysphericalharmonic
[gx, gy, gz] = gravitysphericalharmonic(r'*1000, 'Custom', degree, {my_datafile, my_dfreader}, action);

gravity = [gx; gy; gz]/1000;
% Set the derivatives of the state
dy = [ v
 gravity - 2*cross(omega, v) - cross(omega, cross(omega, r))];

end