function [mappingScore, exploitingScore, new_scores, new_known_map, penalty_of_t] = scientificScore(y, t, spacecraft_data)

data = spacecraft_data.data_guidance;
r_impact = data.r_impact;
r_escape = data.r_escape;
fov1 = data.scientificFov1;
fov2 = data.scientificFov2;

F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;

mapping = spacecraft_data.data_asteroids.mapping;
scores  = spacecraft_data.data_asteroids.features.score;
known_map_features = spacecraft_data.data_asteroids.features.known_map_features;
known_map = spacecraft_data.data_asteroids.mapping.known_map;

new_known_map = known_map;
new_scores = scores;

single_score_mapping = 1 / length(known_map);

r = y(:, 1:3);
n_time = size(y, 1);

mappingScore    = zeros(n_time, 1);
exploitingScore = zeros(n_time, 1);
r_sun_vec       = zeros(n_time, 3);
R_body_cam_all  = zeros(3, 3, n_time);

all_types   = {new_scores.type};
isEmission  = strcmp(all_types, 'emission')';
isPhase     = strcmp(all_types, 'relative')';

pE = find(isEmission, 1, 'first');
pP = find(isPhase,   1, 'first');
emissionTasks = new_scores(pE).angles;
phaseTasks    = new_scores(pP).angles;
idealRange    = new_scores(pP).ideal_range;

C = reshape([new_scores.completeness], 4, []).';

emissionScore = ones(size(emissionTasks,1),1) * (new_scores(pE).score / size(emissionTasks,1));
phaseScore    = ones(size(phaseTasks,1),1)    * (new_scores(pP).score / size(phaseTasks,1));

for i = 1:n_time
    r_sun = cspice_spkgeo(2000433, t(i), 'ECLIPJ2000', 10);
    r_sun = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t(i)) * (-r_sun(1:3));
    r_sun_vec(i, :) = r_sun;
    R_body_cam_all(:, :, i) = body2camera(y(i, 1:3), y(i, 4:6));
end
r_sun_vec = r_sun_vec ./ vecnorm(r_sun_vec, 2, 2);

centroid_vec = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3;

space  = 2;
sigma2 = 0.5;

for i = 1:n_time
    r_now = r(i, :);
    R_cam = R_body_cam_all(:, :, i);
    r_sun = r_sun_vec(i, :);

    r_relative = r_now - centroid_vec;                 
    norm_rrel  = vecnorm(r_relative, 2, 2);            

    emissionAngles = acos( dot(r_relative, N', 2) ./ norm_rrel ); 
    incidenceAngles = acos((r_sun*N)');                           
    phaseAngles     = acos( (r_relative*r_sun') ./ norm_rrel );   

    faceInFov = check_FOV(V, F, r_now, fov1, fov2, R_cam); 

    checkEmissionTasks = emissionAngles > emissionTasks(:,1).' & emissionAngles < emissionTasks(:,2).';
    checkPhaseTasks    = phaseAngles    > phaseTasks(:,1).'    & phaseAngles    < phaseTasks(:,2).';

    checkIncidence = incidenceAngles < mapping.incidence;
    checkEmission  = emissionAngles  < mapping.emission;

    delta        = abs(norm_rrel - idealRange);
    rangeScaling = exp( - max(delta - space, 0).^2 / (2*sigma2) );

    RSmat = rangeScaling(:);             
    deltaC = max(RSmat - C, 0);       

    
    maskE = checkEmissionTasks & isEmission;                          
    maskP = checkPhaseTasks    & isPhase   & checkEmission;   

    contrE = deltaC .* maskE  * emissionScore;  
    contrP = deltaC .* maskP  * phaseScore; 

    gate = known_map_features .* faceInFov .* checkIncidence; 
    exploitingScore(i) = sum( gate .* (contrE + contrP) );

    checkRange = norm_rrel < 50;
    mappable   = checkRange & checkIncidence & checkEmission & faceInFov;
    newlyMapped = (~new_known_map) & mappable;
    mappingScore(i) = sum(newlyMapped) * single_score_mapping;
    new_known_map   = new_known_map | newlyMapped;

    activeMask = (maskE | maskP) & (faceInFov & checkIncidence);  
    newC = max(C, RSmat); 
    C(activeMask) = newC(activeMask);

end

for k = 1:size(C,1)
    new_scores(k).completeness = C(k,:).';
end

allRNorm = vecnorm(y(:, 1:3), 2, 2);
penalty_of_t = ( allRNorm < r_impact | allRNorm > r_escape ) * (-5000);
end
