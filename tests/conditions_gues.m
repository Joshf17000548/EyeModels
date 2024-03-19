function cond = conditions(it, target, participants)

% std = 1.282; % 80%
% std = 1.645; % 90%
% std = 1.96; % 95%
% std = 2.236 % 98%
% std = 2.576 % 99% 

cond.pupil = 2.*ones(size(target,2),1);
cond.head = zeros(size(target,2), 3);
cond.file = 'None';
cond.target = target;
cond.noise = [0, 0];
cond.occlusions = false;

    switch it
        case 1 % None S1
            
        case 2 % Image noise S4
           cond.file = 'Left';
         
        case 3 % Pupil S2
            
            cond.file = 'Right';
            
        case 4 % Head

            cond.file = 'Forward';  
            
        case 5 % Head
            cond.file = 'Backward';  
    end
end

function head = gaze_head(centre, target, scale)

target_distance = target(1, 3)+90;
optic_axis = [centre(:,1), centre(:,2), target_distance*ones(size(target, 1), 1)]-centre;
target_vector = (target + repmat([0, 0, 90], size(target,1),1) - centre);
r = 77.8;
phi = asind(31.5/r);
r1 = r*cosd(phi);

% scale = abs(randn(1,1))
% scale = unifrnd(0,1, [size(target, 1),1]);


deg_matrix = zeros(size(target, 1),3);
for i = 1:size(target, 1)
    R = vrrotvec(target_vector(i,:), optic_axis(i,:));
    anglet = R(4);
    axist = R(1:3);
    deg = axist.*anglet *180/pi;
    deg_matrix(i,:) = deg.*scale(i);
end
deg_matrix = -1.*deg_matrix;

y_x = r1.*cosd(deg_matrix(:,1))-r1;
y_z = -1.*r1*sind(deg_matrix(:,1));
z_x = r.*(cosd(-phi) - cosd(-phi-deg_matrix(:,2)));
z_y = -1*(r.*sind(phi-deg_matrix(:,2))-r*sind(phi));

head = [z_y, y_z, y_x + z_x];

end

function tar = fixation_error(target, std, angle)


%     scale = normrnd(0, 1/std, [size(target,1),2]);
    tar = zeros(size(target,1),3);
    for i = 1:size(target, 1)
        
        t = [0,0,0; target(i,3), target(i,2), target(i,1)]';
        o = [0,0,0; target(i,3),0,0]';
        [~, angle_xy, angle_xz] = quadric.angleRays( t, o);
        
        X(1) = target(i,3)*tand(angle_xz + unifrnd(-0.392, 0.392, [1,1]));
        X(2) = -1.*target(i,3)*tand(angle_xy + unifrnd(-0.392, 0.392, [1,1]));
        
        tar(i,:) = [X(1), X(2), target(i,3)];
%        X
%         R = [0, 0, 0; new_u(2), new_u(1), -new_u(3)]' 
%         X = quadric.intersectRayPlane([0; 0; 1],target(1,:)',R)
    end
end

function [pupil_noise] = feature_extraction(target)


k = unifrnd(0,1, [size(target, 1),1]);
k((k>0.16)) = 2.3260;
k((k>0.12) & (k<=0.16)) = 1.645;
k(k<=0.12) = 0.842;

optic_axis = [zeros(size(target, 1), 1), zeros(size(target, 1), 1), target(:, 3)];

     deg_scale = zeros(size(target, 1), 2);
     pupil_noise = zeros(size(target, 1), 2);
    for j = 1:size(target, 1)
        
        R = vrrotvec(target(j,:), optic_axis(j,:));
        anglet = R(4);
        axist = R(1:3);
        deg = axist.*anglet *180/pi;
        
        deg_scale(j,1) = normrnd(0, -1*deg(1)/18.2, [1, 1]);
        if deg_scale(j,1) < 0
            deg_scale(j,1) = 0;
        end
        pupil_noise(j,:) = normrnd(0, deg_scale(j,1)*1./k(j,1), [1, 2]);
%         pupil_fixate(j,1) = norm(pupil_new - pupil(j,:));
        

    end
end



