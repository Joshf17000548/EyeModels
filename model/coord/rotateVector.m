function rotatedVectror = rotateVector(vector,rot, directionFlag)

if nargin==2
    directionFlag = 'forward';
    R = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);
end

% if nargin==4
%     R = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);
% end
% 
% if isempty(R)
%     R = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);
% end


%% Define the eye rotation matrix
% Assemble a rotation matrix from the head-fixed Euler angle rotations. In
% the head-centered world coordinate frame, positive azimuth, elevation and
% torsion values correspond to rightward, upward and clockwise (as seen
% from the perspective of the subject) eye movements

    R.azi = [cosd(rot(3)) -sind(rot(3)) 0; sind(rot(3)) cosd(rot(3)) 0; 0 0 1];
    R.ele = [cosd(-rot(2)) 0 sind(-rot(2)); 0 1 0; -sind(-rot(2)) 0 cosd(-rot(2))];
    R.tor = [1 0 0; 0 cosd(rot(1)) -sind(rot(1)); 0 sind(rot(1)) cosd(rot(1))];
    R.empty = false;

% We shift the points to each rotation center, rotate, shift back, and
% repeat. We must perform the rotation independently for each Euler angle
% to accomodate having rotation centers that differ by Euler angle.

rotateMatrix = R.tor * R.ele * R.azi;

switch directionFlag
    case 'forward'
        
        rotatedVectror = vector * rotateMatrix;
                
    case 'inverse'
        
        rotatedVectror = vector * rotateMatrix';
                        
end


end

