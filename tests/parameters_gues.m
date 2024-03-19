function [camera, screen, tracing, eyeparam] = parameters_gues()

tracing.spectralDomain = 'nir'; % Wavelength of light used (nm)
% tracing.nPupil = 30; % Number of pupil rays for ray tracing
tracing.nPupil = 30; % Number of pupil rays for ray tracing

ipd = 63;

camera.camearPositionlight.glintSourceRelative = [188.5, -188.5; ...
                                                  150.5, 150.5 ;
                                                  35, 35 ]; % mm, glint positions 

% camera.cameraIntrinsic.matrix = [7.8, 0, 0; 0 7.8 0; 0 0 1]; %[focal_length_x skew offset_x; 0 focal_length_y offset_y; 0   0  1] pixels
% camera.cameraIntrinsic.matrix = [3, 0, 0; 0 3 0; 0 0 1];
% camera.cameraIntrinsic.matrix = [10 0 0; 0 10 0; 0 0 1];
camera.cameraIntrinsic.matrix = [4666.67 0 320; 0 4666.67 240; 0 0 1];
camera.cameraIntrinsic.matrixmm = [35 0 2.4; 0 35 1.8; 0 0 1];

camera.cameraIntrinsic.sensorResolution = [640 480];
camera.cameraIntrinsic.radialDistortionVector = [0, 0];
% camera.cameraPosition.translation = [ipd/2; 0; 0];
camera.cameraPosition.translation = [ipd/2; -150.5; 615];
% camera.cameraPosition.translation = [ipd/2; -250; normrnd(619.45,39.6, [1,1])]; % mm, relative to the apex of the eye
camera.cameraPosition.unrotated_translation = [ipd/2; -150.5; 650]; % mm, relative to the apex of the eye
% camera.cameraPosition.translation = [ipd/2; 0; 200]; % mm, relative to the apex of the eye
camera.cameraPosition.rotation = [atand(150.5/615); 0; 0];
% camera.cameraPosition.rotation = [0; 0; 0];

screen.position = [ipd/2; 0; 650];
screen.dimensions = [377 301];

screen.horTargets = 3;
screen.verTargets = 3;
screen.nTargets = screen.horTargets*screen.verTargets;
[X,Y] = meshgrid(linspace(screen.dimensions(1), 0, screen.horTargets), linspace(screen.dimensions(2), 0, screen.verTargets));
screen.targets = [X(:)'-(screen.dimensions(1)/2); Y(:)'-(screen.dimensions(2)/2); zeros(1, screen.nTargets)]; %World

eyeparam.fovea = [5.45, 2.5];
eyeparam.index = 'nir';

% environment(eye, camera, screen, tracing)

end