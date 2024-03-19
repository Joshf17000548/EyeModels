function [camera, screen, tracing, eyeparam] = parameters_right_blig2014()

tracing.spectralDomain = 'nir'; % Wavelength of light used (nm)
% tracing.nPupil = 30; % Number of pupil rays for ray tracing
tracing.nPupil = 30; % Number of pupil rays for ray tracing

ipd = 63;

camera.camearPositionlight.glintSourceRelative = [0; ...
                                                  50 ;
                                                  40 ]; % mm, glint positions 

% camera.cameraIntrinsic.matrix = [7.8, 0, 0; 0 7.8 0; 0 0 1]; %[focal_length_x skew offset_x; 0 focal_length_y offset_y; 0   0  1] pixels
% camera.cameraIntrinsic.matrix = [3, 0, 0; 0 3 0; 0 0 1];
% camera.cameraIntrinsic.matrix = [10 0 0; 0 10 0; 0 0 1];
camera.cameraIntrinsic.matrix = [3751.43 0 800; 0 3751.43 600; 0 0 1];

camera.cameraIntrinsic.sensorResolution = [1600 1200];
camera.cameraIntrinsic.radialDistortionVector = [0, 0];
% camera.cameraPosition.translation = [ipd/2; 0; 0];
camera.cameraPosition.translation = [ipd/2; -300; 520]; % mm, relative to the apex of the eye
camera.cameraPosition.unrotated_translation = [ipd/2; 0; 600]; % mm, relative to the apex of the eye
% camera.cameraPosition.translation = [ipd/2; 0; 200]; % mm, relative to the apex of the eye
camera.cameraPosition.rotation = [atand(300/520); 0; 0];
% camera.cameraPosition.rotation = [0; 0; 0];

screen.position = [ipd/2; -140; 800];
screen.dimensions = [495 280];

screen.horTargets = 15;
screen.verTargets = 9;
screen.nTargets = screen.horTargets*screen.verTargets;
[X,Y] = meshgrid(linspace(screen.dimensions(1), 0,screen.horTargets), linspace(screen.dimensions(2), 0, screen.verTargets));
screen.targets = [X(:)'-(screen.dimensions(1)/2); Y(:)'-(screen.dimensions(2)/2); zeros(1, screen.nTargets)]; %World

eyeparam.fovea = [5.45, 2.5];
eyeparam.index = 'nir';

% environment(eye, camera, screen, tracing)

end