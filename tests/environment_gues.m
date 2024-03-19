clc;
clear;
close all;

%% Parameters
eyeModels = {'reduced', 'asphericReduced', 'fourSurface', 'asphericFourSurface', 'finite', 'stochastic'};
eyes = 6:6;
% eyes = [1, 5];
positions = [0  , 0, 0;
    0  , 30, 0;
    0  , -30, 0;
    40 , 0, 0;
    -40, 0, 0];

[camera, screen, tracing, eyeparam] = parameters_gues();

cameraCentre = camera.cameraPosition.translation;
screenCentre = screen.position;
screenTargets = screen.targets + screenCentre;

for e = eyes

    if ~strcmp(eyeModels{e}, 'stochastic')
        participants = 1;
    else
        participants = 100;
    end

    for p = 51:participants

        eye = modelEyeParameters('species',eyeModels{e}, 'alpha', eyeparam.fovea, 'spectralDomain', eyeparam.index, 'participant', p);

        features = [];
        features.target = [];
        features.eyeModel = [];
        features.pupilRad = [];
        features.glintImage = [];
        features.PupilCentreImage = [];
        features.headMovement = [];
        features.glintNumber = [];
        features.errors = [];

        for pos = 1:size(positions,1)

            for t = 1:size(screen.targets,2)

                [p, e, pos, t]
                %% Build Environment

                camera.cameraPosition.translation = cameraCentre - positions(pos,:)';
                screen.position = screenCentre - positions(pos,:)';
                screen.targets = screenTargets - positions(pos,:)';

                sceneGeometry = createSceneGeometry(...
                    'eye',eye,...
                    'intrinsicCameraMatrix', camera.cameraIntrinsic.matrix, ...
                    'intrinsicCameraMatrixmm', camera.cameraIntrinsic.matrixmm, ...
                    'sensorResolution', camera.cameraIntrinsic.sensorResolution, ...
                    'radialDistortionVector', camera.cameraIntrinsic.radialDistortionVector, ...
                    'cameraTranslation', camera.cameraPosition.translation, ...
                    'cameraRotation', camera.cameraPosition.rotation,...
                    'cameraGlintSourceRelative', camera.camearPositionlight.glintSourceRelative, ...
                    'screenPosition', screen.position, ...
                    'screenDimensions', screen.dimensions, ...
                    'spectralDomain',eyeparam.index, ...
                    'targets',screen.targets);

                %% Rotate Eye
                [eyepose, fixate_error, desiredFixationPoint] = calcFixationPose(eye, screen.targets(:,t), eye.pupil, true);
                
                %% Plot scene
%                                     plotScene(sceneGeometry)
                %                 renderEyePose(eyepose, sceneGeometry, 'modelEyeSymbolSizeScaler', 0.4, ...
                %                 'nStopPerimPoints',tracing.nPupil, ...
                %                 'addPseudoTorsion',true)

                %% Execute ray tracing

                [pupilEllipseParams, ~, imagePoints, worldPoints, headPoints, eyePoints, cameraPoints, pointLabels, ~, ~, paths, glintPath, trace_error] = ...
                    projectModelEye(eyepose, sceneGeometry, ...
                    'cameraTrans',[0;0;0], ...
                    'fullEyeModelFlag', false, ...
                    'replaceReflectedPoints', true, ...
                    'nStopPerimPoints',tracing.nPupil, ...
                    'addPseudoTorsion',true);
                
                error = [fixate_error, trace_error];
                
                     %% Store results
                features = parseFeatures(features, pupilEllipseParams, imagePoints, pointLabels, convertWorldToEyeCoord(screen.targets(:,t)), eye.pupil, positions(pos,:), error);

            end

            if   (sum(positions(pos,:)) == 0)
                save(['tests/Guestrin/Data/None/' int2str(p) '_'  'scene.mat'], 'sceneGeometry')
            end

        end
        
        if strcmp(eyeModels(e), 'stochastic')
            save(['tests/Guestrin/Data/None/' int2str(p) '_'  'eye.mat'], 'eye')
        end

        save(['tests/Guestrin/Data/None/' int2str(p) '_' eyeModels{e} '.mat'], 'features')

    end
end

