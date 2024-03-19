clc;
clear;
close all;

%% Parameters
eyeModels = {'reduced', 'asphericReduced', 'fourSurface', 'asphericFourSurface', 'finite', 'stochastic'};
eyes = 6:6;

for e = eyes
    
    if ~strcmp(eyeModels{e},'stochastic')
        participants = 1;
    else
        participants = 100;
    end
    
    for p = 93:participants
        
        for side = 1:2
            
            if side == 1.
                laterality = 'Right';
                [camera, screen, tracing, eyeparam] = parameters_right_blig2014();
            else
                laterality = 'Left';
                [camera, screen, tracing, eyeparam] = parameters_left_blig2014();
            end
            
            eye = modelEyeParameters('species',eyeModels{e}, 'alpha', [5.45, 2.5], 'spectralDomain', 'nir', 'participant', p, 'eyeLaterality', laterality);

            cameraCentre = camera.cameraPosition.translation;
            screenCentre = screen.position;
            screenTargets = screen.targets + screenCentre;
            
            features = [];
            features.target = [];
            features.eyeModel = [];
            features.pupilRad = [];
            features.glintImage = [];
            features.PupilCentreImage = [];
            features.headMovement = [];
            features.glintNumber = [];
            features.errors = [];
            
            for t = 1:size(screen.targets,2)
                
                [p, side, e, t]
                %% Build Environment
                
                camera.cameraPosition.translation = cameraCentre ;
                screen.position = screenCentre ;
                screen.targets = screenTargets ;
                
                sceneGeometry = createSceneGeometry(...
                    'eye',eye,...
                    'intrinsicCameraMatrix', camera.cameraIntrinsic.matrix, ...
                    'sensorResolution', camera.cameraIntrinsic.sensorResolution, ...
                    'radialDistortionVector', camera.cameraIntrinsic.radialDistortionVector, ...
                    'cameraTranslation', camera.cameraPosition.translation, ...
                    'cameraRotation', camera.cameraPosition.rotation,...
                    'cameraGlintSourceRelative', camera.camearPositionlight.glintSourceRelative, ...
                    'screenPosition', screen.position, ...
                    'screenDimensions', screen.dimensions, ...
                    'spectralDomain',eyeparam.index);
                
                %% Rotate Eye
                
                [eyepose, fixate_error, desiredFixationPoint] = calcFixationPose(eye, screen.targets(:,t), eye.pupil, true);
                
                %% Plot scene
                %                 plotScene(sceneGeometry, eye, camera, screen, tracing)
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
                
                error = [fixate_error, trace_error]
                
                %% Store results
                
                features = parseFeatures(features, pupilEllipseParams, imagePoints, pointLabels, screen.targets(:,t), eye.pupil, [0, 0, 0], error);
                
            end
            
   
            save(['tests/Blignaut2014/Data/None/' int2str(p) '_' eyeModels{e} '_' int2str(side) '.mat'], 'features')
        end
        
    end
end

