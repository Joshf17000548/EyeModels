function gaze_results = gaze_execute_gues(data)

[camera, screen, tracing, eyeparam] = parameters_gues();
cameraCentre = camera.cameraPosition.translation;
screenCentre = screen.position;

pupil = data.features.PupilCentreImage;
glint = data.features.glintImage;
targets = data.features.target;
% targets = targets(:,)
head = data.features.headMovement;

% % for t = 1:9
% figure(3)
% plot(glint(:,1), glint(:,2), '.','MarkerSize',13)
% % plot(glint(2:2:18,1), glint(2:2:18,2), 'x','MarkerSize',8)
% hold on
% set(gca, 'YDir','reverse')
% xlabel('x (pixels)')
% ylabel('y (pixels)')
% % xlim([70 120])
% % ylim([230, 265])
% % % print(filename_images +'\glint_1', '-dpng', '-r300')


%  figure(4)
% plot(pupil(1:9,1), pupil(1:9,2), '.','MarkerSize',13)
% set(gca, 'YDir','reverse')
% xlabel('x (pixels)')
% ylabel('y (pixels)')
% xlim([50 120])
% ylim([230, 270])
% lgd = legend('Reduced', 'Aspheric reduced', 'Four surface', 'Aspheric four surface', 'Finite','Location','northeast');
% title(lgd,'Eye Models')
% print(filename_images +'\pupil_1', '-dpng', '-r300')
% hold on
% 
% end

if sum(isnan(pupil(:,1))) >0
    causeException = MException('MATLAB:notEnoughInputs','Not enough input arguments.');
    throw(causeException);
end

constant = ones(size(data.features.PupilCentreImage,1),1);
cc = 1:9;

%  
%         index = ((pos-1)*pos+1):pos*9;
%         p = pupil(index, :);
%         index_led = ((pos-1)*pos*2+1):pos*2*9;
%         g = glint(index_led, :);
        
        camera_pos = repmat(cameraCentre', size(pupil, 1),1) - head;
        screen_pos = repmat(screenCentre', size(pupil, 1), 1) - head;
        led_pos = repmat(convertWorldToEyeCoord(camera.camearPositionlight.glintSourceRelative),size(pupil,1), 1);
        led_pos = led_pos + repelem(camera_pos(:, [3 1 2]),2,1);
        
        led_pos = convertEyeToCameraCoord(led_pos,repelem(camera_pos,2,1));
%         plot3(led_pos(:,1), led_pos(:,2), led_pos(:,3), 'ob');

        targets = convertEyeToCameraCoord(targets,repelem(camera_pos,1,1));
        
        glint = convertImageToCameraCoord(glint, atand(camera.cameraPosition.translation(2)/camera.cameraPosition.translation(3)), ...
            camera.cameraIntrinsic.sensorResolution, 7.5*10^-3);

%         t = convertCameraToEyeCoord(glint,repelem(camera_pos,2,1))
%         plot3(t(1,1), t(1,2), t(1,3), '.b');
        
        pupil = convertImageToCameraCoord(pupil, atand(camera.cameraPosition.translation(2)/camera.cameraPosition.translation(3)), ...
            camera.cameraIntrinsic.sensorResolution, 7.5*10^-3);
        

        gaze_results = Guestrin(pupil,glint,targets, camera_pos, screen_pos, led_pos);
%         gaze_results = zeros(45,1);
        gaze_results = [mean(gaze_results(1:9)), mean(gaze_results(10:end))];
        
        
%         error = real(acosd(max(min(dot(G(v, :),targets(v, :), 2)./(vecnorm(G(v, :),2,2).*vecnorm(targets(v, :),2,2)),1),-1)));
%         error_all = real(acosd(max(min(dot(G,targets, 2)./(vecnorm(G,2,2).*vecnorm(targets,2,2)),1),-1)));
        
% %         size(error_shape)
%         reshape(error_all, 9,15)
        
%         error_all(c,:)
%         mean(error)
%   
%         bubblechart(targets(:,1),targets(:,2),error_shape,'red');
%         bubblesize(([min(error_shape) max(error_shape)]+0.1).*8)
%         [f,j]
        
%         gaze_results(f,j) = mean(error);
        
        
%         figure(1)
%         plot(vector(:,2), targets(:,2), 'o')
%         hold on
%         plot(X, Y, '.')
%         hold off
        

