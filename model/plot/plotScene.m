function  plotScene(scene)


boundingBox = [-20 150 -50 50 -50 20];

%% Plot eye
plotOpticalSystem(scene,'newFigure',true);
hold on

%% Plot light sources

% lights = camera.camearPositionlight.glintSourceRelative;
% for i = 1:size(lights, 1)
%     pos = camera.cameraPosition.translation + lights(i,:);
%     S = quadric.translate(quadric.scale(quadric.unitSphere,[2 2 2]), pos);
%     quadric.plotImplicitSurface(S, boundingBox, 'k');
%
% end

%% Plot screen
t=scene.targets
addScreenIcon(scene, [0, 0], convertWorldToEyeCoord(scene.targets));

%% Plot camera
addCameraIcon(scene);

%% Plot targets

%% Plot rays

%% Labels
xlabel('x')
ylabel('y')
zlabel('z')

%%^Coordinates
% print(filename_images+ 'env', '-dpng', '-r300')
% a1 = gca;
% f2 = figure;
% a2 = copyobj(a1,f2);
% hold on
% plot3([0 7], [0 0], [0 0], '-k', 'LineWidth',2)
% plot3([0 0], [0 7], [0 0], '-k','LineWidth',2)
% plot3([0 0], [0 0], [0 7], '-k','LineWidth',2)
% xlim([-25 7])
% ylim([-20 20])
% zlim([-20 20])
% print(filename_images+ 'eye', '-dpng', '-r300')

end

