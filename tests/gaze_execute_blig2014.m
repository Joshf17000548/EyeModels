function gaze_results = gaze_execute_blig2014(calib, data_right, data_left, screen)

targets = (screen.targets)';

pupil_right = data_right.features.PupilCentreImage;
pupil_left = data_left.features.PupilCentreImage;

glint_right = data_right.features.glintImage;
glint_left = data_left.features.glintImage;

scale = vecnorm(pupil_right - pupil_left,2,2);

vector_right = (pupil_right - glint_right)./scale;
vector_left = (pupil_left - glint_left)./scale;


constant = ones(size(pupil_right,1),1);
v = 1:135;


functions = [1, 2, 3, 4, 2, 5, 6, 7, 7, 7, 9, 9; % X functions
             1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 8, 9]; % Y functions
         
% for i = 1:135
%  figure(4)
% plot(glint_right(:,1), glint_right(:,2), '.','MarkerSize',15)
% set(gca, 'YDir','reverse')
% xlabel('x (pixels)')
% ylabel('y (pixels)')
% xlim([585 621])
% ylim([620 637])
% hold on

%  figure(5)
% plot(pupil_right(:,1), pupil_right(:,2), '.','MarkerSize',15)
% set(gca, 'YDir','reverse')
% xlabel('x (pixels)')
% ylabel('y (pixels)')
% xlim([570 635])
% ylim([608 645])
% hold on
% % 
% %  figure(3)
% % plot(vector(:,1), vector(:,2), '.','MarkerSize',8)
% % set(gca, 'YDir','reverse')
% % xlabel('x (pixels)')
% % ylabel('y (pixels)')
% % % print(filename_images +'\vector_1', '-dpng', '-r300')
% % hold on
% end

gaze_results = zeros(12,6);

for j = 1:size(calib, 1)
    
    for f = 1:size(functions, 2)

        ii_calib = calib(j,:) ~= 0;
        c = calib(j, ii_calib);
    
        dist = repmat(screen.position(3), size(targets,1),1);
       
        X_func = str2func(['X' int2str(functions(1, f))]);
        Y_func = str2func(['Y' int2str(functions(1, f))]);
               
        X_right = X_func(constant, vector_right(:, 1:2), v, c, targets(:, 1:2));
        Y_right = Y_func(constant, vector_right(:, 1:2), v, c, targets(:, 1:2));

        X_left = X_func(constant, vector_left(:, 1:2), v, c, targets(:, 1:2));
        Y_left = Y_func(constant, vector_left(:, 1:2), v, c, targets(:, 1:2));

        G = [mean([X_right, X_left], 2), mean([Y_right, Y_left], 2), dist];
        t = [targets(:, 1:2), dist];
        
     
        error = real(acosd(max(min(dot(G(v, :),t(v, :), 2)./(vecnorm(G(v, :),2,2).*vecnorm(t(v, :),2,2)),1),-1)));
%         error_grid = reshape(error, 12,6);
        gaze_results(f,j) = mean(error);
        
%         if((j == 3) && (f ==4))
%             
%             f = figure(4);
%             f.Position = [10 50 800 440];
%             x = linspace(1,15,15);
%             y = linspace(1,9,9);
%             [X,Y] = meshgrid(x,y);
%             bubblechart(X(:), Y(:), error);
%             bubblesize([1 15])
%             blgd = bubblelegend(['Gaze estimation' char(10) 'error (' char(176) ')'],'Location','eastoutside');
%             blgd.LimitLabels = {'0.28','0.01'};
%             xticks(1:15)
%             yticks(1:9)
%             yticklabels({'9','8','7','6','5','4','3','2','1'})
%             xlim([0, 16]);
%             ylim([0, 10]);
%             xlabel('X target')
%             ylabel('Y target')
%             set(gca, 'YDir','reverse')
%         end
%         
%         if((j == 5) && (f ==8))
%             figure(5)
%             
%             x = linspace(1,15,15);
%             y = linspace(1,9,9);
%             [X,Y] = meshgrid(x,y);
%             bubblechart(X(:), Y(:) ,error);
%             bubblesize([1 20])
%             t=1;
%         end
        
    end

end

gaze_results = gaze_results(:);
