clc;
clear;
close all;

eyeModels = {'reduced', 'asphericReduced', 'fourSurface', 'asphericFourSurface', 'finite', 'stochastic'};
% eyeModels = {'reduced'};
% cond = [1, 2, 3, 4, 5];
n_participants = 100;
eyes = 1:6;
% eyes = 1;
% positions = 1:5;

[cam, screen, ~, ~] = parameters_gues();

calib = linspace(1, 9, 9);

gaze_res_total = zeros(size(eyeModels,2), 2);
gaze_res_std = zeros(1, 2);


data = load(['tests/Guestrin/Data/None/1_scene.mat']);
scene = data.sceneGeometry;
% plotScene(scene)

targets = (screen.targets+screen.position);

gues_centre = atand([2.93, 4.88, 4.93, 4.68]./650);
gues_edge = atand([4.25, 5.33, 5.74, 8.12]./650);
gues_res_vec = [mean(gues_centre), mean(gues_edge)];
gues_std = [std(gues_centre), std(gues_edge)].*3.182;

biometry = zeros(11, n_participants);


f = figure(1);
f.Position = [10 50 800 440];
hold on
x = 1:2;
xconf = [x x(end:-1:1)] ;
yconf = [gues_res_vec+gues_std gues_res_vec(end:-1:1)-gues_std(end:-1:1)];
p = fill(xconf,yconf,'blue');
p.FaceColor = [0.84 0.93 0.96];
p.EdgeColor = 'none';
plot(x, gues_res_vec, '--o', 'LineWidth', 1)

f = figure(2);
f.Position = [10 50 800 440];
hold on
x = 1:2;
xconf = [x x(end:-1:1)] ;
yconf = [gues_res_vec+gues_std gues_res_vec(end:-1:1)-gues_std(end:-1:1)];
p = fill(xconf,yconf,'blue');
p.FaceColor = [0.84 0.93 0.96];
p.EdgeColor = 'none';
plot(x, gues_res_vec, '--o', 'LineWidth', 1)

for i = eyes
    
    if ~strcmp(eyeModels{i}, 'stochastic')
        participants = 1;
    else
        participants = n_participants;
    end
    
    gaze = zeros(n_participants,2);
    gaze_res_ii = zeros(2, participants);
    
    gaze_results = zeros(participants,2);
    if ~strcmp(eyeModels{i}, 'stochastic')
        for p = 1:participants
            
            data = load(['Guestrin/Data/None/' int2str(p) '_' eyeModels{i} '.mat' ]);
            
            if strcmp(eyeModels{i}, 'stochastic')
                
                d = load(['Eyes/' int2str(p) '_eye_1.mat']);
                eye.R = d.eye.R;
                
                Rca = 337.2./eye.R(:,2);
                Qca = -1*(eye.R(:,5).^2);
                Rcp = (1337.2-1378)./eye.R(:,6);
                Qcp = -1*(eye.R(:,9).^2);
                CCT = eye.R(:,13);
                ACD = eye.R(:,14);
                LT = eye.R(:,12);
                AL = eye.R(:,15);
                Rla = eye.R(:,10);
                Rlp = -1*eye.R(:,11);
                P = eye.R(:,16);
                
                biometry(:,p) = [Rca; Qca; Rcp; Qcp; Rla; Rlp; CCT; ACD; LT; AL; P];
            end
            
            [i, p];
            
%             gaze_results(p,:) = gaze_execute_gues(data);
%             
            Z_res = zscore(gaze_results);
            ii_outliers = abs(Z_res) > 5;
            %         gaze_results(ii_outliers) = nan;
            sum(ii_outliers);
            
        end
    else
        d = load('tests/Guestrin/gaze data.mat');
        gaze_results = d.gaze_results;
        biometry = d.biometry;
    end

%     save(['tests/Guestrin/gaze data.mat'], 'gaze_results', 'biometry')
    
%     gaze_res_total(i,:) = mean(gaze_results,1);
    d = load('tests/Guestrin/norm data.mat');
    gaze_res_total = d.gaze_res_total;
    
    if ~strcmp(eyeModels{i}, 'stochastic')
        figure(1)
        plot(1:2, gaze_res_total(i, :), '-o', 'LineWidth', 1)
    else
        gaze_res_std = std(gaze_results,1);
    end
    
end

figure(1)
% lgd = legend('95% CI', 'Experimental', 'Reduced', 'Aspheric reduced', 'Four surface', 'Aspheric four surface', 'Finite');
lgd = legend('Experimental 95% CI', 'Experimental mean', 'Reduced (E_1)', 'Aspheric reduced (E_2)', 'Four surface (E_3)', 'Aspheric four surface (E_4)','Location','northwest');
% lgd = legend('Experimental 95% CI', 'Experimental mean', 'Reduced', 'Aspheric reduced', 'Four surface', 'Aspheric four surface');
% title(lgd,'Legend')
xlabel('Head positions')
ylabel(['Gaze estimation error (' char(176) ')'])
xticks([1, 2])
xticklabels({'H1','H2'})
xlim([0.8, 2.2]);
ylim([0, 1.1]);
 print(['C:\Users\LIQID Medical\OneDrive\Thesis\LaTex\figs\complexity_results' '\gues.png'], '-dpng', '-r300')

f = figure(2);
crit = 1.984;
f.Position = [810 50 800 440];
hold on
xconf = [x x(end:-1:1)] ;
avg = gaze_res_total(6,:);
std_dev = gaze_res_std.*1.984;
yconf = [avg+std_dev avg(end:-1:1)-std_dev(end:-1:1)];
p = fill(xconf,yconf,'blue','FaceAlpha',0.5);
p.FaceColor = [1 0.8 0.66];
p.EdgeColor = 'none';
plot(1:2, avg, '--o', 'LineWidth', 1)
lgd = legend('Experimental 95% CI', 'Experimental mean', 'Simulated 95% CI', 'Simulated mean' ,'Location','northwest');
% title(lgd,'Legend')
xlim([0.8, 2.2]);
ylim([0, 1.1]);
xlabel('Head positions')
ylabel(['Gaze estimation error (' char(176) ')'])
xticks([1, 2])
xticklabels({'H1','H2'})


%%% G1
fprintf('G1\nxxxxxxxxxxxxxxxxxx\n')
for i = eyes
    fprintf('\n%s\n', eyeModels{i})
    if strcmp(eyeModels{i}, 'stochastic')
        
        cen = gaze_results(:,1);
        per = gaze_results(:,2);
        crit = 2.6;
        [~,p,~,stats] = ttest2(gues_centre,cen, 'Vartype','unequal');
        fprintf('Centre: t(%d) = %.2f, p = %.2f, h = %d\n', stats.df, stats.tstat, p, abs(stats.tstat)>crit);
        crit = 2.9;
        [~,p,~,stats] = ttest2(gues_edge,per, 'Vartype','unequal');
        fprintf('Peripheral: t(%d) = %.2f, p = %.2f, h = %d\n', stats.df, stats.tstat, p, abs(stats.tstat)>crit)
    else
        crit = 3.182;
        [~,p,~,stats] = ttest(gues_centre, gaze_res_total(i,1));
        fprintf('Centre: t(%d) = %.2f, p = %.2f, h = %d\n', stats.df, stats.tstat, p, abs(stats.tstat)>crit);
        [~,p,~,stats] = ttest(gues_edge, gaze_res_total(i,2));
        fprintf('Peripheral: t(%d) = %.2f, p = %.2f, h = %d\n', stats.df, stats.tstat, p, abs(stats.tstat)>crit)
    end
    
end
fprintf('\nxxxxxxxxxxxxxxxxxx\n')


fprintf('\nG2\nxxxxxxxxxxxxxxxxxx\n\n')
if eyes(end) == 6
    crit = 8.554;
    [~,p,~,stats] = vartest2(cen(:), gues_centre);
    fprintf('Centre: f(%d, %d) = %.2f, p = %.2f, h = %d\n', stats.df1 , stats.df2, stats.fstat, p, stats.fstat>crit);
    crit = 8.554;
    [~,p,~,stats] = vartest2(per(:), gues_edge);
    fprintf('Peripheral: f(%d, %d) = %.2f, p = %.2f, h = %d\n', stats.df1, stats.df2, stats.fstat, p, stats.fstat>crit)
end
fprintf('\nxxxxxxxxxxxxxxxxxx\n\n')

fprintf('\nG3\nxxxxxxxxxxxxxxxxxx\n')
diff = [gues_res_vec(2)-gues_res_vec(1); gaze_res_total(:,2)- gaze_res_total(:,1)];
fprintf('Experimental = %.2f\n', diff(1));
for i = eyes
    fprintf('%s = %.2f\n', eyeModels{i}, diff(i+1));
end
fprintf('\nxxxxxxxxxxxxxxxxxx\n')

fprintf('\nBiometry\nxxxxxxxxxxxx\n')

eye_stats = [mean(biometry,2), std(biometry,[],2)];

% res_vec = gaze_results(:,:)';
% res_z = zscore(gaze_results, [], 1);
% res_vec_z = res_z';


% bio_z = biometry;
bio_z = zscore(biometry,[],2);
% bio_mat = repelem(bio_z', size(res_vec,1), 1);


r_result = zeros(4, 11);
for i = 1:2
    
%     figure(7)
%     hist(gaze_results(:,i))
%     figure(6)
    skewness(gaze_results(:,i));
    median(gaze_results(:,i));
%     res_vec = gaze_results(:,i)';
%     bio_mat = repelem(bio_z', size(res_vec,1), 1);
%     bio_result = res_vec'\bio_z'
    for j = 1:11
        [R,P] = corrcoef(gaze_results(:,i), bio_z(j,:)) ;
%         r_result((i-1)*2+1, j) = R(1,2);
%         r_result((i-1)*2+2, j) = P(1,2);
        slope = bio_z(j,:)'\gaze_results(:,i);
        fprintf('    &    %.2f (p = %.2f) (s = %.2f)', R(1,2), P(1,2), slope)
    end
    fprintf('\\\\ \n')

end

f = figure(6);
f.Position = [10 50 1200 440];
hold on
k = 1;

% p = anovan(gaze_results(:,1),{bio_z(8,:)', bio_z(9,:)', bio_z(10,:)'}, 'model','interaction')
params = {'', 'Q_a_c', '', '', 'R_a_l', '', '', 'ACD', 'LT', 'AL', ''};
% for j = [2]
% 
%     for i = 1:2
%          [R,P] = corrcoef(gaze_results(:,i), bio_z(j,:)) ;
%         ax = subplot(1,2,k);
%         scatter(bio_z(j,:), gaze_results(:,i), '.')
%         h1 = lsline(ax);
%         h1.LineWidth = 2;
% %         slope = bio_z(j,:)'\gaze_results(:,i);
%         mdl = fitlm(bio_z(2,:)', gaze_results(:,1));
%         
%         title([params{j}],'FontSize',16)
%         text(-3.5,0.2,['Error =  ' num2str(mdl.Coefficients.Estimate(1), '%.2f') ' - ' num2str(abs(mdl.Coefficients.Estimate(2)), '%.2f') 'SD'],'FontSize',16)
%         text(-3.5,0.1,['R^2 =  ' num2str(mdl.Rsquared.ordinary, '%.2f')],'FontSize',16)
%         xlabel('SD','FontSize',16)
%         ylabel(['Error (' char(176) ')'] ,'FontSize',16)
%         k = k+1;
%     end
% end

ca = {'H1', 'H2'};
params = {'', 'Q_a_c', '', '', 'R_a_l', '', '', 'ACD', 'LT', 'AL', ''};
bes = {[2], [2], [2], [2], [2], [2], [2]}
    for i = 1:2
%          [R,P] = corrcoef(gaze_results(:,i), bio_z(j,:)) ;
        ax = subplot(1,2,k);
%         scatter(bio_z(j,:), min_val_p(:,i), '.')
%         h1 = lsline(ax);
%         h1.LineWidth = 1;
%         slope = bio_z(j,:)'\gaze_results(:,i);
        mdl = fitlm(bio_z(bes{i},:)', gaze_results(:,i));
        h = plot(mdl,'color', 'b');
        hLeg = legend('example');
        set(hLeg,'visible','off')
        set(h(1), 'Color', 'k')
        set(h(2), 'Color', 'b')
        set(h(3), 'visible','off')
        if length(h)>3
            set(h(4), 'visible','off')
        end
        
        title([ca{i}],'FontSize',14)
        str = num2str(mdl.Coefficients.Estimate(1), 2);
        for q =  1:(length(bes{i}))
            b = bes{i}
            if mdl.Coefficients.Estimate(q+1) > 0
                symb = ' + ' ;
            else
                symb = ' - ' ;
            end
            str = [str symb num2str(abs(mdl.Coefficients.Estimate(q+1)), 2 )  params{b(q)}]
        end
        text(-1.2,0.15,['Error =  ' str],'FontSize',14)
        text(-1.2,0.1,['R^2 =  ' num2str(mdl.Rsquared.ordinary, '%.2f')],'FontSize',14)
        xlabel('Adjusted whole model (SD)','FontSize',14)
        ylabel(['Adjusted error (' char(176) ')'] ,'FontSize',14)
        ylim([0, 0.6])
        xlim([-1.5, 1.5])
        k = k+1;
    end
fprintf('xxxxxxxxxxxxxxxxxx\n')








