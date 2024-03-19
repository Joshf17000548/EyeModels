clc;
clear;
close all;

eyeModels = {'reduced', 'asphericReduced', 'fourSurface', 'asphericFourSurface', 'finite', 'stochastic'};
% eyeModels = {'reduced', 'asphericReduced', 'fourSurface', 'asphericFourSurface', 'finite'};
regression = 1:12;
columns = 1:6;
n_participants = 100;
eyes = 1:6;
% eyes = [1];

blig_results = [1.4, 1.29, 1.15, 1.16, 1.14, 1.11;
    1.37, 1.09, 0.96, 0.94, 0.94, 0.90;
    2.94, 0.91, 0.76, 0.75, 0.73, 0.68;
    4.93, 0.87, 0.72, 0.70, 0.68, 0.64;
    1.50, 1.09, 0.95, 0.96, 0.94, 0.90;
    1.72, 0.93, 0.84, 0.83, 0.81, 0.78;
    3.43, 0.80, 0.62, 0.62, 0.60, 0.56;
    6.17, 2.06, 0.58, 0.56, 0.53, 0.48;
    6.24, 2.09, 0.62, 0.59, 0.57, 0.53;
    7.76, 4.74, 0.70, 0.67, 0.62, 0.55;
    7.32, 3.07, 0.62, 0.60, 0.57, 0.52;
    8.74, 5.53, 0.70, 0.67, 0.62, 0.54];

calib = [1, 9, 68, 127, 135, zeros(1, 130);
    1, 5, 9, 64, 68, 72, 127, 131, 135, zeros(1, 126);
    3, 9, 19, 25, 39, 45, 64, 70, 93, 99, 109, 115, 129, 135, zeros(1, 121);
    1, 9, 14, 20, 26, 39, 43, 58, 64, 72, 78, 93, 97, 110, 116, 122, 127, 135, zeros(1, 117);
    1, 5, 9, 12, 16, 28, 32, 36, 48, 52, 64, 68, 72, 84, 88, 100, 104, 108, 120, 124, 127, 131, 135, zeros(1, 112);
    linspace(1, 135, 135)];


blig_sim = zeros(72, length(eyes));
blig = blig_results(regression, columns);
blig_res_vec = blig(:);
std(blig)
best_reg = zeros(length(eyes),6);
best_reg_error = zeros(length(eyes),6);

ex_best_error  = zeros(6, n_participants);

biometry = zeros(11, n_participants);
gaze_std = zeros(1, 6);

[~, screen, ~, ~] = parameters_right_blig2014();

actual_best_fit = [1.37, 0.87, 0.58, 0.56, 0.53, 0.48];
% actual_best_fit_dev = [0.25, 0.13, 0.12, 0.12, 0.12, 0.19, 0.11, 0.12].*2.11;

f = figure(1);
f.Position = [10 50 800 440];
hold on
x = 1:6;
plot(x, actual_best_fit, '--o', 'LineWidth', 1)

f = figure(2);
f.Position = [10 50 800 440];
hold on
x = 1:6;
plot(x, actual_best_fit, '--o', 'LineWidth', 1)

f = figure(3);
f.Position = [10 50 800 440];
hold on
x = 1:6;
plot(x, actual_best_fit, '--o', 'LineWidth', 1)

for i = eyes
    
    if ~strcmp(eyeModels{i}, 'stochastic')
        participants = 1;
    else
        participants = n_participants;
    end
    
    gaze_results = zeros(length(columns)*length(regression), participants);
    
    
    for p = 1:participants
%       for p = [1,4]
        
        [i, p];
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
        
        data_right = load(['Blignaut2014/Data/None/' int2str(p) '_' eyeModels{i} '_1' '.mat' ]);
        data_left = load(['Blignaut2014/Data//None/' int2str(p) '_' eyeModels{i} '_2' '.mat' ]);
        
        
        [row, column]=find(data_left.features.errors==1);
        
        gaze_results(:,p) = gaze_execute_blig2014(calib, data_right, data_left, screen);
        
    end
    
    blig_sim(:,i) = mean(gaze_results, 2);
    gaze_res_std = std(gaze_results, [], 2);
    
    res = reshape(blig_sim(:,i), [12, 6]);
    res_std = reshape(gaze_res_std, [12, 6]);
    
    [min_val, min_ii] = min(res, [], 1);
    best_reg(i, :) = min_ii;
    best_reg_error(i, :) = min_val;
    
    ii_min = sub2ind([12,6],[2,4,8,8,8,8],[1,2,3,4,5,6]);
% ii_min = sub2ind([12,6],[5,5,5,5,5,5],[1,2,3,4,5,6]);
    
    if ~strcmp(eyeModels{i}, 'stochastic')
        figure(1)
        plot(1:6, min_val, '-o', 'LineWidth', 1)
        figure(2)
        plot(1:6, blig_sim(ii_min,i), '-o', 'LineWidth', 1)
    else
        gaze_std = std(gaze_results,[], 2);
        best_std = gaze_std(ii_min);
        ex_best_error = gaze_results(ii_min,:);
    end

end
% repmat(actual_best_fit, 5, 1) - best_reg_error

figure(1)
lgd = legend('Experimental', 'Reduced (E_1)', 'Aspheric reduced (E_2)', 'Four-surface (E_3)', 'Aspheric four-surface (E_4)', 'Finite (E_5)');
% title(lgd,'Data')
xlabel('Calibration configuration')
ylabel(['Gaze estimation error (' char(176) ')'])
xticks([1, 2, 3, 4, 5, 6, 7])
xticklabels({'C5','C9','C14','C18','C23','C135'})
xlim([0, 7]);
ylim([0, 1.5]);
% print(['C:\Users\Joshua\OneDrive\Thesis\LaTex\figs\complexity_results' '\blig14.png'], '-dpng', '-r300')

figure(2)
lgd = legend('Experimental', 'Reduced (E_1)', 'Aspheric reduced (E_2)', 'Four-surface (E_3)', 'Aspheric four-surface (E_4)', 'Finite (E_5)');
% title(lgd,'Data')
xlabel('Calibration configuration')
ylabel(['Gaze estimation error (' char(176) ')'])
xticks([1, 2, 3, 4, 5, 6, 7])
xticklabels({'C5','C9','C14','C18','C23','C135'})
xlim([0, 7]);
ylim([0, 1.5]);
print(['C:\Users\LIQID Medical\OneDrive\Thesis\LaTex\figs\complexity_results' '\blig14.png'], '-dpng', '-r300')

f = figure(3);
xconf = [x x(end:-1:1)] ;
% avg = best_reg_error(6,[1, 2, 3, 4, 5, 6]);
% std_dev = best_std([1, 2, 3, 4, 5, 6])'.*1.984;
avg = mean(ex_best_error,2)';
std_dev = std(ex_best_error, [], 2)'.*1.984;

yconf = [avg+std_dev avg(end:-1:1)-std_dev(end:-1:1)];
p = fill(xconf,yconf, 'red','FaceAlpha',0.5);
p.FaceColor = [1 0.8 0.66];
p.EdgeColor = 'none';
plot(1:6, avg, '-o', 'LineWidth', 1)
lgd = legend('Experimental mean', 'Simulated 95% CI', 'Simulated mean');
xticks([1, 2, 3, 4, 5, 6, 7])
xticklabels({'C5','C9','C14','C18','C23','C135'})
xlim([0, 7]);
ylim([0, 1.5]);
xlabel('Calibration configuration')
ylabel(['Gaze estimation error (' char(176) ')'])
print(['C:\Users\LIQID Medical\OneDrive\Thesis\LaTex\figs\stochastic_results' '\blig14.png'], '-dpng', '-r900')

fprintf('Bl1\nxxxxxxxxxxxxxxxxxx\n')
diff_err = zeros(72, length(eyes)); 
diff_err_eye = zeros(72, length(eyes)); 
diff_sum = zeros(1, length(eyes)); 
blig_sim = [blig_sim zeros(72,1)];
 for i = eyes
     
     p = [1,1,1,1,1,1,1,1];
     crit = 1.984;
     if strcmp(eyeModels{i}, 'stochastic')
          for c = 1:72
             ex_avg = blig_res_vec(c);
             
             N_sim = 100;
             sim_avg = blig_sim(c,6);
             sim_std = gaze_std(c);

             tval = (sim_avg - ex_avg) / (sim_std/sqrt(N_sim)); 

            fprintf('%s, c%d: t(%d) = %.2f, p = %.2f, h = %d\n', eyeModels{i}, c, N_sim-1, tval, 1, abs(tval)>crit);
         end   
         
     else
         diff_err(:,i) = (blig_res_vec - blig_sim(:, i))./blig_res_vec*100 ;
         diff_err_eye(:,i) = (blig_sim(:, i+1) - blig_sim(:, i))./blig_sim(:, i)*100 ;
%          diff_sum(i) = length(find(abs(diff_err(:,i))<0.1));
     end
 end

% fprintf('%.2f\n', diff_err_eye)
mean(diff_err)
fprintf('\nxxxxxxxxxxxxxxxxxx\n')

fprintf('Bl2\nxxxxxxxxxxxxxxxxxx\n')
best_reg
fprintf('\nxxxxxxxxxxxxxxxxxx\n')

fprintf('Bl3\nxxxxxxxxxxxxxxxxxx\n')
increase = zeros(length(eyes),3);
 for i = eyes
     
     ii_min = sub2ind([12,6],[8,8,8,8,8,8],[1,2,3,4,5,6]);
     err = blig_sim(ii_min,i);
     increase(i,:) = [err(4) - err(3), err(5) - err(3), err(6) - err(3)];
 end
 increase
fprintf('\nxxxxxxxxxxxxxxxxxx\n')

fprintf('\nBiometry\nxxxxxxxxxxxx\n')
%%%% Biometry
% gaze_results = gaze_results(:, [1,2,3,5,6,7,8]);

eye_stats = [mean(biometry,2), std(biometry,[],2)]

% res_vec = gaze_results;
bio_z = zscore(biometry,[],2);

% % res_z = zscore(gaze_results, [], 2)
% for i = 1:6
%     k = 1;
%     res_vec = ex_best_error(i,:)';
%     figure(4)
%     hist(res_vec)
%     for j = 1:11
%                 r_result((i-1)*2+1, j) = R(1,2);
%         r_result((i-1)*2+2, j) = P(1,2);
%         
%         [R,P] = corrcoef(res_vec, bio_z(j,:)) ;
%         ax = subplot(6,2,k);
%         scatter(bio_z(j,:)', res_vec, '.')
%         lsline(ax)
%         
%         slope = bio_z(j,:)'\res_vec;
%         slope = polyfit(, 1)
%          slope = [ones(length(res_vec),1), res_vec]\bio_z(j,:)';
%     
%         title(['r =' num2str(R(1,2),2) ' (p =  ' num2str(P(1,2),2) ')' ' (s =  ' num2str(slope,2) ')'])
%         k = k+1;
%         fprintf('    &    %.2f (p = %.2f) (s = %.2f)', R(1,2), P(1,2), slope)
%     end
%     
%     fprintf('\n')
% end
ii_min = sub2ind([12,6],[2,4,8,8,8,8],[1,2,3,4,5,6]);
res_vec = gaze_results(ii_min(1),:)';
mdl = fitlm(bio_z', res_vec);
disp(mdl);

f = figure(5)
f.Position = [10 50 1000 900];
k = 1;
ca = {'C5', 'C9', 'C14', 'C18', 'C23', 'C135'};
params = {'', 'Q_a_c', '', '', 'R_a_l', '', '', 'ACD', 'LT', 'AL', ''};
xpos = [-3, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5];
ypos1 = [0.13, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
ypos2 = [0.05, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4];
bes = {[2, 8], [2,8], [2,8], [2,8], [2,8], [2,8]};

% for j = 1:11
p= polyfit(bio_z(2,:),gaze_results(ii_min(1),:),2)
plot(bio_z(2,:),gaze_results(ii_min(1),:), '*') 
    for i = 1:6
        res_vec = gaze_results(ii_min(i),:)';
%          [R,P] = corrcoef(min_val_p(:,i), bio_z(j,:)) ;
        ax = subplot(3,2,k);
%         scatter(bio_z(j,:), min_val_p(:,i), '.')
%         h1 = lsline(ax);
%         h1.LineWidth = 1;
%         slope = bio_z(j,:)'\gaze_results(:,i);
        mdl = fitlm(bio_z(bes{i},:)', res_vec);
        h = plot(mdl,'color', 'b');
        hLeg = legend('example');
        set(hLeg,'visible','off')
        set(h(1), 'Color', 'k')
        set(h(2), 'Color', 'b')
        set(h(3), 'visible','off')
        if length(h)>3
            set(h(4), 'visible','off')
        end
        
        title([ca{i}],'FontSize',12)
        str = num2str(mdl.Coefficients.Estimate(1), 2);
        for q =  1:(length(bes{i}))
            b = bes{i};
            if mdl.Coefficients.Estimate(q+1) > 0
                symb = ' + ' ;
            else
                symb = ' - ' ;
            end
            str = [str symb num2str(abs(mdl.Coefficients.Estimate(q+1)), 2 )  params{b(q)}];
        end
        text(xpos(i),ypos1(i),['Error =  ' str],'FontSize',10)
        text(xpos(i),ypos2(i),['R^2 =  ' num2str(mdl.Rsquared.ordinary, '%.2f')],'FontSize',10)
        xlabel('Adjusted whole model (SD)','FontSize',12)
        ylabel(['Adjusted error (' char(176) ')'] ,'FontSize',12)
        ylim([0, 0.6])
        xlim([-4, 1.7])
        k = k+1;
    end
% end
% print(['C:\Users\Joshua\OneDrive\Thesis\LaTex\figs\stochastic_results' '\blig14_bio.png'], '-dpng', '-r600')

% bio_result_z = res_vec_z(:)\bio_mat

% std_result = std_vec'\abs(bio_z)'


% figure(6)
% hold on
% [B,I] = sort(abs(bio_z(2,:)))
% plot(abs(bio_z(2,I)), mean_vec(I), '-')
% figure(7)
% plot(bio_z(2,:), std_vec', 'o')
% 
% R = corrcoef(abs(bio_z(11,:)),std_vec')

fprintf('xxxxxxxxxxxxxxxxxx\n')






