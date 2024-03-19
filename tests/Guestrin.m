
function gaze = Guestrin(pupil, glint ,targets, camera_pos, screen_pos, led_pos)

p0 = [7.8, 4.2, 5, 1.5];
% p0 = [8, 4.4];

func = @(p) calibrate(p,pupil(1:9,:), glint(1:18,:) ,targets(1:9,:), camera_pos(1:9,:), screen_pos(1:9,:), led_pos(1:18,:));

options = optimset('fmincon');
options.Display = 'off';

lb = [3, 2, -10, -5];
ub = [20,15,10,5];

p = fmincon(func,p0,[],[],[],[],lb,ub,[],options);

error = zeros(1,size(pupil,1));
for i = 1:size(pupil,1)
    ii_led = ((i-1)*2+1):i*2;
    error(i) = execute_gaze(p, pupil(i,:), glint(ii_led,:) ,targets(i,:), camera_pos(i,:), screen_pos(i,:), led_pos(ii_led,:));
end

gaze = error;

end

function e = calibrate(calib, pupil, glint, targets, camera_pos, screen_pos, led_pos)


    error = zeros(1,9);
    
    for i = 1:9
        ii_led = ((i-1)*2+1):i*2;
        error(i) = execute_gaze(calib, pupil(i,:), glint(ii_led,:) ,targets(i,:), camera_pos(i,:), screen_pos(i,:), led_pos(ii_led,:));
        if ~isreal(error(i))
            error(i) = 10;
        end
    end
    
    e = sum(error);
    
    if calib(1) > 20
        e = 1000;
    elseif calib(2) > 15
        e = 1000;
    end
end

function [e, c, p, g, oa, va] = execute_gaze(calib, pupil, glint, targets, camera_pos, screen_pos, led_pos)


% pupil_w = convertCameraToEyeCoord(pupil,repelem(camera_pos,1,1));
% plot3(pupil_w(1,1), pupil_w(1,2), pupil_w(1,3), '.m');

% glint_w = convertCameraToEyeCoord(glint,repelem(camera_pos,2,1));
% plot3(glint_w(1,1), glint_w(1,2), glint_w(1,3), '.m');

c = cornea(calib, glint, led_pos, camera_pos);
c_w = convertCameraToEyeCoord(c,repelem(camera_pos,1,1));
% plot3(c_w(1,1), c_w(1,2), c_w(1,3), 'om');

r = refractionPoint(c, pupil, calib(1));
r_w = convertCameraToEyeCoord(r,repelem(camera_pos,1,1));
% plot3(r_w(:,1), r_w(:,2), r_w(:,3), '.r');

[p, f] = calculatePupil(pupil, r, c, calib(2));
p_w = convertCameraToEyeCoord(p,repelem(camera_pos,1,1));
% plot3(p_w(:,1), p_w(:,2), p_w(:,3), 'sk');

gt = (p_w-c_w)./norm(p_w-c_w);
oa = gt./norm(gt);

oa_w = [c_w; c_w+oa.*700];
% plot3(oa_w(:,1), oa_w(:,2), oa_w(:,3), 'b--', 'LineWidth', 1);

va = rotateVector(oa, [0, -calib(4), -calib(3)]);
va_w = [c_w; c_w + va.*700];
% plot3(va_w(:,1), va_w(:,2), va_w(:,3), 'r--', 'LineWidth', 1);

screen_n = [-1 0 0];
tg = dot(screen_n,  [651, -157, -150.5]-c_w, 2 ) ./ dot( va, screen_n, 2 );
g_w = c_w + repmat( tg, 1, 3 ) .* va;

% tg = dot(screen_n,  [0, 150, -35]-c, 2 ) ./ dot( oa, screen_n, 2 );
% g = c + repmat( tg, 1, 3 ) .* oa;

camera_n  = [  0 -0.2377015  0.9713];
t = convertCameraToEyeCoord(camera_n,repelem(camera_pos,1,1));
t2 = [[0 0 0]; camera_n.*100];
t_w = convertCameraToEyeCoord(t2,repelem(camera_pos,1,1));
% plot3(t_w(:,1), t_w(:,2), t_w(:,3), 'b', 'LineWidth', 1);

targets_w = convertCameraToEyeCoord(targets,repelem(camera_pos,1,1));
% plot3(targets_w(:,1), targets_w(:,2), targets_w(:,3), 'om');

% g_w = convertCameraToEyeCoord(g,repelem(camera_pos,1,1));
% plot3(g_w(:,1), g_w(:,2), g_w(:,3), 'xr');


e = atand(sqrt(sum((g_w-targets_w).^2,2))/650);

% e = pupil(self, c, pupil_image)

end

function c = cornea(calib, glint, led_pos, camera_pos)

init = 638;

l1 = led_pos(1,:)./vecnorm(led_pos(1,:),2,2);
l2 = led_pos(2,:)./vecnorm(led_pos(2,:),2,2);
u1 = glint(1,:)./vecnorm(glint(1,:),2,2);
u2 = glint(2,:)./vecnorm(glint(2,:),2,2);

bt = cross(cross(l1, u1), cross(l2, u2));

if(bt(3)<0)
    bt = -1.*bt ;
end

b = bt ./ vecnorm(bt,2,2);
t2 = [[0 0 0]; b.*700];
t_w = convertCameraToEyeCoord(t2,repelem(camera_pos,2,1));
% plot3(t_w(:,1), t_w(:,2), t_w(:,3), 'b', 'LineWidth', 1);

options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt','Display','none','FunctionTolerance', 1*10^-100, 'StepTolerance', 1*10^-10000);
REF = {glint(1, :), b, led_pos(1, :), calib(1), false, camera_pos};
f = @(y) calculateCostCornea(y,REF); % function of dummy variable y
[k_c_res,fval]=fsolve(f,init, options);

% k_c_res = 641.9;
% calculateCostCornea(k_c_res, REF)


c = k_c_res.*b;
c_w = convertCameraToEyeCoord(c,repelem(camera_pos,1,1));
% plot3(c_w(:,1), c_w(:,2), c_w(:,3), '.b');

end

function cost = calculateCostCornea(k_c, var)

ut = var{1};
b = var{2};
led = var{3};
R = var{4};
plot = var{5};
camera_pos = var{6};

u = -1*ut./norm(ut);
c = k_c*b;
B = 2*(dot(-u,c));
C = sum(c.^2) - R^2;
P = [1 B C];
r = roots(P);
kq = min(r);

q = kq*u;
% q_w = convertCameraToEyeCoord(q, camera_pos);
% convertEyeToCameraCoord([-0.6553   -0.1912   -2.3242], camera_pos);
%                                     plot3(q(:,1), q(:,2), q(:,3), '*g');

dis = sqrt(sum((c-q).^2,2));

nn = (q - c)./norm(q - c);
ln = (led - q)./norm(led - q);
lnt = (q - led);


cs1 = dot( nn, ln, 2 );
a = ln - 2 * repmat( cs1, 1, 3 ) .* nn;
a = -1.*a./vecnorm(a,2,2);

%                 R = vrrotvec(a, nn);
%                 angle = R(4)
%
%                 R = vrrotvec(nn, ln);
%                 angle = R(4)
camera_n  =[  0 -0.2377015  0.9713];
t = dot(camera_n, -q, 2 ) ./ dot( a, camera_n, 2 );
rinter = q + repmat( t, 1, 3 ) .* a;

cost = sum((rinter).^2);
if (imag(cost))
    cost = 10000;
end
%             d = sqrt(sum((c).^2,2));
%             if(d < 200)
%                 cost = 10000;
%             end
%
if(k_c < 0)
    cost = 10000;
end
%cost = sum(((2*dot(nn,ln)*nn - ln)-u).^2);
%cost = (dot(a, u))


if plot

    convertCameraToEyeCoord(led, camera_pos)
    t2 = convertCameraToEyeCoord([led; q], camera_pos);

    plot3(t2(:,1), t2(:,2), t2(:,3), '-k', 'LineWidth', 1);
    hold on

    ut = convertCameraToEyeCoord(ut, camera_pos);
    plot3(ut(:,1), ut(:,2), ut(:,3), 'ok');
    plot3(0, 0, 0, 'sk');

    q = convertCameraToEyeCoord(q, camera_pos);
    c = convertCameraToEyeCoord(c, camera_pos);
    a = convertCameraToEyeCoord(a, camera_pos);
    a = a./norm(a);

    plot3(q(:,1), q(:,2), q(:,3), '.k', 'MarkerSize', 20, 'color', [0.173 0.627 0.173]);
    plot3(c(:,1), c(:,2), c(:,3), '.r', 'MarkerSize', 20);

    nn = convertCameraToEyeCoord(nn, camera_pos);
    nn = nn./norm(nn);

    t2 = [c; c+nn.*1];
    plot3(t2(:,1), t2(:,2), t2(:,3), 'b', 'LineWidth', 1);

    t3 = [q; q+a.*k_c];
    plot3(t3(:,1), t3(:,2), t3(:,3), 'k', 'LineWidth', 1, 'color', [0.173 0.627 0.173]);

    xlabel('x')
    ylabel('y')
    zlabel('z')

end


end

% function cost = find_kq(kq, o, u, c, R)
% q = o + kq*(o - u);
% cost = (norm(q-c) - R)^2;
% end

function r = refractionPoint(c, v, R)
v = v./norm(v);
B = 2*(dot(v,c));
C = sum(c.^2) - R^2;

P = [1 B C];
r = roots(P);

kr = min(r);
r = -kr*v;

%             plot3(r(:,1), r(:,2), r(:,3), 'og');
end

function [p, f] = calculatePupil(v, r, c, K)

rn = 1/1.376;
% plot3(v(:,1), v(:,2), v(:,3), '*m');
%rn = self.nrefr ./ rays_out.nrefr; % ratio of in and out refractive indices

v = v./vecnorm(v,2,2);
nrms = (r-c)./vecnorm((r-c),2,2);

cs1 = dot( nrms, v, 2 );

cs2 = sqrt( 1 - rn.^2 .* ( 1 - cs1.^2 ) );
f = repmat( rn, 1, 3 ) .* v - repmat( rn .* cs1 - sign( cs1 ) .* cs2, 1, 3 ) .* nrms;

f = (-1.*f./vecnorm(f,2,2));

B = 2*(dot(f,(r-c)));
C = sum((r-c).^2) - K^2;

P = [1 B C];
root = roots(P);

kf = min(root);

p = r+kf*f;

end


