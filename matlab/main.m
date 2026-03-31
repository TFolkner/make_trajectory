clear;
clc;


%% Начальная и конечная точка на сфере Земли в ГСК
m0 = [45.0, 90.0];
m1 = [0.0, 0.0];

m0 = deg2rad(m0);
m1 = deg2rad(m1);

t = 1; % Время в движении [sec]

W_earth = 7.2921159e-5; % [1/sec]
R_earth = 1; % R Земли [m]
M_IE = [+cos(W_earth * t), -sin(W_earth * t), 0;
        +sin(W_earth * t), +cos(W_earth * t), 0;
        0, 0, 1];

%% main cycle

% Орты нормалей в ГСК
e0 = R_earth * [cos(m0(1)) * cos(m0(2));
      cos(m0(1)) * sin(m0(2));
      sin(m0(1))];


e1 = R_earth * [cos(m1(1)) * cos(m1(2));
      cos(m1(1)) * sin(m1(2));
      sin(m1(1))];


% Орты нормалей в ИСК
e0_ = M_IE * e0;
e1_ = M_IE * e1;


% Угол поворота
alpha = acos (dot(e0_, e1_) / (norm(e0_) * norm(e1_)));
e_c = cross(e0_, e1_) /(norm(e0_) * norm(e1_));


% поворот на половину
C = [0, -e_c(3), e_c(2);
     e_c(3), 0, -e_c(1);
     -e_c(2), e_c(1), 0];
E = eye(3);

H = 0.1;
beta = 0:H:alpha;
e_m = zeros (length(beta), 3);
for i = 1:length(beta)
    A = E * cos(beta(i)) + (1 - cos(beta(i))) * (e_c * e_c') - C * sin(beta(i));
    e_m(i, :) = A' * e0_;
end




%% PLOTS
% для картинок
monitor = 2;
sz = get (0, 'MonitorPositions');
pic_size = [0 0 2000, 1000];
pics_path = "../pics/";
save_graf_trigger = 0;
LW = 1.5;

% Траектория полёта КА
F0 = figure ('Position', pic_size);
[x,y,z] = sphere;
axis equal
surf(R_earth*x, R_earth*y, R_earth*z, 'facecolor','#404040', 'facealpha',.7); 
grid on;
axis equal
hold on; 
plot3 ([0, e0_(1)], [0, e0_(2)], [0, e0_(3)], '-*g', 'linewidth',LW);
hold on
plot3 ([0, e1_(1)], [0, e1_(2)], [0, e1_(3)], '-*b', 'linewidth',LW);
plot3 ([0, e_c(1)], [0, e_c(2)], [0, e_c(3)], '-*r', 'linewidth',LW);

view (90, 0);
title ('(а)')

plot3 (e_m(:, 1), e_m(:, 2), e_m(:, 3), '*y', 'linewidth',LW);

