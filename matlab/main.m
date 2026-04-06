clear;
clc;


%% Начальная и конечная точка на сфере Земли в ГСК
m0 = [55.7558, 37.6173]; % МСК
% m1 = [13.7563, 100.5018]; % Тай
m1 = [51.1, 0.12];% Лондон

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

e0 = M_IE * e0;
e1 = M_IE * e1;


% Угол поворота
alpha = acos (dot(e0, e1) / (norm(e0) * norm(e1)));
e_c = cross(e0, e1); e_c = e_c / norm(e_c);


% поворот 
C = [0,      -e_c(3),  e_c(2);
     e_c(3),  0,       -e_c(1);
     -e_c(2), e_c(1),  0];

H = 0.01;
betta = 0:H:alpha;
e_m = zeros (length(betta), 3);
e_m(1, :) = e0;
A = (eye(3) * cos(H)) + ((1 - cos(H))*(e_c * e_c')) - (C * sin(H));

for i = 2:length(betta)
    e_m(i, :) = A' * e_m(i-1, :)';
end



%% PLOTS
% для картинок
monitor = 2;
sz = get (0, 'MonitorPositions');
pic_size = [0 0 1000, 1000];
pics_path = "../pics/";
save_graf_trigger = 0;
LW = 2;

% Траектория полёта КА
F0 = figure ('Position', pic_size);
[x,y,z] = sphere;
axis equal
surf(R_earth*x, R_earth*y, R_earth*z, 'facecolor','#404040', 'facealpha',.5); 
grid on;
axis equal
hold on; 
plot3 ([0, e0(1)], [0, e0(2)], [0, e0(3)], '-*g', 'linewidth',LW);
hold on
plot3 ([0, e1(1)], [0, e1(2)], [0, e1(3)], '-*b', 'linewidth',LW);
plot3 ([0, e_c(1)], [0, e_c(2)], [0, e_c(3)], '-*r', 'linewidth',LW);

view (90, 0);
title ('(а)')

% plot3 ([0, e_m(1)], [0, e_m(2)], [0, e_m(3)], '-*y', 'linewidth',LW);
plot3 (e_m(:, 1), e_m(:, 2), e_m(:, 3), '-y', 'linewidth',LW);


legend ("E", "Мск", "Лондон", "Нормаль", "маршрут");