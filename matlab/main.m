clear;
clc;


%% Начальная и конечная точка на сфере Земли в ГСК
m0 = [55.7500, 37.6167]; % МСК
m1 = [59.9500, 30.3167];% СПБ

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



%% PLOTS
% для картинок
monitor = 2;
sz = get (0, 'MonitorPositions');
pic_size = [0 0 1000, 500];
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
plot3 (e0(1), e0(2), e0(3), '*g', 'linewidth',LW);
hold on
plot3 (e1(1), e1(2), e1(3), '*b', 'linewidth',LW);
view (50, 0);
title ('(а)')

