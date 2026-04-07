clear;
clc;


%% Блок начальных условий =================================================

initial_data = struct;   % Структура начальных условий
initial_data.m0 = [55.7558, 37.6173]; % МСК
% initial_data.m1 = [51.1, 0.12];% Лондон
initial_data.m1 = [13.7563, 100.5018]; % Бангкок
initial_data.W = 250; % Модуль вектора путевой скорости [м/с ]
initial_data.t = 60 * 60; % Время в движении [sec]
initial_data.h = 0.1; % высота движения
initial_data.Hz = 10; % Герцовка записи данных

const_values = struct; % Структура константных значений
const_values.W_earth = 7.2921159e-5; % [1/sec]
const_values.R_earth = 6371000; % R Земли [m]


output_data = struct; % Структура выходных (итоговых данных)



dt = 1;
M_IE = [+cos(const_values.W_earth * dt), -sin(const_values.W_earth * dt), 0;
        +sin(const_values.W_earth * dt), +cos(const_values.W_earth * dt), 0;
        0, 0, 1];







%% main cycle

initial_data.m0 = deg2rad(initial_data.m0);
initial_data.m1 = deg2rad(initial_data.m1);

% Орты нормалей в ГСК
e0 = [cos(initial_data.m0(1)) * cos(initial_data.m0(2));
      cos(initial_data.m0(1)) * sin(initial_data.m0(2));
      sin(initial_data.m0(1))];


e1 = [cos(initial_data.m1(1)) * cos(initial_data.m1(2));
      cos(initial_data.m1(1)) * sin(initial_data.m1(2));
      sin(initial_data.m1(1))];

% Угол поворота
alpha = acos (dot(e0, e1) / (norm(e0) * norm(e1)));
e_c = cross(e0, e1); e_c = e_c / norm(e_c);

% Длина маршрута:
L = (const_values.R_earth + initial_data.h) * alpha;

% Время в пути:
T = L / initial_data.W;

% Угловая скорость
wc = initial_data.W / (const_values.R_earth + initial_data.h);

% Шаг времени
dt = 1 / initial_data.Hz;

% Общаяя длинна вектора
Ln = int64(T / dt);

% поворот 
C = [0,      -e_c(3),  e_c(2);
     e_c(3),  0,       -e_c(1);
     -e_c(2), e_c(1),  0];


step_beta = wc * dt;
e_m = zeros (length(0:step_beta:alpha), 3);
w_m = zeros (length(0:step_beta:alpha), 3);

A = (eye(3) * cos(step_beta)) + ((1 - cos(step_beta))*(e_c * e_c')) - (C * sin(step_beta));

e_m(1, :) = e0;
for i = 2:length(e_m)
    e_m(i, :) = A' * e_m(i-1, :)';

    w_m(i, :) = cross(e_c, e_m(i, :)); 
    w_m(i, :) = w_m(i, :) / norm (w_m(i, :));
    w_m(i, :) = w_m(i, :) * initial_data.W;
end




%% ВЫвод инфо =============================================================
fprintf ("Путевая скорость - %.4f [m/s] \n", initial_data.W);

fprintf ("Длина маршрута - %.4f [m] \n", L);
fprintf ("               - %.4f [km] \n", L/1000);

fprintf ("Время в пути - %.4f [s]\n", T);
fprintf ("             - %.4f [m]\n", T/60);
fprintf ("             - %.4f [h]\n", T/3600);




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
surf(x, y, z, 'facecolor','#404040', 'facealpha',.5); 
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
plot3 (e_m(:, 1), e_m(:, 2), e_m(:, 3), '-r', 'linewidth',LW);


legend ("E", "Мск", "Лондон", "Нормаль", "маршрут");