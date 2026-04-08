clear;
clc;


%% Блок начальных условий =================================================

initial_data = struct;   % Структура начальных условий
initial_data.m0 = [55.7558, 37.6173]; % МСК
% initial_data.m1 = [51.1, 0.12];% Лондон
initial_data.m1 = [13.7563, 100.5018]; % Бангкок
initial_data.W = 250; % Модуль вектора путевой скорости [м/с ]
initial_data.t = 60 * 60; % Время в движении [sec]
initial_data.h = 1000000; % высота движения
initial_data.Hz = 10; % Герцовка записи данных

const_values = struct; % Структура константных значений
const_values.W_earth = 7.2921159e-5; % [1/sec]
const_values.R_earth = 6371000; % R Земли [m]

output_data = struct; % Структура выходных (итоговых данных)



%% Подготовка дополнительных данных =======================================
% Градусы в радианы
initial_data.m0 = deg2rad(initial_data.m0);
initial_data.m1 = deg2rad(initial_data.m1);

% переводим широты и долготы в декартовы координаты
initial_data.e0 = ...
    (const_values.R_earth +initial_data.h) * [cos(initial_data.m0(1)) * cos(initial_data.m0(2));
                                              cos(initial_data.m0(1)) * sin(initial_data.m0(2));
                                              sin(initial_data.m0(1))];


initial_data.e1 = ...
    (const_values.R_earth +initial_data.h) * [cos(initial_data.m1(1)) * cos(initial_data.m1(2));
                                              cos(initial_data.m1(1)) * sin(initial_data.m1(2));
                                              sin(initial_data.m1(1))];


% Угол поворота
alpha = acos (dot(initial_data.e0, initial_data.e1) / ...
    (norm(initial_data.e0) * norm(initial_data.e1)));

% вектор поворота
e_c = cross(initial_data.e0, initial_data.e1); e_c = e_c / norm(e_c);
e_c = e_c;

% Длина маршрута:
L = (const_values.R_earth + initial_data.h) * alpha;

% Время в пути:
T = L / initial_data.W;

% Угловая скорость
wc = initial_data.W / (const_values.R_earth + initial_data.h);

% Шаг времени
dt = 0:(1 / initial_data.Hz):T;

C = [0,      -e_c(3),  e_c(2);
     e_c(3),  0,       -e_c(1);
     -e_c(2), e_c(1),  0];

% Определяем шаг угловой скорости вращения
step_beta = wc * (1 / initial_data.Hz);

% Создаём итоговые массивы выходных данных
output_data.r_m_gisk = zeros (length(0:step_beta:alpha), 3);
output_data.r_m_isk = zeros (length(0:step_beta:alpha), 3);
output_data.v_m_gisk = zeros (length(0:step_beta:alpha), 3);
output_data.v_m_isk = zeros (length(0:step_beta:alpha), 3);

% Матрица поворота с учётом углового шага
A = (eye(3) * cos(step_beta)) + ((1 - cos(step_beta))*(e_c * e_c')) - (C * sin(step_beta));

% Перевод начальной и конечной точки в ИСК
initial_data.e0_isk =[+cos(const_values.W_earth * 0), -sin(const_values.W_earth * 0), 0;
                      +sin(const_values.W_earth * 0), +cos(const_values.W_earth * 0), 0;
                      0, 0, 1] * initial_data.e0;

initial_data.e1_isk =[+cos(const_values.W_earth * T), -sin(const_values.W_earth * T), 0;
                      +sin(const_values.W_earth * T), +cos(const_values.W_earth * T), 0;
                      0, 0, 1] * initial_data.e1;

%% main cycle

% % Матрица пересчёта ГИСК -> ИСК
% M_IE = [+cos(const_values.W_earth * dt), -sin(const_values.W_earth * dt), 0;
%         +sin(const_values.W_earth * dt), +cos(const_values.W_earth * dt), 0;
%         0, 0, 1];


for count = 1:length(output_data.r_m_gisk)

    if (count == 1)
        output_data.r_m_gisk(count, :) =  initial_data.e0';
        output_data.r_m_isk(count, :) = [+cos(const_values.W_earth * dt(count)), -sin(const_values.W_earth * dt(count)), 0;
                                         +sin(const_values.W_earth * dt(count)), +cos(const_values.W_earth * dt(count)), 0;
                                         0, 0, 1] * output_data.r_m_gisk(count, :)';

        % output_data.v_m_gisk(count, :) = cross(e_c, output_data.r_m_gisk(count, :)); 
        % output_data.v_m_gisk(count, :) = output_data.v_m_gisk(count, :) / norm (output_data.v_m_gisk(count, :));
        % output_data.v_m_gisk(count, :) = output_data.v_m_gisk(count, :) * initial_data.W;
    else
        output_data.r_m_gisk(count, :) = A' * output_data.r_m_gisk(count-1, :)';
        output_data.r_m_isk(count, :) = [+cos(const_values.W_earth * dt(count)), -sin(const_values.W_earth * dt(count)), 0;
                                         +sin(const_values.W_earth * dt(count)), +cos(const_values.W_earth * dt(count)), 0;
                                         0, 0, 1] * output_data.r_m_gisk(count, :)';

        % output_data.v_m_gisk(count, :) = cross(e_c, output_data.r_m_gisk(count, :)); 
        % output_data.v_m_gisk(count, :) = output_data.v_m_gisk(count, :) / norm (output_data.v_m_gisk(count, :));
        % output_data.v_m_gisk(count, :) = output_data.v_m_gisk(count, :) * initial_data.W;
    end

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

% Траектория полёта ЛА
F0 = figure ('Position', pic_size);

subplot (1, 2, 1) % --------------------------------------------- ГИСК
[x,y,z] = sphere;
axis equal
surf(const_values.R_earth * x, ...
     const_values.R_earth * y, ...
     const_values.R_earth * z, ...
     'facecolor','#404040', 'facealpha',.5); 
grid on;
axis equal
hold on;

% вектор начала маршрута
plot3 ([0, initial_data.e0(1)], ...
       [0, initial_data.e0(2)], ...
       [0, initial_data.e0(3)], ...
       '-*g', 'linewidth',LW);

% вектор конца маршрута
plot3 ([0, initial_data.e1(1)], ...
       [0, initial_data.e1(2)], ...
       [0, initial_data.e1(3)], ...
       '-*b', 'linewidth',LW);

% % нормаль вращения
% plot3 ([0, const_values.R_earth * e_c(1)], ...
%        [0, const_values.R_earth * e_c(2)], ...
%        [0, const_values.R_earth * e_c(3)], ...
%        '-*r', 'linewidth',LW);

% Траектория движения ГИСК
plot3 (output_data.r_m_gisk(:, 1), ...
       output_data.r_m_gisk(:, 2), ...
       output_data.r_m_gisk(:, 3), ...
       '-r', 'linewidth',LW);

title ('Траектория движения в ГИСК')
legend ("E", "Мск", "Бангкок", "маршрут");

subplot (1, 2, 2) % --------------------------------------------- ИСК
axis equal
surf(const_values.R_earth * x, ...
     const_values.R_earth * y, ...
     const_values.R_earth * z, ...
     'facecolor','#404040', 'facealpha',.5); 
grid on;
axis equal
hold on;

% вектор начала маршрута в ИСК
plot3 ([0, output_data.r_m_isk(1, 1)], ...
       [0, output_data.r_m_isk(1, 2)], ...
       [0, output_data.r_m_isk(1, 3)], ...
       '-gs', 'linewidth',LW);

% вектор конца маршрута в ИСК
plot3 ([0, output_data.r_m_isk(end, 1)], ...
       [0, output_data.r_m_isk(end, 2)], ...
       [0, output_data.r_m_isk(end, 3)], ...
       '-*b', 'linewidth',LW);

% Траектория движения ИСК
plot3 (output_data.r_m_isk(:, 1), ...
       output_data.r_m_isk(:, 2), ...
       output_data.r_m_isk(:, 3), ...
       '-r', 'linewidth',LW);

title ('Траектория движения в ИСК')
legend ("E", "Мск", "Бангкок", "маршрут");



clear A C dt e_c L Ln T wc LW monitor M_IE x y z sz alpha
clear step_beta save_graf_trigger pics_path pic_size F0 count