clc; 
clear all; 
close all;

E=[1 0 0;0 1 0;0 0 1];

%коэффициенты
A = [-0.32 -2.1; 0.5 -0.74];
B = [-1.5; 0.17];
c = 0.7;
d = 0.15;
T = 0.1;


X0_1 = [0; 5; 4];  

%вычисления для 1 части (устойчивость и управляемость)
flag1 = 1;

root = roots(poly(A));

disp("корни:");
disp(root);

for i = 1:length(root)
    if real(root(i)) > 0  
        flag1 = 0;
        disp("система не устойчива : есть положительный корень")
        break;
    end
end

if flag1 == 1
    disp("система устойчива : все корни отрицательны")
end    

P = [B A*B];
disp("ранг матрицы управляемости : ");
disp(rank(P));
disp("размерность матрицы P : ");
disp(size(P));

disp("N (количество переменных) : ");
disp(size(A, 1));

if rank(P) == size(A, 1)
    disp("система управляема");
else 
    disp("система не управляема");
end  

%графики 1 части

u = 4;

time = 0:0.001:10; % от 0 до 10 с шагом 0.1
t_span = [0 14]; 

[~, Y] = ode45(@(t, y) A*y + B*u, time, [0; 5]);

w_vals = Y(:, 1);
a_vals = Y(:, 2);

% Создаем новый рисунок
figure(1);

% Построение графика для w
subplot(2,1,1);
plot(time, w_vals, 'b-', 'LineWidth', 2);  
xlabel('Время');                             
ylabel('w');
title('значение скорости');
grid on;                                     

% Построение графика для a
subplot(2,1,2);
plot(time, a_vals, 'r-', 'LineWidth', 2);   
xlabel('Время');                            
ylabel('a');
title('значение угла');
grid on;                                     

%часть 2
A = [-0.32 -2.1 B(1); 0.5 -0.74 B(2);0 0 -1/T];
B = [0; 0;1/T];

T_m=([1 0 0;0 c d;0 0 1])^-1;

A_1=(T_m^-1)*A*T_m;
B_1=(T_m^-1)*B;

psi = 0.7;
omega = 3;

a1_1 = 2*psi*omega + 1/T;
a2_1 = omega^2+2*psi*omega*1/T;
a3_1 = 1/T * omega^2;

P = [B_1 A_1*B_1 A_1^2*B_1];

q=A_1*A_1*A_1+a1_1*A_1*A_1+a2_1*A_1+a3_1*E;

K_oc=-[0 0 1]*P^-1*q;

A_cl = A_1 + B_1 * K_oc;

% Исследование устойчивости замкнутой системы
eigenvals_cl = eig(A_cl);

% Проверка устойчивости
if all(real(eigenvals_cl) < 0)
    disp('Замкнутая система устойчива');
else
    disp('Замкнутая система НЕ устойчива');
end

% Проверка управляемости замкнутой системы
P_cl = [B_1 A_cl*B_1 A_cl^2*B_1];
if rank(P_cl) == size(A_cl, 1)
    disp('Замкнутая система управляема');
else
    disp('Замкнутая система НЕ управляема');
end


t_span = 0:0.01:10;
g = 4;  % задающее воздействие от ручки


[t, X] = ode45(@(t, x) closed_loop_system(t, x, A_1, B_1, K_oc, g), t_span, X0_1);

figure(2)
plot(t,X)

figure(3)

subplot(3,1,1);
plot(t,X(:,1))
xlabel('Время');                            
ylabel('w');
title('угловая скорость');
grid on;    

subplot(3,1,2);
plot(t,X(:,2))
xlabel('Время');                            
ylabel('NY');
title('NY');
grid on;    

subplot(3,1,3);
plot(t,X(:,3))
xlabel('Время');                            
ylabel('delta');
title('delta');
grid on;    
%часть 3 наблюдатель (ИСПРАВЛЕНО)
disp("_____________________________");
disp("часть 3 наблюдатель")
disp("_____________________________");

C = [1 0 0];
omega_n = 15;
psi_n = 0.5;

% ИСПРАВЛЕНО: правильная матрица наблюдаемости
H = [C' , (A_1')*C' , ((A_1')^2)*C'];

if rank(H) == 3
    disp("система наблюдаема");
    disp(['rank(H) = ', num2str(rank(H))]);
else 
    disp("система не наблюдаема")
end

An = A_1;
Bn = B_1;
Cn = C;

a1_n = 2*psi_n*omega_n + 1/T;
a2_n = omega_n^2 + 2*psi_n*omega_n*1/T;
a3_n = 1/T * omega_n^2;


q = (An')^3 + a1_n*(An'^2) + a2_n*An' + a3_n*E;

L = ([0 0 1] * inv(H) * q)';

disp('Коэффициенты наблюдателя L:');
disp(L);


disp("_____________________________");
disp("часть 4 динамика + наблюдатель")
disp("_____________________________");


g = 4;  

% Проверка устойчивости матрицы системы с наблюдателем
m1 = [A_1, B_1*K_oc; L*C, An + Bn*K_oc - L*Cn];
eig_m1 = eig(m1)



if all(real(eig_m1) < 0)
    disp('система с наблюдателем устойчива');
else
    disp('система с наблюдателем не  устойчива');
end

t_span = 0:0.01:2;


[t, X] = ode45(@(t, x) observer_system(t, x, A_1, B_1, K_oc, L, C, g), t_span, [X0_1; 0; 0; 0]);

figure(4)
hold on;
plot(t, X(:,1), 'b-', 'LineWidth', 1.5);
plot(t, X(:,2), 'r-', 'LineWidth', 1.5);
plot(t, X(:,3), 'g-', 'LineWidth', 1.5);
plot(t, X(:,4), 'b--', 'LineWidth', 1.5);
plot(t, X(:,5), 'r--', 'LineWidth', 1.5);
plot(t, X(:,6), 'g--', 'LineWidth', 1.5);
hold off;

xlabel('Время (с)');
ylabel('Состояния');
title('Состояния системы и оценки наблюдателя');
legend('\omega_z', 'n_y', '\delta', '\omega_z^{оценка}', 'n_y^{оценка}', '\delta^{оценка}');
grid on;




function dx = closed_loop_system(~, x, A_1, B_1, K_oc, g)
    % x = [w; NY; delta]
    u = K_oc * x + g;  % управление вычисляется из состояний
    dx = A_1 * x + B_1 * u;
end

function dx = observer_system(~, x, A_1, B_1, K_oc, L, C, g)
    % x = [w; NY; delta; w_hat; NY_hat; delta_hat]
    
    % Реальный объект
    X_real = x(1:3);
    % Оценки наблюдателя
    X_hat = x(4:6);
    
    % Управление по оценкам
    u = K_oc * X_hat + g;
    
    % Производные для реального объекта
    dX_real = A_1 * X_real + B_1 * u;
    
    % Измеряемый выход
    Y = C * X_real;
    
    % Производные для наблюдателя
    dX_hat = A_1 * X_hat + B_1 * u + L * (Y - C * X_hat);
    
    % Объединение
    dx = [dX_real; dX_hat];
end