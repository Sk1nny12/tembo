clc;
clear all;
close all;

%% Исходные данные (вариант c9 = -0.07)
c1 = 0.6;
c2 = 1.8;
c3 = 2;
c4 = 1;
c5 = 0.3;
c6 = 2.5;
c9 = -0.07;
g  = 9.81;

% Вспомогательные коэффициенты знаменателя
a1 = c1 + c4 + c5;        % = 1.9
a0 = c2 + c1*c4;          % = 2.4

% Собственная частота и демпфирование объекта
omega_a = sqrt(a0);
xi_a    = a1/(2*omega_a);

disp('--- Параметры объекта ---');
fprintf('a1 = %.4f, a0 = %.4f\n', a1, a0);
fprintf('omega_alpha = %.4f рад/с, xi_alpha = %.4f\n', omega_a, xi_a);
if a0 > 0
    disp('omega_alpha^2 = c2 + c1*c4 > 0 -> самолёт устойчив по перегрузке');
else
    disp('omega_alpha^2 <= 0 -> самолёт неустойчив по перегрузке');
end

%% Часть 1. Разомкнутое короткопериодическое движение
% Знаменатель D(s) = s^2 + a1*s + a0
D = [1 a1 a0];

% Передаточные функции по координатам (см. модель Lab1_)
% W_omega_z(s) = -[(c3-c5*c9)*s + (c3*c4-c2*c9)] / D(s)
num_wz = -[(c3 - c5*c9), (c3*c4 - c2*c9)];
W_wz   = tf(num_wz, D);

% W_alpha(s) = -[c9*s + (c3 - c1*c9)] / D(s)
num_al = -[c9, (c3 - c1*c9)];
W_al   = tf(num_al, D);

% W_ny_alpha(s) = (-c6/g)*[-c9*s^2 - c9*(c1+c5)*s + (c3*c4 - c2*c9)] / D(s)
num_ny = (-c6/g) * [-c9, -c9*(c1+c5), (c3*c4 - c2*c9)];
W_ny   = tf(num_ny, D);

% W_theta(s) = W_omega_z(s) * 1/s
W_th = W_wz * tf(1, [1 0]);

% W_H(s) = W_omega_z(s) * (c6*c4) / [(s+c4)*s^2]
W_H  = W_wz * tf(c6*c4, conv([1 c4], [1 0 0]));

disp('--- Передаточные функции (c9 = -0.07) ---');
W_wz, W_al, W_ny

% Переходные процессы при ступенчатом δв = 1
t1 = 0:0.01:10;

[y_wz, ~] = step(W_wz, t1);
[y_al, ~] = step(W_al, t1);
[y_ny, ~] = step(W_ny, t1);

figure(1);
plot(t1, y_wz, 'b-', t1, y_al, 'r--', t1, y_ny, 'k-.', 'LineWidth', 2);
xlabel('t, с'); ylabel('Координаты');
title('Рис.2. Реакция \omega_z, \alpha, n_{y\alpha} на \delta_в = 1, c_9 = -0.07');
legend('\omega_z (1)', '\alpha (2)', 'n_{y\alpha} (3)');
grid on;

% Траекторные координаты ϑ и H (имеют интеграторы)
t2 = 0:0.01:14;
[y_th, ~] = step(W_th, t2);
[y_H,  ~] = step(W_H,  t2);

figure(2);
yyaxis left;
plot(t2, y_th, 'b-', 'LineWidth', 2); ylabel('\vartheta');
yyaxis right;
plot(t2, y_H, 'r--', 'LineWidth', 2); ylabel('H');
xlabel('t, с');
title('Рис.3. Реакция \vartheta (1) и H (2), c_9 = -0.07');
grid on;

%% Случай c9 = 0 (без вертикальной составляющей ветра)
c9_0 = 0;

num_wz0 = -[(c3 - c5*c9_0), (c3*c4 - c2*c9_0)];
W_wz0   = tf(num_wz0, D);

num_al0 = -[c9_0, (c3 - c1*c9_0)];
W_al0   = tf(num_al0, D);

num_ny0 = (-c6/g) * [-c9_0, -c9_0*(c1+c5), (c3*c4 - c2*c9_0)];
W_ny0   = tf(num_ny0, D);

W_th0   = W_wz0 * tf(1, [1 0]);
W_H0    = W_wz0 * tf(c6*c4, conv([1 c4], [1 0 0]));

[y_wz0, ~] = step(W_wz0, t1);
[y_al0, ~] = step(W_al0, t1);
[y_ny0, ~] = step(W_ny0, t1);

figure(3);
plot(t1, y_wz0, 'b-', t1, y_al0, 'r--', t1, y_ny0, 'k-.', 'LineWidth', 2);
xlabel('t, с'); ylabel('Координаты');
title('Рис.4. Реакция \omega_z, \alpha, n_{y\alpha} на \delta_в = 1, c_9 = 0');
legend('\omega_z (1)', '\alpha (2)', 'n_{y\alpha} (3)');
grid on;

[y_th0, ~] = step(W_th0, t2);
[y_H0,  ~] = step(W_H0,  t2);

figure(4);
yyaxis left;
plot(t2, y_th0, 'b-', 'LineWidth', 2); ylabel('\vartheta');
yyaxis right;
plot(t2, y_H0, 'r--', 'LineWidth', 2); ylabel('H');
xlabel('t, с');
title('Рис.5. Реакция \vartheta (1) и H (2), c_9 = 0');
grid on;

% Совмещённые переходные функции c9 = -0.07 (1) и c9 = 0 (2)
figure(5);
subplot(3,1,1);
plot(t1, y_wz, 'b-', t1, y_wz0, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('\omega_z');
title('Рис.6. \omega_z: c_9 = -0.07 (1) и c_9 = 0 (2)');
legend('c_9 = -0.07', 'c_9 = 0'); grid on;

subplot(3,1,2);
plot(t1, y_al, 'b-', t1, y_al0, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('\alpha');
title('Рис.7. \alpha: c_9 = -0.07 (1) и c_9 = 0 (2)');
legend('c_9 = -0.07', 'c_9 = 0'); grid on;

subplot(3,1,3);
plot(t1, y_ny, 'b-', t1, y_ny0, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('n_{y\alpha}');
title('Рис.8. n_{y\alpha}: c_9 = -0.07 (1) и c_9 = 0 (2)');
legend('c_9 = -0.07', 'c_9 = 0'); grid on;

% Показатели качества (Peak, Steady, σ, tпп, tн) для c9 = -0.07
disp('--- Показатели качества (c9 = -0.07) ---');
print_quality('omega_z', t1, y_wz);
print_quality('alpha  ', t1, y_al);
print_quality('ny_a   ', t1, y_ny);

%% Часть 2. СУУ с демпфером тангажа
disp(' ');
disp('=== Часть 2. СУУ с демпфером тангажа ===');

% Условие управляемости xi_d = 0.707 даёт квадратное уравнение
% (c1+c4+c5+Kwz*c3)^2 = 2*(c1*c4 + c2 + Kwz*c3*c4)
% => c3^2*Kwz^2 + 2*c3*(c1+c5)*Kwz + (c1+c4+c5)^2 - 2*(c1*c4+c2) = 0
% Для данных параметров: 4*Kwz^2 + 3.6*Kwz - 1.19 = 0
A_q = c3^2;
B_q = 2*c3*(c1+c5);
C_q = (c1+c4+c5)^2 - 2*(c1*c4 + c2);

K_roots = roots([A_q B_q C_q]);
fprintf('Корни уравнения для K_wz: %.4f, %.4f\n', K_roots(1), K_roots(2));

% Выбираем положительный корень
Kwz = K_roots(K_roots > 0);
if isempty(Kwz)
    Kwz = max(K_roots);
end
fprintf('Принятое K_wz = %.4f\n', Kwz);

% Проверка по системе (2): требуемая частота и демпфирование
two_xi_Omega = c1 + c4 + c5 + Kwz*c3;
Omega_d_sq   = c1*c4 + c2 + Kwz*c3*c4;
Omega_d      = sqrt(Omega_d_sq);
xi_d         = two_xi_Omega/(2*Omega_d);
fprintf('Omega_d = %.4f рад/с, xi_d = %.4f\n', Omega_d, xi_d);

% Коэффициент АРУ K1 (DC-усиление замкнутого контура по nyα = 1)
W_wz_dc = dcgain(W_wz);
W_ny_dc = dcgain(W_ny);
W_pr_dc = 1;     % привод имеет единичное усиление по постоянному входу

K1 = (1 - Kwz*W_pr_dc*W_wz_dc) / (W_pr_dc*W_ny_dc);
fprintf('Коэффициент АРУ K1 = %.4f (|K1| = %.4f)\n', K1, abs(K1));

% Привод: W_pr(s) = omega_pr^2 / (s^2 + 2*xi_pr*omega_pr*s + omega_pr^2)
omega_pr = 20;
xi_pr    = 0.707;
W_pr     = tf(omega_pr^2, [1, 2*xi_pr*omega_pr, omega_pr^2]);

%% Замкнутая система СУУ: x -> K1 -> W_pr -> самолёт; ОС по omega_z через Kwz
% Внутренний контур: с обратной связью по omega_z
% W_inner_wz(s) = W_pr*W_wz / (1 - Kwz*W_pr*W_wz)
% Знак "-" перед Kwz взят с учётом, что сигнал ОС вычитается на сумматоре
% (см. рис.9); в модели Wωz<0, поэтому минус обеспечивает отрицательную ОС.

build_loop = @(Wpr) struct( ...
    'wz', feedback(Wpr*W_wz, -Kwz), ...
    'al', feedback_branch(Wpr, W_wz, W_al, Kwz), ...
    'ny', feedback_branch(Wpr, W_wz, W_ny, Kwz) );

% а) Идеальный привод, демпфер выключен (Kwz = 0)
G_ideal_nodamp_wz = W_wz;                    % контур разорван
G_ideal_nodamp_al = W_al;
G_ideal_nodamp_ny = K1 * W_ny;               % K1 учитываем во всех режимах

% в) Идеальный привод, демпфер включён
L_ideal = build_loop(tf(1,1));
G_ideal_damp_wz = K1 * L_ideal.wz;
G_ideal_damp_al = K1 * L_ideal.al;
G_ideal_damp_ny = K1 * L_ideal.ny;

% г) Реальный привод, демпфер включён
L_real = build_loop(W_pr);
G_real_damp_wz = K1 * L_real.wz;
G_real_damp_al = K1 * L_real.al;
G_real_damp_ny = K1 * L_real.ny;

t3 = 0:0.005:6;

% Рис.10-12: без демпфера (1) и с демпфером (2), W_pr = 1
[y1_wz, ~] = step(K1*W_wz, t3);     % без демпфера
[y2_wz, ~] = step(G_ideal_damp_wz, t3);
[y1_al, ~] = step(K1*W_al, t3);
[y2_al, ~] = step(G_ideal_damp_al, t3);
[y1_ny, ~] = step(K1*W_ny, t3);
[y2_ny, ~] = step(G_ideal_damp_ny, t3);

figure(6);
subplot(3,1,1);
plot(t3, y1_wz, 'b-', t3, y2_wz, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('\omega_z');
title('Рис.10. \omega_z: без демпфера (1) и с демпфером (2), W_{пр} = 1');
legend('K_{\omega_z} = 0', 'K_{\omega_z} \neq 0'); grid on;

subplot(3,1,2);
plot(t3, y1_al, 'b-', t3, y2_al, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('\alpha');
title('Рис.11. \alpha: без демпфера (1) и с демпфером (2), W_{пр} = 1');
legend('K_{\omega_z} = 0', 'K_{\omega_z} \neq 0'); grid on;

subplot(3,1,3);
plot(t3, y1_ny, 'b-', t3, y2_ny, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('n_{y\alpha}');
title('Рис.12. n_{y\alpha}: без демпфера (1) и с демпфером (2), W_{пр} = 1');
legend('K_{\omega_z} = 0', 'K_{\omega_z} \neq 0'); grid on;

% Рис.13-15: идеальный (1) и реальный (2) привод, демпфер включён
[y3_wz, ~] = step(G_real_damp_wz, t3);
[y3_al, ~] = step(G_real_damp_al, t3);
[y3_ny, ~] = step(G_real_damp_ny, t3);

figure(7);
subplot(3,1,1);
plot(t3, y2_wz, 'b-', t3, y3_wz, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('\omega_z');
title('Рис.13. \omega_z: идеальный (1) и реальный (2) привод, K_{\omega_z} \neq 0');
legend('W_{пр} = 1', 'W_{пр}(s)'); grid on;

subplot(3,1,2);
plot(t3, y2_al, 'b-', t3, y3_al, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('\alpha');
title('Рис.14. \alpha: идеальный (1) и реальный (2) привод, K_{\omega_z} \neq 0');
legend('W_{пр} = 1', 'W_{пр}(s)'); grid on;

subplot(3,1,3);
plot(t3, y2_ny, 'b-', t3, y3_ny, 'r--', 'LineWidth', 2);
xlabel('t, с'); ylabel('n_{y\alpha}');
title('Рис.15. n_{y\alpha}: идеальный (1) и реальный (2) привод, K_{\omega_z} \neq 0');
legend('W_{пр} = 1', 'W_{пр}(s)'); grid on;

% Показатели качества СУУ
disp(' ');
disp('--- Показатели качества СУУ ---');
disp('Wпр=1, Kwz=0:');
print_quality('ny_a   ', t3, y1_ny);
print_quality('omega_z', t3, y1_wz);
print_quality('alpha  ', t3, y1_al);
disp('Wпр=1, Kwz!=0:');
print_quality('ny_a   ', t3, y2_ny);
print_quality('omega_z', t3, y2_wz);
print_quality('alpha  ', t3, y2_al);
disp('Wпр(s), Kwz!=0:');
print_quality('ny_a   ', t3, y3_ny);

%% Выводы
disp(' ');
disp('=== Выводы ===');
disp('1. Самолёт устойчив по перегрузке (omega_a^2 = c2+c1*c4 > 0).');
fprintf('2. Собственная частота omega_a = %.3f рад/с, демпфирование xi_a = %.3f.\n', omega_a, xi_a);
fprintf('3. Требуемый коэффициент демпфера K_wz = %.4f обеспечивает xi_d = %.3f, Omega_d = %.3f рад/с.\n', Kwz, xi_d, Omega_d);
fprintf('4. Коэффициент АРУ K1 = %.3f обеспечивает единичный градиент по перегрузке.\n', K1);
disp('5. Демпфер тангажа снижает перерегулирование и почти вдвое сокращает время переходного процесса.');
disp('6. Динамика реального привода (omega_pr = 20 рад/с) практически не ухудшает качество.');

%% Вспомогательные функции
function L = feedback_branch(Wpr, Wwz, Wout, Kwz)
    % Замкнутый по omega_z контур: вход -> Wpr -> объект (выход Wout),
    % ОС по omega_z через Kwz. Возвращает ПФ от входа до Wout.
    L = Wpr * Wout / (1 - Kwz*Wpr*Wwz);
end

function print_quality(name, t, y)
    [peak, idx] = max(abs(y));
    peak = peak * sign(y(idx));
    steady = y(end);
    if abs(steady) > 1e-9
        sigma = (peak - steady)/steady * 100;
    else
        sigma = NaN;
    end
    % Время регулирования по 5%
    band = 0.05 * abs(steady);
    tpp  = NaN;
    for k = length(y):-1:1
        if abs(y(k) - steady) > band
            if k < length(y)
                tpp = t(k+1);
            end
            break;
        end
    end
    % Время нарастания 0..95%
    tn = NaN;
    target = 0.95 * steady;
    for k = 1:length(y)
        if (steady >= 0 && y(k) >= target) || (steady < 0 && y(k) <= target)
            tn = t(k);
            break;
        end
    end
    fprintf('  %s: Peak=%+.4f, Steady=%+.4f, sigma=%6.2f %%, tпп=%5.2f с, tн=%5.2f с\n', ...
        name, peak, steady, sigma, tpp, tn);
end
