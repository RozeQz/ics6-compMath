% Исходные данные
% 1-я строка: компоненты вектора a=(a_2, a_3, ..., a_n) (подддиагональ заполняется со второго элемента)
% 2-я строка: компоненты вектора b=(b_1, b_2, ..., b_n) (ддиагональ заполняется полностью)
% 3-я строка: компоненты вектора c=(c_1, c_2, ..., c_(n-1) ) (наддиагональ заполняется без последнего элемента)
% 4-я строка: компоненты вектора d=(d_1, d_2, ..., d_n) правых частей

a = [2; -1; 1; 1; 0; 1];
b = [104; 98; 88; -141; 81; 73; 67];
c = [1; 0; -1; 1; 1; 1];
d = [22; 20; 19; -26; 16; 15; 13];

A = diag(a, -1) + diag(b, 0) + diag(c, 1);
F = d;
X = A \ F;

disp('Исходная матрица A:');
disp(A)

a = [0; a]; % Добавляем элемент в начало вектора a
c = [c; 0]; % Добавляем элемент в конец вектора с
 
n = size(A, 1);

% Прямая прогонка
[alpha, beta, y] = Direct_Run(a, b, c, d);

fprintf('Выполнение достаточных условий применимости метода прогонки: ');
Check_Sufficient_Condition(alpha, beta, y);

% Обратная прогонка
fprintf('\nОбратная прогонка:');
x = Reverse_Run(alpha, beta);

%Вывод значений
fprintf('\nАльфа: ');
fprintf('%6.3f\t', alpha);
fprintf('\nБета:  ');
fprintf('%6.3f\t', beta);
fprintf('\nХ:     ');
fprintf('%6.3f\t', x);
 
%Невязка решения
fprintf('\n\n');
R = d - A * x;
disp('Невязка:');
disp(R);

%Норма невязки
norm_1 = norm(R, 1)    ;                       
disp('Единичная норма невязки:');
disp(norm_1)
norm_inf = norm(R, inf);  
disp('Бесконечная норма невязки:');
disp(norm_inf)
 
%Устойчивость к возмущениям
d1=d+0.01;

% Прямая прогонка
[alpha, beta, y] = Direct_Run(a, b, c, d1);

fprintf('Выполнение достаточных условий применимости метода прогонки: ');
Check_Sufficient_Condition(alpha, beta, y);

% Обратная прогонка
x_perturbed = Reverse_Run(alpha, beta);

%Вывод значений
fprintf('\nУстойчивость к малым возмущениям:\n');

fprintf('\nАльфа: ');
fprintf('%6.3f\t', alpha);
fprintf('\nБета:  ');
fprintf('%6.3f\t', beta);
fprintf('\nХ:     ');
fprintf('%6.3f\t', x_perturbed);

%Сравнение решений
fprintf('\n\n');
disp('Решение исходной системы:');
disp(x);
disp('Решение возмущенной системы:');
disp(x_perturbed);

% Проверка выполнения достаточных условий применимости метода прогонки
function Check_Sufficient_Condition(a, b, c)
    res = 'TRUE';
    n = size(a, 1);
    for i= 1:n
        if abs(c(i)) < abs(a(i)) + abs(b(i))
            res = 'FALSE';
        end
    end
    disp(res);
end

% Прямая прогонка
function [alpha, beta, y] = Direct_Run(a, b, c, d)
    n = size(a, 1);
    alpha = zeros(n, 1);
    beta = zeros(n, 1);
    y = zeros(n, 1);

    y(1) = b(1);
    alpha(1) = -c(1) / y(1);
    beta(1) = d(1) / y(1);
    for i = 2:n
        y(i) = b(i) + a(i) * alpha(i - 1);
        alpha(i) = -c(i) / y(i);
        beta(i) = (d(i) - a(i) * beta(i - 1)) / y(i);
    end
end

% Обратная прогонка
function x = Reverse_Run(alpha, beta)
    n = size(alpha, 1);
    x = zeros(n, 1);
    x(n) = beta(n);
    for i = 1:(n-1)
        x(n-i) = alpha(n-i) * x(n-i + 1) + beta(n-i);
    end
end