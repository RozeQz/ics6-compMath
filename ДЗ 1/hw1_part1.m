%Исходные данные
A1_EXT = [-62.8000	-6.9700	7.7300	0.0000	-221.3100; 
          0.0000		49.4000 	-3.5600	-4.5700	-523.8800;
          -7.1100		-5.9700	97.6000	-7.5800	1121.6300; 	
          6.8100		-1.2900	-7.9500	42.2000	294.6700];

ANS1 = [6; -9; 12; 8];

A2_EXT = [-363.94200	-21.30800	4463.83200	738.52800	175.22000;
          -78.20700	-4.57100	959.24700	158.69700	37.97000;
          -0.07200	0.00000	0.88600	0.14400	0.18800;
          -181.53800	-10.65400	2226.58800	368.39800	86.48000];

ANS2 = [10; 32; 2; -6];

A1 = A1_EXT(:, 1:(end-1));
B1 = A1_EXT(:, end);

A2 = A2_EXT(:, 1:(end-1));
B2 = A2_EXT(:, end);

%Метод Гаусса:
disp('Решение СЛАУ №1 реализованным методом Гаусса:');
X1 = Gauss(A1, B1);
disp(X1);

disp('Решение СЛАУ №1 встроенным методом Гаусса:');
X1_MATLAB = A1 \ B1;
disp(X1_MATLAB);

disp('Решение СЛАУ №2 реализованным методом Гаусса:');
X2 = Gauss(A2, B2);
disp(X2);

disp('Решение СЛАУ №2 встроенным методом Гаусса:');
X2_MATLAB = A2 \ B2;
disp(X2_MATLAB);

%Невязка СЛАУ:
disp('Невязка СЛАУ №1:');
R1 = B1 - A1 * X1;
disp(R1);

disp('Невязка СЛАУ №2:');
R2 = B2 - A2 * X2;
disp(R2);

%Норма невязки:
disp('Нормы невязки (1 и inf) СЛАУ №1:');
R1_1_norm = norm(R1, 1);
R1_inf_norm = norm(R1, inf);
disp(R1_1_norm);
disp(R1_inf_norm);

disp('Нормы невязки (1 и inf) СЛАУ №2:');
R2_1_norm = norm(R2, 1);
R2_inf_norm = norm(R2, inf);
disp(R1_1_norm);
disp(R1_inf_norm);

%Абсолютная погрешность норм:
disp('Абсолютная погрешность нормы невязки (1 и inf) СЛАУ №1:');
abs_delta_R1_1_norm = norm(ANS1-X1, 1);
abs_delta_R1_inf_norm = norm(ANS1-X1, inf);
disp(abs_delta_R1_1_norm);
disp(abs_delta_R1_inf_norm);

disp('Абсолютная погрешность нормы невязки (1 и inf) СЛАУ №2:');
abs_delta_R2_1_norm = norm(ANS2-X2, 1);
abs_delta_R2_inf_norm = norm(ANS2-X2, inf);
disp(abs_delta_R2_1_norm);
disp(abs_delta_R2_inf_norm);

%Относительная погрешность норм:
disp('Относительная погрешность нормы невязки (1 и inf) СЛАУ №1:');
rel_delta_R1_1_norm =  norm(ANS1-X1, 1) / norm(ANS1, 1);
rel_delta_R1_inf_norm =  norm(ANS1-X1, inf) / norm(ANS1, inf);
disp(rel_delta_R1_1_norm);
disp(rel_delta_R1_inf_norm);

disp('Относительная погрешность нормы невязки (1 и inf) СЛАУ №2:');
rel_delta_R2_1_norm =  norm(ANS2-X2, 1) / norm(ANS2, 1);
rel_delta_R2_inf_norm =  norm(ANS2-X2, inf) / norm(ANS2, inf);
disp(rel_delta_R2_1_norm);
disp(rel_delta_R2_inf_norm);

%Обратные матрицы:
disp('Матрица, обратная первой матрице:');
A1_inv = Inv_Gauss(A1);
disp(A1_inv);
disp('Матрица, обратная второй матрице:');
A2_inv = Inv_Gauss(A2);
disp(A2_inv);

%Проверка равенства A^-1 * A = E:
disp('Проверка равенства A^-1 * A = E:');
disp(A1_inv * A1);
disp(A2_inv * A2);

%Числа обусловленности матрицы A1:
disp('Числа обусловленности (1 и inf) первой матрицы:');
cond1_1 = cond(A1, 1);
cond1_inf = cond(A1, inf);
disp(cond1_1);
disp(cond1_inf);

%Числа обусловленности матрицы A2:
disp('Числа обусловленности (1 и inf) второй матрицы:');
cond2_1 = cond(A2, 1);
cond2_inf = cond(A2, inf);
disp(cond2_1);
disp(cond2_inf);

% Реализованный метод Гаусса
function x = Gauss(A, B)
    % A - матрица коэффициентов
    % B - вектор свободных членов
    % x - вектор неизвестных
    
    n = size(A, 1);
    
    % прямой ход метода Гаусса
    for k = 1:n-1
        for i = k+1:n
            coef = A(i,k) / A(k,k);
            A(i,k:n) = A(i,k:n) - coef * A(k,k:n);
            B(i) = B(i) - coef * B(k);
        end
    end
    
    % обратный ход метода Гаусса
    x(n) = B(n) / A(n,n);
    for k = n-1:-1:1
        x(k) = (B(k) - A(k,k+1:n)*x(k+1:n)') / A(k,k);
    end
    x = x';
end

% Реализованный метод Гаусса для нахождения обратной матрицы
function A_inv = Inv_Gauss(A)
    n = size(A, 1);

    % находим обратную матрицу методом Гаусса
    E = eye(n);
    for k = 1:n
        A_inv(:,k) = Gauss(A, E(:,k))';
    end

end

