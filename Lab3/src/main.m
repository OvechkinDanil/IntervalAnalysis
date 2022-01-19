% Коэффициенты ИСЛАУ
A = [infsup(0.5, 1.5), infsup(1.5, 3.5); 1, -3; infsup(-0.5, 0.5), 0];
b = [infsup(3, 7); 0; infsup(-1, 1)];
% Распознающий функционал
Tol = @(A,b,x) min(rad(b) - mag(mid(b) - A * x));

[tolMax,argMax,envs,ccode] = tolsolvty(inf(A), sup(A), inf(b), sup(b));
tolMax
argMax
%  Допусковое множество решений интервальной линейной системы пусто 


n = 100;
levels = 40;
drawContour(A,b,n,levels)

%% Коррекция вектора b
e = [infsup(-1, 1); infsup(-1, 1); infsup(-1, 1)];
C = 1.5 * abs(tolMax);
b1 = b + C * e;
b1
drawContour(A,b1,n,levels)

[tolMax1,argMax1,~,~] = tolsolvty(inf(A), sup(A), inf(b1), sup(b1));
tolMax1
argMax1
ive1 = ive(A, b1);
rve1 = rve(A, tolMax1);
iveBox = [midrad(argMax1(1), ive1);midrad(argMax1(2), ive1)];
rveBox = [midrad(argMax1(1), rve1);midrad(argMax1(2), rve1)];
plotintval([iveBox, rveBox], 'n');

%% Коррекция матрицы А
koef = 1.5;

%% Получил её из python-скрипта
A_new = [infsup(0.991, 1.0), infsup(2.498, 2.5); 1, -3; infsup(-0.001, 0.001), 0];
drawContour(A_new,b,n,levels);

[tolMax2,argMax2,~,~] = tolsolvty(inf(A_new), sup(A_new), inf(A_new), sup(A_new));
tolMax2
argMax2
ive2 = ive(A_new, b);
rve2 = rve(A_new, tolMax2);
iveBox2 = [midrad(argMax2(1), ive2);midrad(argMax2(2), ive2)];
rveBox2 = [midrad(argMax2(1), rve2);midrad(argMax2(2), rve2)];
plotintval([iveBox2, rveBox2], 'n');


%% Положение максимума Tol
iterations = 10;
figure
A2 = A;
for i = 1:iterations
    A2 = A2 ./ 2;
    [~,argMax,~,~] = tolsolvty(inf(A2), sup(A2), inf(b), sup(b));
    plot(argMax(1), argMax(2), '*b');
    hold on
end
grid on

%Положение максимумов Tol
A2 = A;
line1 = [1, 1; 0, 0; 0, 0]
line2 = [0, 0; 1, 1; 0, 0]
line3 = [0, 0; 0, 0; 1, 1]
figure
drawTolMax(A2, b, line1, iterations)
figure
drawTolMax(A2, b, line2, iterations)
figure
drawTolMax(A2, b, line3, iterations)

