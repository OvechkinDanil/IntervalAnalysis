% Коэффициенты ИСЛАУ
%A = [infsup(0, 0.5), infsup(1.5, 2); 1, -infsup(3, 3); infsup(-0.5, 0.5), 0];
%b = [infsup(3, 5); 0; infsup(-1, 1)];

A = [infsup(0.5, 1.5), infsup(1.5, 3.5); 1, -infsup(3, 3); infsup(-0.5, 0.5), 0];
b = [infsup(3, 7); 0; infsup(-1, 1)];
% Распознающий функционал
Tol = @(A, b, x) min(rad(b) - mag(mid(b) - A * x));

iveVal = ive(A, b)

rveVal = rve(A, b)

[origtolmax,origargmax,envs,ccode] = tolsolvty(inf(A), sup(A), inf(b), sup(b));
origtolmax
origargmax
n = 100;
levels = n / 3;
x = linspace(-5, 8, n);
y = linspace(-5, 8, n);
[xx, yy] = meshgrid(x, y);
zz = zeros([size(xx, 1), size(xx, 2)]);
for i=1:size(xx, 1)
    for j=1:size(yy,1)
        zz(i, j) = Tol(A, b, [xx(i, j); yy(i, j)]);
    end
end
figure
contour(xx, yy, zz, levels)
hold on
plot(origargmax(1), origargmax(2), 'r*')
hold on 
colorbar
title('Tol($x, \mathbf{A}, \mathbf{b}$)','interpreter','latex')

e = [infsup(-1, 1); infsup(-1, 1); infsup(-1, 1)];
K = 1.5;

b1 = b + K * e;
[btolmax,bargmax,envs,ccode] = tolsolvty(inf(A), sup(A), inf(b1), sup(b1));
btolmax
bargmax
zz = zeros([size(xx, 1), size(xx, 2)]);
for i=1:size(xx, 1)
    for j=1:size(yy,1)
        zz(i, j) = Tol(A, b1, [xx(i, j); yy(i, j)]);
    end
end
figure
contour(xx, yy, zz, levels)
hold on
plot(bargmax(1), bargmax(2), 'r*')
hold on 
colorbar
title('Tol($x, \mathbf{A}, \hat{b}$)','interpreter','latex')

[V,P1,P2,P3,P4] = EqnTolR2(inf(A), sup(A), inf(b1), sup(b1));
iveVal = ive(A, b1)
rveVal = rve(A, b1)

E = 1 * [1 1;0 1;1 0];

A1 = [infsup(0.991, 1.0), infsup(2.498, 2.5); 1, -3; infsup(-0.001, 0.001), 0];

zz = zeros([size(xx, 1), size(xx, 2)]);
for i=1:size(xx, 1)
    for j=1:size(yy,1)
        zz(i, j) = Tol(A1, b, [xx(i, j); yy(i, j)]);
    end
end
[aetolmax,aeargmax,envs,ccode] = tolsolvty(inf(A1), sup(A1), inf(b), sup(b));
aetolmax
aeargmax
figure
contour(xx, yy, zz, levels)
hold on
plot(aeargmax(1), aeargmax(2), 'r*')
hold on 
colorbar
title('Tol($x, \hat{A}, \mathbf{b}$)','interpreter','latex')
[V,P1,P2,P3,P4] = EqnTolR2(inf(A1), sup(A1), inf(b), sup(b));
iveValForMatr = ive(A1, b)
rveValForMatr = rve(A1, b)
