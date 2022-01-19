function condVal = getMinCond(A)
condVal = inf;
m = size(A, 1);
n = size(A, 2);
NN = 20;
M1 = ones(m, n);
M2 = ones(m, n);
for k=1:NN
    EPM = randi([0, 1], m, n);
    for i=1:m
        for j=1:n
            if EPM(i, j) == 0
                M1(i, j) = inf(A(i, j));
                M2(i, j) = sup(A(i, j));
            else
                M2(i, j) = inf(A(i, j));
                M1(i, j) = sup(A(i, j));
            end
        end
    end
    c1 = cond(M1, 2);
    c2 = cond(M2, 2);
    condVal = min([condVal, c1, c2]);
end
end