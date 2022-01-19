function [y] = rve(A, tolmax)
y = getMinCond(A) * tolmax;
end
