function [y] = ive(A, b)
[tolmax,argmax,~,~] = tolsolvty(inf(A), sup(A), inf(b), sup(b));
y = getMinCond(A) * norm(argmax) * tolmax / norm(mag(b));
end
