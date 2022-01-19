function [] = drawTolMax(A2,b, line, iterations)
for i = 1:iterations
    C = rad(A2)./2 .*line;
    A2 = infsup(inf(A2)+C, sup(A2)-C);
    [~,argMax,~,~] = tolsolvty(inf(A2), sup(A2), inf(b), sup(b));
    plot(argMax(1), argMax(2), '*r'); 
    hold on    
end
grid on
end