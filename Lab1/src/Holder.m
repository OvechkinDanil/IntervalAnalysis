function z = Himmel(x, y)
z = -abs(sin(x).*cos(y).*exp(abs(1-(sqrt(x.^2 + y.^2) / pi))));
end