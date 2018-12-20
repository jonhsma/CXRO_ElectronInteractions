function yq = lineQuery(x,y,xq)
    yq = (xq-x(1))./(x(2)-x(1)).*(y(2)-y(1))+y(1);
end