function deriv=dydx(x,y)

numerator=(y(3:end)-y(1:end-2));
denominator=(x(3:end)-x(1:end-2));

deriv=numerator./denominator;