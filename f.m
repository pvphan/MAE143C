clc;

function y = f (x)
	x1 = sym('x1');
	x2 = sym('x2');
	y1 = -2*x1^2 + 3*x1*x2   + 4*sin(x2) - 6;
	y2 =  3*x1^2 - 2*x1*x2^2 + 3*cos(x1) + 4;

	[x, fval, info] = fsolve (@f, [1; 2])
end