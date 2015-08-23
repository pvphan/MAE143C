clc; format long g;

za = 0.9434; z1 = 16.68977; k2 = 0.5901;
s_plus = -0.5 + 0.866*i;
quant = (s_plus^2 + s_plus*(za-z1) - z1*za)/s_plus^2;
K2_plus = 1/(k2*quant)

s_plus = -0.5 - 0.866*i;
quant = (s_plus^2 + s_plus*(za-z1) - z1*za)/s_plus^2;
K2_minus = 1/(k2*quant)
