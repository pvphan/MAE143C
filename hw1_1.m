function [ ret ] = hw1_1 ()
clear; clc; close all;
% L = 0.065; % m battery at top
% I_bc = 0.000341; % kg*m^2 battery at top
L = 0.0586; % m. 		      Distance from wheel axis to MIP body (middle config)
I_bc = 0.000312; % kg*m^2.    Moment of inertia of MIP body about its CG
% L = 0.0489; % m battery at bottom
% I_bc = 0.000269 % kg*m^2 battery at bottom
I_bw = 0.00116; % kg*m^2.     Moment of inertia of MIP body about wheel axis
m_b = 0.249; % kg. 			  Total mass of BeagleMIP body
m_w = 0.0262; % kg. 		  Total mass of both wheels
R_w = 0.0352; % m. 			  Radius of wheels
I_w = m_w*(R_w^2)/2; % kg*m^2.Moment of both wheels about thier CM
g_const = 9.81; % m/s^2

% Motor Parameters
motor.K=0.0525;     % Motor torque constant [measured in N*m/A = V/(rad/s)] for the 50:1 Pololu Motors
motor.R=7.22;       % Motor resistance [measured in ohms]
motor.V_max=7.9;    % Power applied to the motor driver / H-bridge [measured in volts]
motor.eff = 0.9     % Motor/ gearbox efficiency
d1=motor.eff*motor.K^2/motor.R; % {d1,d2} are constants in the motor model  d1*omega(t) + tau(t) = -d2*u(t)
d2=motor.eff*motor.K*motor.V_max/motor.R;   % [see (19.11) in NR, taking tau=torque output, u\in[-1,1]=motor driver input, L=0].
                    % Note that omega(t)=d(phi(t)-theta(t))/dot is the speed of the motor on MIP,
                    % where theta is the angle of the MIP body and phi measures the rotation of the wheels.

c1 = m_b*R_w*L; c2a = I_bc*m_b*L^2; c2b = m_b*g_const*L; c3 = I_w + R_w^2*(m_b+m_w);
t1 = c2a*c3-c1^2; t2 = c1 + c3; t2hat = c1 + c2a; 

disp('Constants in actual plant:')
k1 = d2*(c1+c3)/t1
a2 = d1*(2*c1+c2a+c3)/t1
a1 = c2b*c3/t1
a0 = d1*c2b/t1
k2 = (c1+c2a)/t2
z1 = sqrt(c2b/(t2*k2))

disp('Poles and Zeros of G1')
P_G1 = roots([1 a2 -a1 a0])
Z_G1 = roots([k1 0])

numG1 = [k1 0];
denG1 = PolyConv([1 P_G1(1)],[1 P_G1(2)],[1 P_G1(3)]);

disp('Chosen D1 params:')
pa = 10 - P_G1(3)
s_plus = -5 + 8.66*i
as = (s_plus + pa)*(s_plus + P_G1(3)); bs = k1;
K1 = -(as/bs)
numD1 = K1*PolyConv([1 P_G1(1)], [1 P_G1(2)]); % Pole-zero cancelations with plant G1
denD1 = PolyConv([1 0], [1 pa]);

figure(1), g.K = logspace(-1.5, 1.5, 200); g.axes=[-16 16 -16 16];
RLocus(numG1,denG1,numD1,denD1,g)

figure(2), g.omega=logspace(-1,1,500); g.line=0;
g.style='k-'; Bode(bs,as,g), hold on, h=1;
g.style='b-'; [bz,az]=C2DTustin(bs,as,h);   Bode(bz,az,g,h)
g.style='r-'; [bz,az]=C2DTustin(bs,as,h,1); Bode(bz,az,g,h)
end %hw1_1