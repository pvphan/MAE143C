function [ ret ] = hw2_3 ()
clear; close all;
% BATTERY AT TOP: 
	% L = 0.065; I_bc = 0.000341;
% BATTERY AT MID:
	L = 0.0586; % m. 		      Distance from wheel axis to MIP body (middle config)
	I_bc = 0.000312; % kg*m^2.    Moment of inertia of MIP body about its CG
% BATTERY AT BOT:
	% L = 0.0489; I_bc = 0.000269; 

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
d1=motor.eff*motor.K^2/motor.R; 
					% {d1,d2} are constants in the motor model  d1*omega(t) + tau(t) = -d2*u(t)
d2=motor.eff*motor.K*motor.V_max/motor.R;   
					% [see (19.11) in NR, taking tau=torque output, u\in[-1,1]=motor driver input, L=0].
                    % Note that omega(t)=d(phi(t)-theta(t))/dot is the speed of the motor on MIP,
                    % where theta is the angle of the MIP body and phi measures the rotation of the wheels.
 clc;
c1 = m_b*R_w*L; c2a = I_bc*m_b*L^2; c2b = m_b*g_const*L; c3 = I_w + R_w^2*(m_b+m_w);
t1 = c2a*c3-c1^2; t2 = c1 + c3; t2hat = c1 + c2a; 
Ts = 0.005; % seconds

% disp('Constants in actual plant:')
k1 = d2*(c1+c3)/t1;
a2 = d1*(2*c1+c2a+c3)/t1;
a1 = c2b*c3/t1;
a0 = d1*c2b/t1;
k2 = (c1+c2a)/t2;
z1 = sqrt(c2b/(t2*k2));
disp('');
d1; d2;

% disp('Poles and Zeros of G1')
P_G1 = roots([1 a2 -a1 -a0]);
Z_G1 = roots([k1 0]);

numG1 = [k1 0];
denG1 = real(PolyConv([1 P_G1(1)],[1 P_G1(2)],[1 P_G1(3)]));
G1 = tf(numG1,denG1)

% disp('Chosen D1 params:')
tr = 0.018 % sec % 0.485! % controls b of difference equation
damp_ratio = 8.66/5 % 8.66/5!
omega_c = 1.8/tr %2! % controls a of difference equation
pa = omega_c - P_G1(3);
s_plus = (omega_c/2)*(-1 + damp_ratio*i);
a = (s_plus + pa)*(s_plus + P_G1(3)); b = k1;
K1 = -(a/b); 
% lead comp
	% pLead = 10*pa; zLead = pa; KLead = pLead/zLead;
	% numD1_lead = [1 zLead]; % ADJUST
	% denD1_lead = [1 pLead]; % ADJUST
	% % pLag = 15; zLag = 100; KLag = 1;
	% % numD1_lag = [1 zLag]; % ADJUST
	% % denD1_lag = [1 pLag]; % ADJUST
	% % numD1_LL = KLead*KLag*PolyConv(numD1_lead, numD1_lag);
	% % denD1_LL = PolyConv(denD1_lead, denD1_lag);
	% numD1 = real(K1*KLead*PolyConv([1 P_G1(1)], [1 P_G1(2)], numD1_lead)); % Pole-zero cancelations with plant G1
	% denD1 = real(PolyConv([1 0], [1 pa], denD1_lead));

numD1 = real(K1*PolyConv([1 P_G1(1)], [1 P_G1(2)])); % Pole-zero cancelations with plant G1
denD1 = real(PolyConv([1 0], [1 pa]));

% CT, DT forms of controler D1
% D1_matlab = tf(numD1,denD1); D1_DT = c2d(D1_matlab,Ts, 'tustin'); [numD1_DT_m,denD1_DT_m] = tfdata(D1_DT,'v'), disp('')
[numD1_DT,denD1_DT] = C2DTustin(numD1,denD1,Ts,omega_c)
% Get coefficients for Difference Equation
% disp('Coefficients b0 ... bk:'), for i=1:length(numD1_DT), D1_b(i) = numD1_DT(length(numD1_DT) - i + 1); D1_b(i), end
% disp('Coefficients a0 ... ak:'),for i=1:length(denD1_DT), D1_a(i) = denD1_DT(length(denD1_DT) - i + 1); D1_a(i), end
% note a0 would be DT_a(1) due to MATLABs index 1 system

disp('')
figure(1), g.K = logspace(-1.5, 1.5, 200); g.axes=[-(3/4)*omega_c (1/4)*omega_c -omega_c/2 omega_c/2 ];
RLocus(numG1,denG1,numD1,denD1,g)

bs = PolyConv(numG1,numD1); 
as = PolyConv(denG1,denD1);
GD = tf(bs,as);

figure(2); clf; g.omega=logspace(-1,1,500); g.line=0; 
g.style='k-'; Bode(bs,as,g), hold on, h=Ts;
g.style='b-'; [bz,az]=C2DTustin(numD1,denD1,h);   Bode(bz,az,g,h)
g.style='r-'; [bz,az]=C2DTustin(bs,as,h,omega_c); Bode(bz,az,g,h)

numH1 = real(PolyConv(numG1,numD1));
denH1 = real(PolyAdd(PolyConv(denG1,denD1),PolyConv(numG1,numD1)));

% Get Prefactor of H1
i = 0; P = 0;
for i=0:length(numH1)-1
	if (numH1(length(numH1)-i) ~= 0) && (denH1(length(denH1)-1) ~= 0)
		P = denH1(length(denH1)-i)/numH1(length(numH1)-i);
		break;
	end
end
H1 = tf(numH1,denH1);
PH1 = P*H1;
PH1_DT = c2d(PH1,Ts,'tustin', omega_c);
figure(3), step(PH1, 'r', PH1_DT, 'b');

% Part 2

end %hw2_3