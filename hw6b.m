% The Mobile Inverted Pendulum (MIP) Control Problem, solved in the Successive Loop Closure
% (SLC) framework using classical control techniques in the preferred SLC formulation
% (with theta closed on inner loop, and phi closed on outer loop), with simplified
% constants as suggested in the UCSD MAE143b S2 2013 final exam (version 1)
% Code modified by Paul Vinh Phan from Prof Thomas Bewley, MAE 143C.

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
t1 = abs(c2a*c3-c1^2); t2 = c1 + c3; t2hat = c1 + c2a; % why is t1 negative?
Ts = 0.005; % seconds

% disp('Constants in actual plant:')
k1 = d2*(c1+c3)/t1,
a2 = d1*(2*c1+c2a+c3)/t1;
a1 = c2b*c3/t1;
a0 = d1*c2b/t1;
k2 = (c1+c2a)/t2,
z1 = sqrt(c2b/(t2*k2)),
disp('');
P_G1 = roots([1 a2 -a1 -a0]), Z_G1 = roots([k1 0]), disp('')
% display with proper signsd

% Set up plant
numG1=[k1 0]; denG1=PolyConv([1 -P_G1(1)],[1 -P_G1(2)],[1 -P_G1(3)]);


% Set up various controllers for problem #1:

%%%%%%% Choose parameters %%%%%%%
	disp('Inner Loop Chosen Parameters:')
	tr = 0.18 % Rise Time
	damp_ratio = 0.5 % Damping Ratio	
	disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Results in: ')

	omega_n = 1.8/tr
	Mp = exp(-pi*damp_ratio/sqrt(1-damp_ratio^2))
	omega_d = omega_n*sqrt(1-damp_ratio^2);
	sigma = damp_ratio*omega_n;
	ts = 4.6/sigma 
	pa_initial = (2*sigma + P_G1(1));
	pa = 1*pa_initial % try tweaking
	s_plus = -sigma + omega_d*i;
	a = (s_plus - P_G1(1))*(s_plus + pa); b = k1;
	K1_initial = -real(a/b);

%%%%%%%% Choose gain K1   %%%%%%%
	K1 = 1*K1_initial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('')

numD1 = K1*PolyConv([1 -P_G1(3)],[1 -P_G1(2)]);  denD1 =PolyConv([1 pa], [1 0]);
alpha = 3.5

% Maybe it'll make sense with Pade approx
numD1a= alpha*K1*PolyConv([1 -P_G1(3)],[1 -P_G1(2)]);  denD1a=PolyConv([1 alpha*-P_G1(2)],[1 0]);
numD1b= K1*PolyConv([1 -P_G1(3)],[1 -P_G1(2)]);  denD1b=PolyConv([1 pa], [1 .05]);

disp('Coefficients for difference equation:')
[numD1_DT,denD1_DT] = C2DTustin(numD1,denD1,Ts,omega_n); disp('')
[numD1a_DT,denD1a_DT] = C2DTustin(numD1a,denD1a,Ts,omega_n), disp('')
[numG1_DT,denG1_DT] = C2DTustin(numG1,denG1,Ts,omega_n); disp('')
% Prof's code for above:
	%% numD1 = 3*PolyConv([1 10],[1 1]);  denD1 =PolyConv([1 20], [1 0]);
	%% numD1a=40*PolyConv([1 10],[1 1]);  denD1a=PolyConv([1 100],[1 0]);
	%% numD1b= 3*PolyConv([1 10],[1 1]);  denD1b=PolyConv([1 20], [1 .05]);

% Problem #1a
figure(1); g.K=logspace(-1.5,1.5,400);
% g.axes=[-1.5 1.5 -1.5 1.5]; RLocus(numG1_DT,denG1_DT,numD1a_DT,denD1a_DT,g); % print -depsc figs/final_143b_13_1a.eps
g.axes=[-25 25 -25 25]; RLocus(numG1,denG1,numD1a,denD1a,g); % D1a

% Problem #1b
figure(2); g.omega=logspace(-1.5,2.5,400); g.line=0;
g.style='k-';  Bode(numG1,denG1,g); hold on;
g.style='b-';  Bode(PolyConv(numG1,numD1),PolyConv(denG1,denD1),g);
g.style='r--'; Bode(PolyConv(numG1,numD1a),PolyConv(denG1,denD1a),g);
g.style='m-.'; Bode(PolyConv(numG1,numD1b),PolyConv(denG1,denD1b),g);
subplot(2,1,1); axis([10^(-1.5) 10^(2.5) 10^(-2.5) 10^(1.5)])
subplot(2,1,2); axis([10^(-1.5) 10^(2.5) -180 -90]); % print -depsc figs/final_143b_13_1b.eps

% Problem #1c
numH1 =PolyConv(numG1,numD1);  denH1 =PolyAdd(PolyConv(numG1,numD1), PolyConv(denG1,denD1));
numH1a=PolyConv(numG1,numD1a); denH1a=PolyAdd(PolyConv(numG1,numD1a),PolyConv(denG1,denD1a));
numH1b=PolyConv(numG1,numD1b); denH1b=PolyAdd(PolyConv(numG1,numD1b),PolyConv(denG1,denD1b));

% Problem #1d
P=denH1(end-1)/numH1(end-1);  Pa=denH1a(end-1)/numH1a(end-1),  Pb=P;

% % Problem #1e
g.T=1.2; g.h=Ts; g.styleu='k-'; figure(3)
g.styley='b-';  ResponseTF(P*numH1,  denH1, 1,g); hold on; axis([0 1.2 0 1.4]);
g.styley='r--'; ResponseTF(Pa*numH1a,denH1a,1,g); hold on;
g.styley='m-.'; ResponseTF(Pb*numH1b,denH1b,1,g); % print -depsc figs/final_143b_13_1e.eps

%% Set up various controllers for problem #2: %%

%%%%%%% Choose parameters %%%%%%%
	disp('Outter Loop Chosen Parameters:')
	tr = 1.8 % Rise Time
	damp_ratio = 0.7 % Damping Ratio	
	disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Results in: ')

	omega_n = 1.8/tr
	Mp = exp(-pi*damp_ratio/sqrt(1-damp_ratio^2))
	omega_d = omega_n*sqrt(1-damp_ratio^2);
	sigma = damp_ratio*omega_n;
	ts = 4.6/sigma;
	za = 0.9434; z1 = 16.68977; k2 = 0.5901;
	s_plus = -sigma + omega_d*i;
	quant = (s_plus^2 + s_plus*(za-z1) - z1*za)/s_plus^2;
%%%%%%% Choose gain K2   %%%%%%%
	K2 = real(1/(k2*quant))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('')

numG2 = -k2*PolyConv([1 -z1],[1 z1]); denG2=[1 0 0];
numD2 = K2*[1 za]; denD2 =[1 z1];
alpha = 6,
numD2a = (alpha)*K2*[1 (1/alpha)*z1]; denD2a =[1 z1];

disp('Coefficients for difference equation:'), disp('')
[numD2_DT,denD2_DT] = C2DTustin(numD2,denD2,Ts,omega_n), disp('')
[numD2a_DT,denD2a_DT] = C2DTustin(numD2a,denD2a,Ts,omega_n), disp('')
[numG2_DT,denG2_DT] = C2DTustin(numG2,denG2,Ts,omega_n); disp('')

% Problem #2a - Root Locus Outter Loop
figure(4); g.K=logspace(-1.5,3,800); g.axes=[-2*z1 2*z1 -2*z1 2*z1];
RLocus(numG2,denG2,numD2a,denD2a,g); % print -depsc figs/final_143b_13_2a.eps
% figure(5); g.K=logspace(-1.5,3,800); g.axes=[-1.5 1.5 -1.5 1.5];
% RLocus(numG2_DT,denG2_DT,numD2_DT,denD2_DT,g); % print -depsc figs/final_143b_13_2a.eps

% Problem #2b - Bode Outter Loop
figure(5); g.omega=logspace(-2.5,2.5,400); g.line=0;
g.style='k-';  Bode(numG2,denG2,g); hold on;
g.style='b-';  Bode(PolyConv(numG2,numD2), PolyConv(denG2,denD2),g);
g.style='r--'; Bode(PolyConv(numG2,numD2a),PolyConv(denG2,denD2a),g);
subplot(2,1,1); axis([10^(-2.5) 10^(2.5) 10^(-2) 10^(4)]);
subplot(2,1,2); axis([10^(-2.5) 10^(2.5) -180 -90]); % print -depsc figs/final_143b_13_2b.eps

% Problem #2c
numH2 =PolyConv(numG2,numD2);  denH2 =PolyAdd(PolyConv(numG2,numD2), PolyConv(denG2,denD2));
numH2a=PolyConv(numG2,numD2a); denH2a=PolyAdd(PolyConv(numG2,numD2a),PolyConv(denG2,denD2a));

% Problem #2e
g.T=10; g.h=Ts; g.styleu='k-'; figure(6)
P2=denH2(end-1)/numH2(end-1); P2a=denH2a(end-1)/numH2a(end-1);
g.styley='b-';  ResponseTF(P2*numH2,denH2,1,g)
hold on, g.styley='r--'; ResponseTF(numH2a,denH2a,1,g), % print -depsc figs/final_143b_13_2e.eps

% % Problem #4b
% numH2full=PolyConv(numG2,Pb*numH1b,numD2);
% denH2full=PolyAdd(PolyConv(numG2,Pb*numH1b,numD2),PolyConv(denG2,denH1b,denD2));
% numF1=100^2; denF1=[1 2*.707*100 100^2];
% numF2=10^2;  denF2=[1 2*.707*10  10^2 ];
% numH1bF=PolyConv(numG1,numD1b,denF1);
% denH1bF=PolyAdd(PolyConv(numG1,numD1b,numF1),PolyConv(denG1,denD1b,denF1));
% numH2fullF=PolyConv(numG2,Pb*numH1bF,numD2,denF2);
% denH2fullF=PolyAdd(PolyConv(numG2,Pb*numH1bF,numD2,numF2),PolyConv(denG2,denH1bF,denD2,denF2));
% figure(7), g.styley='b-';  ResponseTF(numH2full, denH2full, 1,g), hold on
%            g.styley='r--'; ResponseTF(numH2fullF,denH2fullF,1,g),
% % print -depsc figs/final_143b_13_4b.eps

% % Problem #5
% g.figs=8; g.eps=.1; g.R=1000;
% g.figL=9; Nyquist(PolyConv(numG1,numD1),PolyConv(denG1,denD1),g); axis equal
% % print -depsc figs/final_143b_13_5a.eps

% g.figL=10; Nyquist(PolyConv(numG1,numD1/3),PolyConv(denG1,denD1),g); axis equal
% % print -depsc figs/final_143b_13_5b.eps

% % Problem #6
% PolyInertia([1 10 100])