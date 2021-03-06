% script <a href="matlab:ResponseTFdtTest">ResponseTFdtTest</a>
% Test <a href="matlab:help ResponseTFdt">ResponseTFdt</a> by plotting the impulse, step, and ramp response of an oscillatory system.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 17.1.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap17">Chapter 17</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.

clear; close all; r=.9; theta=pi/10; a0=r^2; a1=-2*r*cos(theta);
numG=1+a1+a0; denG=[1 a1 a0]; g.T=60; g.h=1; g.styleu='r--'; g.styley='b*-'; 
figure, ResponseTFdt(numG,denG,0,g), title('Impulse response') 
figure, ResponseTFdt(numG,denG,1,g), title('Step response')
figure, ResponseTFdt(numG,denG,2,g), title('Ramp response')
figure, ResponseTFdt(numG,denG,3,g), title('Response to quadratic')
figure, ResponseTFdt(numG,denG,4,g), title('Response to cubic')

% end script ResponseTFdtTest