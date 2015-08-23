function ResponseTFdt(numG,denG,type,g)
% function ResponseTFdt(numG,denG,type,g)
% Using its partial fraction expansion, compute the response Y(z)=G(z)*U(z) of
% a discrete-time SISO linear system to an impulse (type=0), step (type=1),
% or polynomial (type>=2) input.  The derived type g groups together convenient
% plotting parameters: g.T is the interval over which response is plotted,
% g.h is the timestep, and {g.styleu,g.styley} are the linestyles used.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 17.4.3.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap17">Chapter 17</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.
% Depends on <a href="matlab:help PartialFractionExpansion">PartialFractionExpansion</a>, <a href="matlab:help Fac">Fac</a>, <a href="matlab:help PolylogarithmNegativeInverse">PolylogarithmNegativeInverse</a>.
% Verify with: <a href="matlab:help ResponseTFdtTest">ResponseTFdtTest</a>.

switch type, case 0, numU=1; denU=1; case 1, numU=[1 0]; denU=[1 -1]; 
             otherwise, [numU,denU]=PolylogarithmNegativeInverse(type-1,1); end
[ua,ud,upp,un]=PartialFractionExpansion(numU,denU);
[ya,yd,ypp,yn]=PartialFractionExpansion(PolyConv(numU,numG),PolyConv(denU,denG));
k=[0:g.T/g.h]; t=k*g.h; y=zeros(size(k)); u=zeros(size(k));
for i=1:yn, a=yd(i)/(Fac(ypp(i)-1)*ya(i)^ypp(i));
  b=a*ones(size(k)); for j=1:ypp(i)-1, b=b.*(k-j); end
  if ypp(i)>0, y(2:end)=y(2:end)+b(2:end).*ya(i).^k(2:end); else, y(1)=y(1)+yd(i); end
end
for i=1:un, a=ud(i)/(Fac(upp(i)-1)*ua(i)^upp(i));
  b=a*ones(size(k)); for j=1:upp(i)-1, b=b.*(k-j); end
  if upp(i)>0, u(2:end)=u(2:end)+b(2:end).*ua(i).^k(2:end); else, u(1)=u(1)+ud(i); end
end, plot(t,real(y),g.styley), if type>0, hold on; plot(t,real(u),g.styleu), hold off; end
end % function ResponseTFdt