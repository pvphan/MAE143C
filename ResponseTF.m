function ResponseTF(numG,denG,type,g)
% function ResponseTF(numG,denG,type,g)
% Using its partial fraction expansion, compute the response Y(s)=G(s)*U(s) of
% a continuous-time SISO linear system to an impulse (type=0), step (type=1),
% or polynomial (type>=2) input.  The derived type g groups together convenient
% plotting parameters: g.T is the interval over which response is plotted,
% g.h is the timestep, and {g.styleu,g.styley} are the linestyles used.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 17.3.3.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap17">Chapter 17</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.
% Depends on <a href="matlab:help PartialFractionExpansion">PartialFractionExpansion</a>, <a href="matlab:help Fac">Fac</a>.  Verify with: <a href="matlab:help ResponseTFtest">ResponseTFtest</a>.

numU=Fac(type-1); denU=1; for i=1:type, denU=[denU 0]; end,  numG=numG/denG(1);
[up,ud,uk]=PartialFractionExpansion(numU,denU);              denG=denG/denG(1);
[yp,yd,yk]=PartialFractionExpansion(PolyConv(numU,numG), PolyConv(denU,denG));
km=1000; h=g.T/km; t=[0:km]*h; for k=1:km+1
  if type>0, u(k)=real(sum(ud.*(t(k).^(uk-1).*exp(up*t(k))))); else, u(k)=0; end
             y(k)=real(sum(yd.*(t(k).^(yk-1).*exp(yp*t(k)))));
end, plot(t,y,g.styley), if type>0, hold on; plot(t,u,g.styleu), hold off; end
end % function ResponseTF