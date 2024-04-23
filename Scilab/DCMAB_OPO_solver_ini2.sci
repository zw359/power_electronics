// Algorithm from paper [Yu2012]
// DCMAB mode solver  same as OPO mode
// This mode happens both below and above the Fr
// @ Very light load
function [x1, correctMode] = DCMAB_OPO_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0)
Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 
Fr = 1/(2*%pi*sqrt(Lr*Cr));
Lambda = Lr/Lm;
Gamma = %pi/Fn;

Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;

// general variables 
// use the direct equations
// mc(0)        x(1)
// mc(alpha)    x(2)
// mc(beta)     x(3)
// mc(gamma)    x(4)
// iLr(0)       x(5)
// iLr(alpha)   x(6)
// iLr(beta)    x(7)
// iLr(gamma)   x(8)
// iLm(0)       x(9)
// iLm(alpha)   x(10)
// iLm(beta)    x(11)
// iLm(gamma)   x(12)
// M            x(13)
// alpha        x(14)
// beta         x(15)

//  mc(alpha)   O mode  0-- alpha
func1 = 'res(1) = -x(2) + (x(1)-1/x(13)) .* cos(Kx*x(14))  + (x(5)/Kx) .* sin(Kx*x(14)) + 1/x(13)';
//  iLr(alpha)  O mode  0-- alpha
func2 = 'res(2) = -x(6)  + (-x(1) +1/x(13)) *Kx .* sin(Kx*x(14))  + x(5) .* cos(Kx*x(14))';
//  iLm(alpha)  O mode  0-- alpha
func3 = 'res(3) = -x(10) + x(6)';
//  mc(beta)    P mode alpha -- beta
func4 = 'res(4) =  -x(3) + (x(2)-1/x(13)+1) .* cos(x(15)-x(14)) + x(6).*sin(x(15)-x(14)) +1/x(13) -1';
// iLr(beta)   P mode alpha -- beta
func5 = 'res(5) =  -x(7) + (-x(2) +1/x(13) -1 ) .* sin(x(15)-x(14)) + x(6) .*cos(x(15)-x(14))' ;
// iLm(beta)   P mode alpha -- beta
func6 = 'res(6) =-x(11) + x(10) + Lambda * (x(15)-x(14))';
//  mc(gamma)    O mode beta -- gamma
func7 = 'res(7) =  -x(4) + (x(3)-1/x(13)) .* cos(Kx*(Gamma-x(15)))  + (x(7)/Kx) .* sin(Kx*(Gamma-x(15))) + 1/x(13)' ;
//  iLr(gamma)  O mode beta -- gamma
func8 = 'res(8) = -x(8) + (-x(3) +1/x(13))*Kx .* sin(Kx*(Gamma-x(15)))  + x(7) * cos(Kx*(Gamma-x(15)))'; 
//  iLm(gamma)   O mode.  beta -- gamma
func9 = 'res(9) = -x(12) + x(8)';
//  mc(0) + mc(gamma) =0
func10 = 'res(10) = x(1)+x(4)';
// iLr(0) + iLr(gamma) =0
func11 = 'res(11) = x(5) +x(8)';
// iLm(0) + iLm(gamma) =0
func12 = 'res(12) = x(9) + x(12)';
// iLr(beta) - iLm(beta) =0
func13 ='res(13) =  x(7) -x(11)' ;
//  mm(alpha) =1
//  mc(alpha) = 1/M - Lambda -1
func14 = 'res(14) = -x(2)+ 1/x(13) - Lambda - 1';
// iout  voltage balance balance
func15 = 'res(15) = ( (-x(2) +1/x(13) -1) .*( 1-cos(x(15)-x(14))) + x(6) .* sin(x(15) - x(14)) -x(10).*(x(15)-x(14)) - 0.5*Lambda*(x(15)-x(14)).^2 ) * rL -Gamma'; 

deff('res=DCMAB_mode(x)',[func1; func2; func3; func4; func5; func6; func7; func8; func9; func10; func11; func12; func13; func14; func15]);

[xsol1, res1, info] =fsolve(x0, DCMAB_mode); 

mm2_0 = abs((-xsol1(1) + 1/xsol1(13) )/(1+Lambda));
mm2_alpha = abs((-xsol1(2) +1/xsol1(13))/(1+Lambda));
mm2_gamma = abs((-xsol1(4) +1/xsol1(13))/(1+Lambda));
mm2_beta = abs((-xsol1(3) +1/xsol1(13))/(1+Lambda));

dt = 0.5;
// iLr(beta)   P mode alpha -- beta
I_lr_alpha_d =  (-xsol1(2) +1/xsol1(13) -1 ) .* sin((xsol1(15)-xsol1(14))*dt) + xsol1(6) .*cos((xsol1(15)-xsol1(14))*dt) ;
// iLm(beta)   P mode alpha -- beta
I_lm_alpha_d =  xsol1(10) + Lambda * ((xsol1(15)-xsol1(14))*dt);

DCMAB_OPO =%F;
if (I_lr_alpha_d > I_lm_alpha_d ) & (xsol1(15)<Gamma) & (mm2_0<1) & (mm2_gamma<1) &( xsol1(14)< xsol1(15) ) & (xsol1(14)>0) & (xsol1(15)>0) & (info ==1) & (xsol1(13)>0 )then
    DCMAB_OPO = %T;
end
correctMode = DCMAB_OPO;
x1=xsol1;
endfunction
