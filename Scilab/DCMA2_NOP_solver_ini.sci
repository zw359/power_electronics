// Algorithm from paper [Yu2012]
// DCMA2 mode solver  same as NOP mode
// this a very narrow working mode
// set Rload at 2.3
// mm2_alpha <1
// mm2_beta >=1
// I_lr_0 < I_lm_0
function [x1, correctMode]= DCMA2_NOP_solver_ini(Cr, Lr, Lm, N, Rload,Fn, x0)
Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 
Fr = 1/(2*%pi*sqrt(Lr*Cr));
Lambda = Lr/Lm;
Gamma = %pi/Fn;
Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;

// general variables for 3-section operation
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
// Beta         x(15)

//  mc(alpha)   N mode  0-- alpha
func1 = 'res(1) = -x(2) + (x(1)-1/x(13)-1).*cos(x(14)) + x(5).*sin(x(14)) +1/x(13) +1';
//  iLr(alpha)  N mode  0-- alpha
func2 = 'res(2) = -x(6) + (-x(1) +1/x(13) +1 ) .* sin(x(14)) + x(5) .*cos(x(14))';
//  iLm(alpha)  N mode  0-- alpha
func3 = 'res(3) = -x(10) + x(9) - Lambda * x(14)';
//  mc(beta)    O mode alpha -- beta
func4 = 'res(4) =  -x(3) + (x(2)-1/x(13)) .* cos(Kx*(x(15)-x(14)))  + (x(6)/Kx) .* sin(Kx*(x(15)-x(14))) + 1/x(13)';
// iLr(beta)   O mode alpha -- beta
func5 = 'res(5) =  -x(7) + (-x(2) +1/x(13)) *Kx .* sin(Kx*(x(15)-x(14)))  + x(6) .* cos(Kx*(x(15)-x(14)))' ;
// iLm(beta)   O mode alpha -- beta
func6 =  'res(6) = -x(11) + x(7)'; 
//  mc(gamma)    P mode beta -- gamma
func7 = 'res(7) =  -x(4) + (x(3) -1/x(13) +1) .* cos(Gamma-x(15))  + x(7) * sin(Gamma-x(15)) + 1/x(13) -1' ;
//  iLr(gamma)  P mode beta -- gamma
func8 = 'res(8) = -x(8) + (-x(3) +1/x(13) -1) .* sin(Gamma-x(15))  + x(7) * cos(Gamma-x(15))'; 
//  iLm(gamma)   P mode.  beta -- gamma
func9 = 'res(9) = -x(12)  + x(11) + Lambda * (Gamma-x(15))';
//  mc(0) + mc(gamma) =0
func10 = 'res(10) = x(1)+x(4)';
// iLr(0) + iLr(gamma) =0
func11 = 'res(11) = x(5) +x(8)';
// iLm(0) + iLm(gamma) =0
func12 = 'res(12) = x(9) + x(12)';
// iLr(alpha) - iLm(alpha) =0
func13 =' res(13) =  x(6) -x(10)' ;
//  mm(beta) = 1
//  mc(beta) = 1/M - Lambda -1
func14 = 'res(14) = -x(3)+ 1/x(13) - Lambda - 1';
//equation iout balance
func15 = 'res(15) = ( x(9).*(x(14))- 0.5*Lambda*x(14).^2 - (-x(1)+1/x(13)+1).*(1-cos(x(14))) - x(5).* sin(x(14))  + (-x(3)+1/x(13)-1).* (1-cos(Gamma-x(15))) + x(7).*sin(Gamma-x(15)) - x(11) .* (Gamma-x(15)) -0.5 * Lambda *(Gamma - x(15)).^2    ) * rL -Gamma'; 
deff('res=DCMA2_NOP_mode(x)',[func1; func2; func3; func4; func5; func6; func7; func8; func9; func10; func11; func12; func13; func14; func15]);

/* Initial conditions disabled here. they will come as input variables
x1_0=0
x2_0=0
x3_0=0
x4_0=0
x5_0=0
x6_0=0
x7_0=0
x8_0=0
x9_0=0
x10_0=0
x11_0=0
x12_0=0
x13_0=1
x14_0=0
x15_0=0
x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0; x7_0; x8_0; x9_0; x10_0; x11_0; x12_0; x13_0; x14_0; x15_0];
*/
[xsol1, res1,info] =fsolve(x0, DCMA2_NOP_mode); 

mc_0= xsol1(1);
mc_alpha = xsol1(2);
mc_beta = xsol1(3);
mc_gamma = xsol1(4);
I_lr_0 = xsol1(5);
I_lr_alpha = xsol1(6);
I_lr_beta = xsol1(7);
I_lr_gamma = xsol1(8);
I_lm_0 = xsol1(9);
I_lm_alpha = xsol1(10);
I_lm_beta = xsol1(11);
I_lm_gamma = xsol1(12);
M = xsol1(13);
alpha = xsol1(14);
Beta  = xsol1(15);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
Vc_alpha = xsol1(2) * Vbase;

// ====== detect the operation modes ===========
//  "O" mode -->  mm2 < 1,  
//  nomalized output voltage m2  = 1
mm2_0 = abs((-mc_0 + 1/M )/(1+Lambda));
mm2_alpha = abs((-mc_alpha +1/M)/(1+Lambda)); //"O" mode
mm2_beta = abs((-mc_beta +1/M)/(1+Lambda));
mm2_gamma = abs((-mc_gamma +1/M)/(1+Lambda));

// the above methods [Yu2012] don't work reliably
// use the methods [Wei2021]
// take the point half way between 0 and alpha
dt=0.5;
//  iLr(alpha)  N mode  0-- alpha
I_lr_0_d = (-xsol1(1) +1/xsol1(13) +1 ) .* sin(xsol1(14)*dt) + xsol1(5) .*cos(xsol1(14)*dt);
//  iLm(alpha)  N mode  0-- alpha
I_lm_0_d = xsol1(9) - Lambda * xsol1(14)*dt;

DCMA2_NOP =%F;
if (I_lr_0_d < I_lm_0_d ) & (alpha>0) & (Beta>alpha) &  (Beta<Gamma) & (mm2_alpha <1) then
    DCMA2_NOP = %T;
end

correctMode = DCMA2_NOP;
x1=xsol1;

endfunction
