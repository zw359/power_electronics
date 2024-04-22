// Algorithm from paper [Yu2012]
// DCMAB mode solver  same as OPO mode
// This mode happens both below and above the Fr
// @ Very light load
function [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB_OPO_solver(Cr, Lr, Lm, N, Rload,Fn)
/*
Vin=280;
Rload =5;
N = 16;

Cr = 25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;
*/
Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 
Fr = 1/(2*%pi*sqrt(Lr*Cr));

Fs = 100E3;
Fs = 90E3;
//Fs = Fr;

//F = Fs/Fr;
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

// Initial condition
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
x13_0=2 // M > 1 
x14_0=Gamma *0.3 // critial alpha
x15_0=Gamma *0.6  // critical beta  x15_0 > x14_0


x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0; x7_0; x8_0; x9_0; x10_0; x11_0; x12_0; x13_0; x14_0; x15_0];
[xsol1, res1, info] =fsolve(x0, DCMAB_mode); 

mc_0= xsol1(1);
mc_alpha = xsol1(2);
I_lr_0 = xsol1(5);
I_lr_alpha = xsol1(6);
M = xsol1(13);
iLr_0 = xsol1(5);
iLm_0 =xsol1(9);
alpha = xsol1(14);
beta1 =xsol1(15);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
Vc_alpha = xsol1(2) * Vbase;
I_lr_0= xsol1(3) * Ibase;
I_lr_alpha= xsol1(4) * Ibase;

mm2_0 = abs((-mc_0 + 1/M )/(1+Lambda));
mm2_alpha = abs((-mc_alpha +1/M)/(1+Lambda));
mm2_gamma = abs((-xsol1(4) +1/M)/(1+Lambda));
mm2_beta = abs((-xsol1(3) +1/M)/(1+Lambda));

dt = 0.5;
// iLr(beta)   P mode alpha -- beta
I_lr_alpha_d =  (-xsol1(2) +1/xsol1(13) -1 ) .* sin((xsol1(15)-xsol1(14))*dt) + xsol1(6) .*cos((xsol1(15)-xsol1(14))*dt) ;
// iLm(beta)   P mode alpha -- beta
I_lm_alpha_d =  xsol1(10) + Lambda * ((xsol1(15)-xsol1(14))*dt);

DCMAB_OPO =%F;
if (I_lr_alpha_d > I_lm_alpha_d ) & (xsol1(15)<Gamma) & (mm2_0<1) & (mm2_gamma<1) &( xsol1(14)< xsol1(15) ) & (xsol1(14)>0) & (xsol1(15)>0) & /*(info ==1) &*/ (xsol1(13)>0 )then
    DCMAB_OPO = %T;
end
correctMode = DCMAB_OPO;

/*
printf('DCMAB_OPO mode check %s \n', DCMAB_OPO)
printf('|mm2_0| \t |mm2_alpha| \t |mm2_gamma|\n')
printf('%f \t %f \t %f\n', abs(mm2_0), abs(mm2_alpha), abs(mm2_gamma));

printf('x1_0=%f\n', xsol1(1));
printf('x2_0=%f\n', xsol1(2));
printf('x3_0=%f\n', xsol1(3));
printf('x4_0=%f\n', xsol1(4));
printf('x5_0=%f\n', xsol1(5));
printf('x6_0=%f\n', xsol1(6));

printf('res1=%f\n', res1);
printf('Vo = %f\n', Vo);
*/
endfunction
