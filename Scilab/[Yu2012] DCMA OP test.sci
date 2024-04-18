// Algorithm from paper [Yu2012]
// DCMB2 mode solver  same as PON mode

Vin=280;
Rload = 5;
N = 16;

Cr = 25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;

Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 

Fr = 1/(2*%pi*sqrt(Lr*Cr));

Fs = Fr;
Fs = 180E3;
//Fs = Fr;

F = Fs/Fr;
Lambda = Lr/Lm;
Gamma = %pi/F;

Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;

// general variables  TWO sections
// use the direct equations
// mc(0)        x(1)
// mc(alpha)    x(2)
// mc(gamma)    x(3)
// iLr(0)       x(4)
// iLr(alpha)   x(5)
// iLr(gamma)   x(6)
// iLm(0)       x(7)
// iLm(alpha)   x(8)
// iLm(gamma)   x(9)
// M            x(10)
// alpha        x(11)

func1 = 'res(1) = -x(2) + (x(1)-1/x(10)+1).*cos(x(11)) + x(4).*sin(x(11)) +1/x(10) -1';

//  iLr(alpha)  P mode  0-- alpha
func2 = 'res(2) = -x(5)  +  (-x(1) +1/x(10) -1 ) .* sin(x(11)) + x(4) .*cos(x(11))';

//  iLm(alpha)  P mode  0-- alpha
func3 = 'res(3) = -x(8) + x(7) + Lambda * x(11)';

//  mc(gamma)    N mode alpha -- gamma
func4 = 'res(4) =  -x(3) + (x(2) -1/x(10) -1) .* cos(Gamma-x(11))  + x(5) * sin(Gamma-x(11)) + 1/x(10) +1' ;

//  iLr(gamma)  N mode alpha -- gamma
func5 = 'res(5) = -x(6) +  (-x(2) +1/x(10) +1) .* sin(Gamma-x(11))  + x(5) * cos(Gamma-x(11))'; 

//  iLm(gamma)   N mode alpha -- gamma
func6 = 'res(6) = -x(9) + x(8) - Lambda * (Gamma-x(11))';


//  mc(0) + mc(gamma) =0
func7 = 'res(7) = x(1)+x(3)';

// iLr(0) + iLr(gamma) =0
func8 = 'res(8) = x(4) +x(6)';

// iLm(0) + iLm(gamma) =0
func9 = 'res(9) = x(7) + x(9)';

// iLr(alpha) - iLm(alpha) =0
func10 ='res(10) =  x(5) -x(8)' ;

// iout  voltage balance balance
func11 ='res(11) =  (( -x(1)+1/x(10)-1 ) .* (1-cos(x(11)))  + x(4) .*sin(x(11)) -x(7).*x(11) -0.5*Lambda*x(11).^2  + x(8) .*(Gamma-x(11)) - 0.5*Lambda*( Gamma -x(11)).^2 + ( -x(2) + 1/x(10) +1 ) .* ( cos(Gamma - x(11))-1) - x(5) .* sin(Gamma -x(11))) * rL - Gamma'; 


deff('res=DCMA2_mode(x)',[func1; func2; func3; func4; func5; func6; func7; func8; func9; func10; func11]);


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
x10_0=1
x11_0=%pi


x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0; x7_0; x8_0; x9_0; x10_0; x11_0];
xsol1 =fsolve(x0, DCMA2_mode); 
res1 = DCMA2_mode(xsol1) ;

mc_0= xsol1(1);
mc_alpha = xsol1(2);
mc_gamma = xsol1(3);
I_lr_0 = xsol1(4);
I_lr_alpha = xsol1(5);
I_lr_gamma = xsol1(6);
I_lm_0 =xsol1(7);
I_lm_alpha =xsol1(8)
I_lm_gamma =xsol1(9)
M = xsol1(10);
Alpha = xsol1(11);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = mc_0 * Vbase;
Vc_alpha = mc_alpha * Vbase;
I_lr_0= I_lr_0 * Ibase;
I_lr_alpha= I_lr_alpha * Ibase;

mm2_0 = (-mc_0 + 1/M )/(1+Lambda);
mm2_alpha = (-mc_alpha +1/M)/(1+Lambda);

mm2_gamma = (-mc_gamma - 1/M)/(1+Lambda);


printf('CCM mode check \n')
printf('|mm2_0| \t |mm2_alpha| \t |mm2_gamma|\n')
printf('%f \t %f \t %f\n', abs(mm2_0), abs(mm2_alpha), abs(mm2_gamma));

printf('x1_0=%f\n', xsol1(1));
printf('x2_0=%f\n', xsol1(2));
printf('x3_0=%f\n', xsol1(3));
printf('x4_0=%f\n', xsol1(4));
printf('x5_0=%f\n', xsol1(5));
printf('x6_0=%f\n', xsol1(6));
printf('x7_0=%f\n', xsol1(7));
printf('x8_0=%f\n', xsol1(8));
printf('x9_0=%f\n', xsol1(9));
printf('x10_0=%f\n', xsol1(10));
printf('x11_0=%f\n', xsol1(11));
printf('x12_0=%f\n', xsol1(12));
printf('x13_0=%f\n', xsol1(13));
printf('x14_0=%f\n', xsol1(14));
printf('x5_0=%f\n', xsol1(15));

printf('res1=%f\n', res1);
printf('Vo = %f\n', Vo);

