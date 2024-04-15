// Algorithm from paper [Yu2012]
// DCMB mode solver
// NOT work yet
Vin=280;
Rload = 0.7;
N = 16;

Cr = 16E-9;
Lr = 112.5902E-6
Lm = 134.5184E-6;



Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 

Fr = 1/(2*%pi*sqrt(Lr*Cr));

Fs = Fr;
Fs = 0.6812*Fr;
Fs = 100E3;
//Fs = Fr;

F = Fs/Fr;
Lambda = Lr/Lm;
Gamma = %pi/F;

Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;
// 6 unknow variables
// mc(0)        x(1)
// mc(alpha)    x(2)
// iLr(0)       x(3)
// iLr(alpha)   x(4)
// M            x(5)
// alpha        x(6)

// equation mc(alpha)   P mode  0-- alpha
func1 = 'res(1) = -x(2) + (x(1)-1/x(5)+1).*cos(x(6)) + x(3).*sin(x(6)) +1/x(5) -1';

// equation iLr(alpha)  P mode  0-- alpha
func2 = 'res(2) = -x(4) + ( -x(1) +1/x(5) -1 ) .* sin(x(6)) + x(3) .*cos(x(6))';

// equation iLm(alpha)  P mode  0-- alpha
func3 = 'res(3) = -x(4) + x(3) + Lambda * x(6)';

// equation (18a)
func4 = 'res(4) = x(1) + (x(2) -1/x(5)) .* cos(Kx *(Gamma-x(6))) + (1/Kx) * x(4) .* sin(Kx*(Gamma-x(6))) + 1/x(5)';

// equation (18b)
func5 = 'res(5) = x(3) +  Kx *(-x(2) +1/x(5)).* sin(Kx*(Gamma -x(6))) + x(4) .* cos(Kx*(Gamma -x(6)))' ;

// equation (18d)
func6 = 'res(6) =  (rL/Gamma) * ((-x(1) + 1/x(5) -1) .* (1- cos(x(6)))  + x(3) .* sin(x(6)) -x(3) .* x(6) - 0.5 * Lambda * x(6).^2) -1';




deff('res=DCMB_mode(x)',[func1; func2; func3; func4; func5; func6]);

// Initial condition
x1_0=0
x2_0=0
x3_0=0
x4_0=0
x5_0=1
x6_0=%pi

x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0];
xsol1 =fsolve(x0, DCMB_mode); 
res1 = DCMB_mode(xsol1) ;

mc_0= xsol1(1);
mc_alpha = xsol1(2);
I_lr_0 = xsol1(3);
I_lr_alpha = xsol1(4);
M = xsol1(5);
Alpha = xsol1(6);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
Vc_alpha = xsol1(2) * Vbase;
I_lr_0= xsol1(3) * Ibase;
I_lr_alpha= xsol1(4) * Ibase;

mm2_0 = (-mc_0 + 1/M )/(1+Lambda);
mm2_alpha = (-mc_alpha +1/M)/(1+Lambda);
mm2_gamma = (mc_0 -1/M)/(1+Lambda);

printf('CCM mode check \n')
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

