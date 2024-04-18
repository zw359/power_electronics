// Algorithm from paper [Hu2016]
Vo=7;
Po=600;
Rload = Vo^2/Po;

Io = Vo/Rload;
Cr = 16E-9;
Vin_min = 280;
Vin=280;
Fsw_min = 100E3;

N = 16;



Cr = 16E-9;
Lr = 112.5902E-6
Lm = 134.5184E-6;

Fr = 1/(2*%pi*sqrt(Lr*Cr));

Fs = 118.59E3;
//Fs = 130E3;
Fs = 150E3;

F = Fs/Fr;
Lambda = Lr/Lm;
Gamma = %pi/F;

Im_alpha = -Lambda*Gamma/2;  // i_lm(Alpha)
Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;
// 6 unknow variables
// mc(0)        x(1)
// mc(alpha)    x(2)
// iLr(0)       x(3)
// iLm(0)       x(4)
// M            x(5)
// alpha        x(6)

// equation (25a)
func1 =  'res(1) = x(1) + (x(2) -1/x(5) +1) .* cos(Gamma-x(6)) + Im_alpha .* sin(Gamma-x(6)) + 1/x(5)-1';

// equation (22a)
func2 =  'res(2) = -x(2) + (x(1) -1/x(5) -1).*cos(x(6)) + x(3) .* sin(x(6)) + 1/x(5) +1' ;

// equation (22d)
func3 =  'res(3) = x(4)- Lambda * x(6) - Im_alpha';

// equation (25b)
func4 =  'res(4) = x(3) + (-x(2) +1/x(5) -1) .* sin(Gamma-x(6)) + Im_alpha .* cos(Gamma-x(6))';

// equation (22c)
func5 =  'res(5) = -Im_alpha + (-x(1) +1/x(5) +1) .* sin(x(6)) + x(3) .* cos(x(6))';

// euqation (25c)
func6 =  'res(6) = rL * (1/Gamma) * ( x(6) .* x(4) - 0.5* Lambda * x(6).^2 - (-x(1) +1/x(5) +1) .* (1-cos(x(6)))  - x(3).*sin(x(6)) + (-x(2)+1/x(5)-1) .* (1-cos(Gamma-x(6))) +Im_alpha * sin(Gamma-x(6)) - Im_alpha*(Gamma-x(6)) -0.5*Lambda*(Gamma-x(6)).^2) -1' ; 


deff('res=CCMA_mode(x)',[func1; func2; func3; func4; func5; func6]);

// Initial condition
x1_0=-2.144084
x2_0=-2.144474
x3_0=-1.315849
x4_0=-1.314373
x5_0=0.999818
x6_0=0.000296
//Theta=2.924462
//Lambda=0.819139

x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0];
xsol1 =fsolve(x0, CCMA_mode); 
res1 = CCMA_mode(xsol1) ;

mc_0= xsol1(1);
mc_alpha = xsol1(2);
I_lr_0 = xsol1(3);
I_lm_0 = xsol1(4);
M = xsol1(5);
Alpha = xsol1(6);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
Vc_alpha = xsol1(2) * Vbase;
I_lr_0= xsol1(3) * Ibase;
I_lm_0= xsol1(4) * Ibase;

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

