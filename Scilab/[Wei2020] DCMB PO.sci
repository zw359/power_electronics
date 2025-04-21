// Algorithm from paper [Yu2012]
// Different design method from [Wei2020]]
// DCMB mode solver
// NOT work yet
Vin_min = 280;
Vin_max = 360;
Vo = 48; // fixed Vo
Vo_min = 46;
Vo_max = 60;
Po = 192;
Rload = Vo^2/Po;

// this could be designed differently
// use Vin nominal
N = Vin_max/(2*Vo);

//Omega0= 1/sqrt(Lr*Cr);
//Omega1= 1/sqrt((Lr+Lm)*Cr);
//Kx= Omega1/Omega0; 

Cr = 60E-9;

M_max = 2*N*Vo/Vin_min;
Fr = 1/(2*%pi*sqrt(Lr*Cr));
Fr = 113E3;
Fs = 0.7*Fr;
Fs=70E3;
//Fs = Fr;

Lr =1/ (Fr^2 *4 *%pi^2 *Cr);
F = Fs/Fr;
//Lambda = Lr/Lm;
//Kx =  sqrt(Lambda/(1+Lambda));
Gamma = %pi/F;



Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;

// Lambda know
// 6 unknow variables
// mc(0)        x(1)
// mc(alpha)    x(2)
// iLr(0)       x(3)
// iLr(alpha)   x(4)
// M_max        x(5)   this changed to KNOWN Constant
// Lambda       x(5)   this changed to UNKNOWN
// alpha        x(6)

// equation (15a)
func1 = 'res(1) = -x(2) + (x(1)-1/M_max+1).*cos(x(6)) + x(3).*sin(x(6)) +1/M_max -1';

// equation (15c)
func2 = 'res(2) = -x(4) + ( -x(1) +1/M_max -1 ) .* sin(x(6)) + x(3) .*cos(x(6))';

// equation (15d)
func3 = 'res(3) = -x(4) + x(3) + x(5) .* x(6)';                                                                         

// equation (18a)
func4 = 'res(4) = x(1) + (x(2) -1/M_max) .* cos(sqrt(x(5)./(1+x(5))) .*(Gamma-x(6))) + (1/sqrt(x(5)./(1+x(5)))) .* x(4) .* sin(sqrt(x(5)./(1+x(5))).*(Gamma-x(6))) + 1/M_max';

// equation (18b)
func5 = 'res(5) = x(3) +  sqrt(x(5)./(1+x(5))) .*(-x(2) +1/M_max).* sin(sqrt(x(5)./(1+x(5))).*(Gamma -x(6))) + x(4) .* cos(sqrt(x(5)./(1+x(5))).*(Gamma -x(6)))' ;

// equation (18d)
func6 = 'res(6) =  (rL/Gamma) * ((-x(1) + 1/M_max -1) .* (1- cos(x(6)))  + x(3) .* sin(x(6)) -x(3) .* x(6) - 0.5 * x(5) .* x(6).^2) -1';




deff('res=DCMB_mode(x)',[func1; func2; func3; func4; func5; func6]);

// Initial condition
x1_0=0
x2_0=0
x3_0=0
x4_0=0
x5_0=0.2
x6_0=%pi

x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0];
xsol1 =fsolve(x0, DCMB_mode); 
res1 = DCMB_mode(xsol1) ;

mc_0= xsol1(1);
mc_alpha = xsol1(2);
I_lr_0 = xsol1(3);
I_lr_alpha = xsol1(4);
Lambda = xsol1(5);
Alpha = xsol1(6);

Vbase = 0.5*Vin*M_max;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
Vc_alpha = xsol1(2) * Vbase;
I_lr_0= xsol1(3) * Ibase;
I_lr_alpha= xsol1(4) * Ibase;

mm2_0 = (-mc_0 + 1/M_max )/(1+Lambda);
mm2_alpha = (-mc_alpha +1/M_max)/(1+Lambda);
mm2_gamma = (mc_0 -1/M_max)/(1+Lambda);

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


// Further process from [Wei2020]

