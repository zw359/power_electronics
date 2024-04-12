// Algorithm from paper [Yu2012]
// DCMB mode solver
// NOT work yet

Vin=280;
Rload = 1;
N = 16;
Vo=13.3;


Cr =25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;


Cr =35E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;
Lm = Lr *4.8;

Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 

Fr = 1/(2*%pi*sqrt(Lr*Cr));

Fs = Fr;
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


func1 = 'res= Cr * ( 8*N^2*Rload*Vin*Vo*Fs - 4*N*Rload*Vin^2*Fs + N*Rload*Vin^2*Fs*x^2 + 4*N*Rload*Vin^2*Fs*cos(x) +2*N^2*Rload*Vin*Vo*Fs*x^2 -8*N^2*Rload*Vin*Vo*cos(x) +N*Rload*Vin^2*Fs* x^2*cos(x)+2*N^2*Rload*Vin*Vo*Fs*x^2*cos(x) -8*N^2*Rload*Vin*Vo*Fs*x*sin(x))- (4*N*Vo^2 -2*Vin*Vo +N*Vo^2*x^2 -4*N*Vo^2*cos(x) +2*Vin*Vo*cos(x) +N*Vo^2*x^2*cos(x) -4*N*Vo^2*x*sin(x))';
  
  

deff('res=PO_mode_cal(x)', [func1]);

// Initial condition

x0 = 2;
xsol1 =fsolve(x0, PO_mode_cal); 
res1 = PO_mode_cal(xsol1) ;

printf('x1_0=%f\n', xsol1);
printf('res1=%f\n', res1);

