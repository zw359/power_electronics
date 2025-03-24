/* Algorithm taken from paper [Wei2020] DCMB (PO) mode solver
There are many erors in the equations [13, 16, 20 28 30, 32 33 34]
They are all corrected in this script
The correctiveness verified by calculation and PSIM simulation
This "single non-linear equation" method claimed to be easier 
than solve systems of non-linear equations
*/

Vin=280;
Vo=48;
Po = 192;
Rload = Vo^2/Po;
N = Vin/(2*Vo);
N = 3.75;

//three parameters given
Fs=70E3;
Cr =60E-9;
K = 7.7365642;

// relation between Cr and theta 
// use "x" as theta for easy input
// These equation reaches less accurate solutions, compared with [Yu2012]
func1 = 'res= Cr * ( 8*N^2*Rload*Vin*Vo*Fs - 4*N*Rload*Vin^2*Fs + N*Rload*Vin^2*Fs*x^2 + 4*N*Rload*Vin^2*Fs*cos(x) +2*N^2*Rload*Vin*Vo*Fs*x^2 -8*N^2*Rload*Vin*Vo*cos(x) +N*Rload*Vin^2*Fs* x^2*cos(x)+2*N^2*Rload*Vin*Vo*Fs*x^2*cos(x) -8*N^2*Rload*Vin*Vo*Fs*x*sin(x))- (4*N*Vo^2 -2*Vin*Vo +N*Vo^2*x^2 -4*N*Vo^2*cos(x) +2*Vin*Vo*cos(x) +N*Vo^2*x^2*cos(x) -4*N*Vo^2*x*sin(x))';
  

k1=(Vin - Qi/Cr)/2-Vin+N*Vo;
k2 = ((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1);
k3 = ((Vin - Qi/Cr)/2-Vin+N*Vo) * cos(theta) + (((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1))*sin(theta) - N*Vo;
k4 = sqrt(K+1) * (-((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta) +(((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1))*cos(theta))

func1 = ' res = (K+1)*(((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1))^2 + (-((Vin - Qi/Cr)/2-Vin+N*Vo)-Vin+N*Vo)^2 - (((Vin - Qi/Cr)/2-Vin+N*Vo) * cos(theta) + (((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1))*sin(theta) - N*Vo)^2 - (sqrt(K+1) * (-((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta) +(((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1))*cos(theta)))^2 ';

deff('res=PO_mode_cal(theta)', [func1]);

// Initial condition

x0 = %pi;// use a quater of switching period
x0= 3.4866026;
xsol1 =fsolve(x0, PO_mode_cal); 
res1 = PO_mode_cal(xsol1) ;

printf('x1_0=%f\n', xsol1);  
printf('res1=%f\n', res1);
theta = xsol1;

Eo =Vo^2/(Rload*Fs);
Qi =Eo/Vin;
Vcr_0= (Vin - Qi/Cr)/2;
Vcr_0_1 = (Cr*Rload*Vin^2*Fs-Vo^2)/(2*Rload*Fs*Cr*Vin);

k1= Vcr_0 - Vin + N*Vo;
k2 = ((Vo*N*theta/K) + k1*sin(theta))/(cos(theta)-1);
k3= k1 * cos(theta) + k2*sin(theta) - N*Vo;
k4= sqrt(K+1) * (-k1*sin(theta) +k2*cos(theta));


k11=(Vin - Qi/Cr)/2-Vin+N*Vo;
k22 = ((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1);
k33 = ((Vin - Qi/Cr)/2-Vin+N*Vo) * cos(theta) + (((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1))*sin(theta) - N*Vo;
k44 = sqrt(K+1) * (-((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta) +(((Vo*N*theta/K) + ((Vin - Qi/Cr)/2-Vin+N*Vo)*sin(theta))/(cos(theta)-1))*cos(theta))



phi =  -acos((-k1-Vin+N*Vo)/sqrt(k3^2 + k4^2)) + atan(k4/k3) +%pi;

Fr = (Fs/%pi)*(theta + phi*sqrt(K+1));
Lr = 1/Cr/(Fr*2*%pi)^2;

