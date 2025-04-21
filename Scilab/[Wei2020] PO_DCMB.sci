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
K = 10;

Eo =Vo^2/(Rload*Fs);
Qi =Eo/Vin;
Vcr_0= (Vin - Qi/Cr)/2;
Vcr_0_1 = (Cr*Rload*Vin^2*Fs-Vo^2)/(2*Rload*Fs*Cr*Vin);


// relation between Cr and theta 
// use "x" as theta for easy input
// These equation reaches less accurate solutions, compared with [Yu2012]
func1 = 'res= Cr * ( 8*N^2*Rload*Vin*Vo*Fs - 4*N*Rload*Vin^2*Fs + N*Rload*Vin^2*Fs*x^2 + 4*N*Rload*Vin^2*Fs*cos(x) +2*N^2*Rload*Vin*Vo*Fs*x^2 -8*N^2*Rload*Vin*Vo*cos(x) +N*Rload*Vin^2*Fs* x^2*cos(x)+2*N^2*Rload*Vin*Vo*Fs*x^2*cos(x) -8*N^2*Rload*Vin*Vo*Fs*x*sin(x))- (4*N*Vo^2 -2*Vin*Vo +N*Vo^2*x^2 -4*N*Vo^2*cos(x) +2*Vin*Vo*cos(x) +N*Vo^2*x^2*cos(x) -4*N*Vo^2*x*sin(x))';
  

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
Wr = 2*%pi*Fr;
Wm = Wr/sqrt(K+1);

Lr = 1/Cr/(Fr*2*%pi)^2;
Lm = Lr * K;

//-------function  iL(t) ----
function y = pri_current(x)
    x1 = x(0 <= x & x <= theta/Wr);
    y(find(0 <= x & x <= theta/Wr)) = Cr * Wr*(-k1*sin(Wr*x1) + k2*cos(Wr*x1));
    
    x2 = x(theta/Wr<x & x<=theta/Wr+phi/Wm);
    y(find(theta/Wr<x & x<=theta/Wr+phi/Wm)) = Cr * Wm *(-k3*sin(x2*Wm-theta/sqrt(K+1)) + k4*cos(x2*Wm-theta/sqrt(K+1)));
    
    x3 = x(theta/Wr+phi/Wm<x & x<=2*theta/Wr+phi/Wm);
    y(find(theta/Wr+phi/Wm<x & x<=2*theta/Wr+phi/Wm)) = -(Cr * Wr*(-k1*sin(x3*Wr-theta-phi*sqrt(K+1)) + k2*cos(x3*Wr-theta-phi*sqrt(K+1))));
    
    x4 = x(2*theta/Wr +phi/Wm<x & x<=2*theta/Wr+2*phi/Wm);
    y(find(2*theta/Wr +phi/Wm<x & x<=2*theta/Wr+2*phi/Wm)) = -(Cr * Wm *(-k3*sin(x4*Wm-2*theta/sqrt(K+1)-phi) + k4*cos(x4*Wm-2*theta/sqrt(K+1)-phi)));
endfunction

x=linspace(0, 2*theta/Wr+2*phi/Wm, 400);
y = pri_current(x);
clf; plot(x,y)


// primay current "P"section. [Wei2020] equation 5
function y = pri_current1(t)
     y= (Cr * Wr*(-k1*sin(Wr*t) + k2*cos(Wr*t)))^2 ;
endfunction
// primary current "0" section. [Wei 2020] equation 31
function y = pri_current2(t)
     y= (Cr * Wm *(-k3*sin(Wm*(t-theta/Wr)) + k4*cos(Wm*(t-theta/Wr))))^2;
endfunction

// secondary current "P" rms value
function y = sec_current1(t)
    y = ((Cr * Wr*(-k1*sin(Wr*t) + k2*cos(Wr*t))) - (Cr*Wr*k2  + t* N*Vo/Lm))^2*N^2;
endfunction

// secondary current "P" avg value
function y = sec_current2(t)
    y = ((Cr * Wr*(-k1*sin(Wr*t) + k2*cos(Wr*t))) - (Cr*Wr*k2  + t* N*Vo/Lm))*N;
endfunction


// DC gain function based on FHA
function y = M_fha(h,Q, fn)
  //  y = 1/sqrt((1+1/h-1/fn^2/h)^2 + Q^2*(1/fn -fn)^2);
    y = 1./sqrt((1+1/h-1./(fn^2 *h))^2 + Q^2*(1./fn -fn)^2);
endfunction

h = K;
Rsc = 0.1;
Rsec = 8*Rsc/%pi^2;
Rpri = N^2*Rsc;
Q = sqrt(Lr/Cr)/Rpri;
fn = 270E3/Fr;

fn = linspace(0.5,2.5, 300);
//clf; plot(fn, 1./sqrt((1+1/h-1./(fn^2 *h))^2 + Q^2*(1./fn -fn)^2));
clf; plot(fn, M_fha(h,Q,fn));
//M_sc = M_fha(h,Q, fn);

// primary rms calculation
T = theta/Wr+phi/Wm;    // half switching period
time_step = 10E-9;
Integ1= intg(0, theta/Wr, pri_current1); 
Integ2 = intg(theta/Wr, theta/Wr+phi/Wm, pri_current2);
Integ3 = intg(0, theta/Wr, sec_current1);
pri_current_rms = sqrt((Integ1 + Integ2)/T);// caculation based on half Ts
sec_current_rms = sqrt(Integ3/T);
sec_current_avg = intg(0, theta/Wr, sec_current2)/T;




