// Algorithm from paper [Yu2012]
// DCMA2 mode solver  same as NOP mode
// this a very narrow working mode
// set Rload at 2.3
// mm2_alpha <1
// mm2_beta >=1
// I_lr_0 < I_lm_0

Vin = 390;
N =4;
Rload = 50;
Cr = 16.4E-9;
Lr = 35E-6;
Lm = 105E-6;
Fr = 1/(2*%pi*sqrt(Lr*Cr));
Fs = Fr;
Fs = 300E3;
Fn = Fs/Fr;
Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 
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
deff('res=DCMA2_mode(x)',[func1; func2; func3; func4; func5; func6; func7; func8; func9; func10; func11; func12; func13; func14; func15]);

// Initial conditions
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

//===== solve the Non-linear equation systems
[xsol1, res1,info]= fsolve(x0, DCMA2_mode); 
//res1 = DCMA2_mode(xsol1) ;

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
I_lr_gamma = I_lr_gamma*Ibase;

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

printf('DCMA2_NOP mode check %s \n', string(DCMA2_NOP));
printf('|mm2_0| \t |mm2_alpha| \t |mm2_beta| \t |mm2_gamma|\n')
printf('%f \t %f \t %f \t %f\n', abs(mm2_0), abs(mm2_alpha), abs(mm2_beta), abs(mm2_gamma));
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
printf('x15_0=%f\n', xsol1(15));
printf('res1=%f\n', res1);
printf('Vo = %f\n', Vo);

W0 = Omega0;
W1 = Omega1;

// primay current i_lr "N"section  0 -- alpha
function i_lr = Lr_current_N(t)
     i_lr= (-mc_0 + 1/M +1)*sin(W0*t) + I_lr_0*cos(W0*t) ;
endfunction

function i_lr_rms = Lr_current_N_rms(t)
     i_lr= Lr_current_N(t);
     i_lr_rms  = i_lr^2;
endfunction

// primay current i_lr "O"section, alpha -- Beta
function i_lr = Lr_current_O(t)
     i_lr= (-mc_alpha + 1/M )*Kx*sin(Kx*(W0*t-alpha))+ I_lr_alpha* cos(Kx*(W0*t-alpha));
endfunction

function i_lr_rms = Lr_current_O_rms(t)
     i_lr= Lr_current_O(t);
     i_lr_rms  = i_lr^2;
endfunction

// primay current i_Lr "P"section, Beta -- gamma
function i_lr = Lr_current_P(t)
     i_lr =   (-mc_beta +1/M -1)* sin(W0*t-Beta) + I_lr_beta* cos(W0*t-Beta);
endfunction

function i_lr_rms = Lr_current_P_rms(t)
     i_lr= Lr_current_P(t);
     i_lr_rms  = i_lr^2;
endfunction


function i_sec = Sec_current_N(t)
     i_m = I_lm_0 - Lambda *W0*t;
     i_lr= (-mc_0 + 1/M +1)*sin(W0*t) + I_lr_0*cos(W0*t) ;
     i_sec = abs(i_lr -i_m)*N;
endfunction

function i_sec = Sec_current_P(t)
     i_m = I_lm_beta + Lambda *(W0*t-Beta);
     i_lr =   (-mc_beta +1/M -1)* sin(W0*t-Beta) + I_lr_beta* cos(W0*t-Beta);
     i_sec = abs(i_lr - i_m)*N;
endfunction


function i_sec_rms = Sec_current_N_rms(t)
     i_sec = Sec_current_N(t)
     i_sec_rms = i_sec ^2
endfunction

function i_sec_rms = Sec_current_P_rms(t)
     i_sec = Sec_current_P(t)
     i_sec_rms = i_sec ^2
endfunction

Integ1= intg(0, alpha/W0, Lr_current_N_rms); 
Integ2 = intg(alpha/W0, Beta/W0 , Lr_current_O_rms);
Integ3 = intg(Beta/W0, Gamma/W0, Lr_current_P_rms);
T= 1/Fs/2;
Lr_current_rms = Ibase* sqrt((Integ1 + Integ2 + Integ3)/T);
disp(Lr_current_rms);

Integ5= intg(0, alpha/W0, Sec_current_N_rms); 
Integ6 = intg(Beta/W0, Gamma/W0, Sec_current_P_rms);
Sec_current_rms = Ibase* sqrt((Integ5 + Integ6)/T);
disp(Sec_current_rms);

Integ7= intg(0, alpha/W0, Sec_current_N); 
Integ8 = intg(Beta/W0, Gamma/W0, Sec_current_P);
Sec_current_avg = Ibase*(Integ7 +Integ8)/T;
disp(Sec_current_avg);

Io = Vo/Rload;
Sec_ac_current_rms = sqrt(Sec_current_rms^2 - Io^2);
disp(Sec_ac_current_rms);

Iin = Vo^2/Rload/Vin;
Pri_ac_current_rms = sqrt(Lr_current_rms^2 - Iin^2)/sqrt(2);
disp(Pri_ac_current_rms);

t = linspace(0, alpha/W0, 100);
clf;plot(t,Ibase*pri_current_N(t));

t = linspace(alpha/W0, Beta/W0, 100);
plot(t,Ibase*pri_current_O(t));

t = linspace(Beta/W0, Gamma/W0, 100);
plot(t,Ibase*pri_current_P(t));
