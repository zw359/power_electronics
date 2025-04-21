funcprot(0);
exec('DCMB_PO_solver.sci', -1)
exec('CCMB_PN_solver.sci', -1)
exec('DCMB2_PON_solver.sci', -1)
exec('DCMAB2_ONO_solver_ini.sci', -1)
exec('DCMAB_OPO_solver_ini.sci', -1)
exec('DCMAB2_ONO_solver.sci', -1)
exec('DCMAB_OPO_solver.sci', -1)
exec('DCMA2_NOP_solver_ini.sci',-1)

// circuit example taken from [Musa2013]
Vin = 390;
N =4;
Rload = 20;
Cr = 16.4E-9;
Lr = 35E-6;
Lm = 105E-6;
Fr = 1/(2*%pi*sqrt(Lr*Cr));
Fs = 300E3;
Fn = Fs/Fr;

// use default 15 variables
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

x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0; x7_0; x8_0; x9_0; x10_0; x11_0; x12_0; x13_0; x14_0; x15_0]
[xsol1, correctMode]= DCMA2_NOP_solver_ini(Cr, Lr, Lm, N, Rload,Fn, x0);

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
disp(Vo,correctMode, I_lr_gamma)
