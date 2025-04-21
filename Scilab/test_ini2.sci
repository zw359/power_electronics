
funcprot(0);
exec('DCMB_PO_solver.sci', -1)
exec('CCMB_PN_solver.sci', -1)
exec('CCMB_PN_solver_ini2.sci', -1)

exec('DCMB2_PON_solver.sci', -1)
exec('DCMB2_PON_solver_ini.sci', -1)

exec('DCMAB_OPO_solver_ini.sci', -1)
exec('DCMAB_OPO_solver.sci', -1)

exec('DCMAB2_ONO_solver_ini.sci', -1)
exec('DCMAB2_ONO_solver.sci', -1)



Vin=280;

Rload =0.1; // PN
Rload =0.2; // PON
Rload =1;   // PO
Rload =1;   //OPO

N = 16;                                                                                                                                                                                                                 
Cr =25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;

Fn = 0.5


//M_0=1 // M <1 
//alpha_0=Gamma *0.3 // critial alpha
//beta1_0=Gamma *0.6  // critical beta  x15_0 > x14_0

// The initial condisions are reserved for 3-section modes
// there four time values: 0, alpha, beta, gamma
x1_0=0;  // mc(0)
x2_0=0;  // mc(alpha)
x3_0=0;  // mc(beta)
x4_0=0;  // mc(gamma)
x5_0=0;  // iLr(0)
x6_0=0;  // iLr(alpha)
x7_0=0;  // iLr(beta)
x8_0=0;  // iLr(gamma)
x9_0=0;  // iLm(0)
x10_0=0; // iLm(alpha)
x11_0=0; // iLm(beta)
x12_0=0; // iLm(gamma)
x13_0= 1.2 // M <1 
x14_0=%pi*0.3 // critial alpha
x15_0=%pi*0.6  // critical beta  x15_0 > x14_0

x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0; x7_0; x8_0; x9_0; x10_0; x11_0; x12_0; x13_0; x14_0; x15_0];

correctMode = %F;

if correctMode ==%F  then
    [x1,correctMode] = CCMB_PN_solver_ini2(Cr, Lr, Lm, N, Rload,Fn,x0);
    printf('CCMB_PN %s \n', correctMode);
end



if correctMode == %F  then
    [x1, correctMode] = DCMB2_PON_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMB2_PON %s \n', correctMode);
end


if correctMode == %F then
    [x1, correctMode] = DCMB_PO_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMB_PO %s\n', correctMode);
end


if correctMode == %F then
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB_OPO_solver_ini(Cr, Lr, Lm, N, Rload,Fn, M_0, alpha_0, beta_0 );
    printf('DCMAB_OPO %s \n', correctMode); 
end


if correctMode == %F then
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB2_ONO_solver_ini(Cr, Lr, Lm, N, Rload,Fn,  M_0, alpha_0, beta_0);
    printf('DCMAB2_ONO %s \n', correctMode); 
end


Vbase = 0.5*Vin*M;
Zbase = sqrt(Lr/Cr);
Ibase = Vbase/Zbase;

Vo = Vbase/N;
printf('%f', Vo)
