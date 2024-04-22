
funcprot(0);
exec('DCMB_PO_solver.sci', -1)
exec('CCMB_PN_solver.sci', -1)
exec('DCMB2_PON_solver.sci', -1)
exec('DCMAB_OPO_solver_ini.sci', -1)
exec('DCMAB_OPO_solver.sci', -1)
exec('DCMAB2_ONO_solver_ini.sci', -1)

Vin=280;

Rload =0.1; // PN
Rload =0.2; // PON
Rload =1;   // PO
Rload =6;   //OPO

N = 16;                                                                                                                                                                                                                 
Cr =25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;

Fn = 0.5


M_0=2 // M <1 
//alpha_0=Gamma *0.3 // critial alpha
//beta1_0=Gamma *0.6  // critical beta  x15_0 > x14_0

correctMode = %F;

if correctMode ==%F  then
    [mc_0, iLr_0, iLm_0, alpha, M, correctMode] = CCMB_PN_solver(Cr, Lr, Lm, N, Rload,Fn);
    printf('CCMB_PN %s \n', correctMode);
end

if correctMode == %F  then
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMB2_PON_solver(Cr, Lr, Lm, N, Rload,Fn);
    printf('DCMB2_PON %s \n', correctMode);
end

if correctMode == %F then
    [mc_0, iLr_0, iLm_0, alpha, M, correctMode] = DCMB_PO_solver(Cr, Lr, Lm, N, Rload,Fn);
    printf('DCMB_PO %s\n', correctMode);
end

if correctMode == %F then
    //[mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB_OPO_solver_ini(Cr, Lr, Lm, N, Rload,Fn, alpha_0, beta1_0, M_0);
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB_OPO_solver(Cr, Lr, Lm, N, Rload,Fn);
    printf('DCMAB_OPO %s \n', correctMode); 
end

if correctMode == %F then
    //[mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB2_ONO_solver_ini(Cr, Lr, Lm, N, Rload,Fn, alpha_0, beta1_0, M_0);
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB2_ONO_solver(Cr, Lr, Lm, N, Rload,Fn);
    printf('DCMAB2_ONO %s \n', correctMode); 
end

alpha_0 = alpha;
beta1_0 = beta1;
M_0 =M;

Vbase = 0.5*Vin*M;
Zbase = sqrt(Lr/Cr);
Ibase = Vbase/Zbase;

Vo = Vbase/N;
printf('%f', Vo)
