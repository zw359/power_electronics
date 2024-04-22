funcprot(0);
exec('DCMB_PO_solver.sci', -1)
exec('CCMB_PN_solver.sci', -1)
exec('DCMB2_PON_solver.sci', -1)
exec('DCMAB2_ONO_solver_ini.sci', -1)
exec('DCMAB_OPO_solver_ini.sci', -1)
exec('DCMAB2_ONO_solver.sci', -1)
exec('DCMAB_OPO_solver.sci', -1)


function [M_final, alpha_final, beta1_final]= LLC_M_calc(Cr, Lr, Lm, N, Rload,Fn, alpha_0, beta1_0, M_0 )

/*
Vin=280;

Rload =0.1; // PN
Rload =0.2; // PON
Rload =1;   // PO
Rload =10;   //OPO

N = 16;
Cr =25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;

Fn =  0.6131112;
*/
correctMode = %F;
alpha=0; 
beta1=0; 
M =1;


if correctMode ==%F  then
    [mc_0, iLr_0, iLm_0, alpha, M, correctMode] = CCMB_PN_solver(Cr, Lr, Lm, N, Rload,Fn);
    //printf('CCMB_PN %s \n', correctMode);
end



if correctMode == %F  then
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMB2_PON_solver(Cr, Lr, Lm, N, Rload,Fn);
    //printf('DCMB2_PON %s \n', correctMode);
end


if correctMode == %F then
    [mc_0, iLr_0, iLm_0, alpha, M, correctMode] = DCMB_PO_solver(Cr, Lr, Lm, N, Rload,Fn);
    //printf('DCMB_PO %s\n', correctMode);
end


if correctMode == %F then
    //[mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB_OPO_solver_ini(Cr, Lr, Lm, N, Rload,Fn, alpha_0, beta1_0, M_0 );
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB_OPO_solver(Cr, Lr, Lm, N, Rload,Fn);
    //printf('DCMAB_OPO %s \n', correctMode); 
end


if correctMode == %F then
    //[mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB2_ONO_solver_ini(Cr, Lr, Lm, N, Rload,Fn, alpha_0, beta1_0, M_0);
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB2_ONO_solver(Cr, Lr, Lm, N, Rload,Fn);
    //printf('DCMAB2_ONO %s \n', correctMode); 
end

M_final =M;
alpha_final = alpha;
beta1_final = beta1;
endfunction
