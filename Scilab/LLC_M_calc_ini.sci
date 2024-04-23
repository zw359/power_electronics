funcprot(0);
exec('DCMB_PO_solver.sci', -1)
exec('DCMB_PO_solver_ini2.sci', -1)
exec('CCMB_PN_solver.sci', -1)
exec('DCMB2_PON_solver.sci', -1)
exec('DCMB2_PON_solver_ini.sci', -1)
exec('DCMAB2_ONO_solver_ini.sci', -1)
exec('DCMAB_OPO_solver_ini2.sci', -1)
exec('DCMAB2_ONO_solver.sci', -1)
exec('DCMAB_OPO_solver.sci', -1)


function [x1, correctMode]= LLC_M_calc_ini(Cr, Lr, Lm, N, Rload,Fn,x0 )

correctMode = %F;

if correctMode ==%F  then
    [x1,correctMode] = CCMB_PN_solver_ini2(Cr, Lr, Lm, N, Rload,Fn,x0);
    printf('CCMB_PN %s \n', correctMode);
end



if correctMode == %F  then
    [x1, correctMode] = DCMB2_PON_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0)
    printf('DCMB2_PON %s \n', correctMode);
end


if correctMode == %F then
    [x1, correctMode] = DCMB_PO_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMB_PO %s\n', correctMode);
end


if correctMode == %F then
    [x1, correctMode] = DCMAB_OPO_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMAB_OPO %s \n', correctMode); 
end


if correctMode == %F then
    [mc_0, iLr_0, iLm_0, alpha, beta1, M, correctMode] = DCMAB2_ONO_solver_ini(Cr, Lr, Lm, N, Rload,Fn,  M_0, alpha_0, beta_0);
    printf('DCMAB2_ONO %s \n', correctMode); 
end

endfunction
