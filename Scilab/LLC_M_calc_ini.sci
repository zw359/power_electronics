funcprot(0);
exec('DCMB_PO_solver.sci', -1)
exec('DCMB_PO_solver_ini2.sci', -1)
exec('CCMB_PN_solver.sci', -1)
exec('CCMB_PN_solver_ini2.sci', -1)
exec('DCMB2_PON_solver.sci', -1)
exec('DCMB2_PON_solver_ini2.sci', -1)
exec('DCMAB2_ONO_solver_ini.sci', -1)
exec('DCMAB_OPO_solver_ini2.sci', -1)
exec('DCMAB2_ONO_solver.sci', -1)
exec('DCMAB_OPO_solver.sci', -1)

exec('DCMA_OP_solver_ini.sci', -1)
exec('DCMA2_NOP_solver_ini.sci', -1)
exec('CCMA_NP_solver_ini.sci', -1)


function [x1, correctMode]= LLC_M_calc_ini_below(Cr, Lr, Lm, N, Rload,Fn,x0 )

correctMode = %F;


if correctMode == %F  then
    [x1, correctMode] = DCMB2_PON_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMB2_PON %s \n', correctMode);
end


if correctMode == %F then
    [x1, correctMode] = DCMB_PO_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMB_PO %s\n', correctMode);
end


if correctMode ==%F  then
    [x1,correctMode] = CCMB_PN_solver_ini2(Cr, Lr, Lm, N, Rload,Fn,x0);
    printf('CCMB_PN %s \n', correctMode);
end


if correctMode == %F then
    [x1, correctMode] = DCMAB_OPO_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMAB_OPO %s \n', correctMode); 
end


if correctMode == %F then
    [x1, correctMode] = DCMAB2_ONO_solver_ini(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMAB2_ONO %s \n', correctMode); 
end

endfunction



function [x1, correctMode]= LLC_M_calc_ini_above(Cr, Lr, Lm, N, Rload,Fn,x0 )

correctMode = %F;


if correctMode == %F then
    [x1, correctMode] = DCMAB_OPO_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMAB_OPO %s \n', correctMode); 
end


if correctMode == %F then
    [x1, correctMode] = DCMAB2_ONO_solver_ini(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMAB2_ONO %s \n', correctMode); 
end

if correctMode == %F then
    [x1,correctMode]=  DCMA_OP_solver_ini(Cr, Lr, Lm, N, Rload,Fn,x0)
    printf('DCMA_OP %s \n', correctMode); 
end

if correctMode == %F then
    [x1, correctMode]= DCMA2_NOP_solver_ini(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('DCMA2_NOP %s \n', correctMode); 
end

if correctMode == %F then
     [x1, correctMode] = CCMA_NP_solver_ini(Cr, Lr, Lm, N, Rload,Fn, x0);
    printf('CCMA_NP %s \n', correctMode); 
end

endfunction
