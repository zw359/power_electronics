// This function gives an approximate switching frequency Fn
// of a half-bridge LLC DC/DC converters
// Normalized switching frequency Fn = Fs/ Fr;
// Given inductance ratio h = Lm / Lr
// Given reciprocal of inductnace ratio h,  a = 1/h = Lr/Lm
// Given Q = Sqrt(Lr/Cr)/Rpri
// Given Rpri
// Given M = n*Vo/(0.5*Vin) for half bridge

function Fn = LLC_Gain_FHA_Fn(h,Q,M)
    // construct the Non-linear Equation
    // 
    func1 = "res = M * sqrt( (1+ 1/h- 1/(Fn^2 * h))^2 + Q^2 * (1/Fn -Fn)^2 ) -1";
    deff('res=LLC_Gain_FHA(Fn)', [func1]);
    xsol1 =fsolve(1, LLC_Gain_FHA); 
    res1 = LLC_Gain_FHA(xsol1);
    printf('Fn=%f\n', xsol1);  
    printf('res=%f\n', res1);
    Fn = xsol1;
endfunction
