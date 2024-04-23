// Algorithm from paper [Yu2012]
// DCMB2 mode solver  same as PON mode
function [x1, correctMode] = DCMB2_PON_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0)

Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 
Fr = 1/(2*%pi*sqrt(Lr*Cr));
Lambda = Lr/Lm;
Gamma = %pi/Fn;
Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;

// general variables 
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
// beta         x(15)

//  mc(alpha)   P mode  0-- alpha
func1 = 'res(1) = -x(2) + (x(1)-1/x(13)+1).*cos(x(14)) + x(5).*sin(x(14)) +1/x(13) -1';
//  iLr(alpha)  P mode  0-- alpha
func2 = 'res(2) = -x(6) + (-x(1) +1/x(13) -1 ) .* sin(x(14)) + x(5) .*cos(x(14))';
//  iLm(alpha)  P mode  0-- alpha
func3 = 'res(3) = -x(10) + x(9) + Lambda * x(14)';
//  mc(beta)    O mode alpha -- beta
func4 = 'res(4) =  -x(3) + (x(2)-1/x(13)) .* cos(Kx*(x(15)-x(14)))  + (x(6)/Kx) .* sin(Kx*(x(15)-x(14))) + 1/x(13)';
// iLr(beta)   O mode alpha -- beta
func5 = 'res(5) =  -x(7) + (-x(2) +1/x(13)) *Kx .* sin(Kx*(x(15)-x(14)))  + x(6) .* cos(Kx*(x(15)-x(14)))' ;
//  mc(gamma)    N mode beta -- gamma
func6 = 'res(6) =  -x(4) + (x(3) -1/x(13) -1) .* cos(Gamma-x(15))  + x(7) * sin(Gamma-x(15)) + 1/x(13) +1' ;
//  iLr(gamma)  N mode beta -- gamma
func7 = 'res(7) = -x(8) + (-x(3) +1/x(13) +1) .* sin(Gamma-x(15))  + x(7) * cos(Gamma-x(15))'; 
//  iLm(gamma)   N mode.  beta -- gamma
func8 = 'res(8) = -x(12)  + x(11) - Lambda * (Gamma-x(15))';
//  mc(0) + mc(gamma) =0
func9 = 'res(9) = x(1)+x(4)';
// iLr(0) + iLr(gamma) =0
func10 = 'res(10) = x(5) +x(8)';
// iLm(0) + iLm(gamma) =0
func11 = 'res(11) = x(9) + x(12)';
// iLr(alpha) - iLm(alpha) =0
func12 =' res(12) =  x(6) -x(10)' ;
// iLr(beta) - iLm(beta) =0
func13 = 'res(13) = x(7) - x(11)';
//  mm(beta) =-1
//  mc(beta) = 1/M + Lambda+1
func14 = 'res(14) = -x(3)+ 1/x(13) + Lambda + 1';
//equation iout balance
func15 = 'res(15) = (( -x(1)+1/x(13)-1 ) .* (1-cos(x(14)))  + x(5) .*sin(x(14)) -x(9).*x(14) -0.5*Lambda*x(14).^2 + x(11) .*(Gamma-x(15)) - 0.5*Lambda*( Gamma -x(15)) .^2 + ( -x(3) + 1/x(13) +1 ) .* (cos(Gamma - x(15))-1) - x(7) .* sin(Gamma -x(15))  ) * rL -Gamma'; 
deff('res=DCMB2_PON_mode(x)',[func1; func2; func3; func4; func5; func6; func7; func8; func9; func10; func11; func12; func13; func14; func15]);

[xsol1, res1, info] =fsolve(x0, DCMB2_PON_mode); 


Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;
Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
Vc_alpha = xsol1(2) * Vbase;


mm2_0 = abs((-xsol1(1) + 1/M )/(1+Lambda));
mm2_alpha = abs((-xsol1(2) +1/M)/(1+Lambda));
mm2_gamma = abs((-xsol1(3) +1/M)/(1+Lambda));

dt = 0.5;
I_lr_0_d = (-xsol1(1) + 1/xsol1(13) -1).*sin(xsol1(14)*dt) + xsol1(5).*cos(xsol1(14)*dt);
I_lm_0_d = xsol1(9) + Lambda * xsol1(14)*dt;

//  iLr(gamma)  N mode beta -- gamma
I_lr_beta_d = (-xsol1(3) +1/xsol1(13) +1) .* sin((Gamma-xsol1(15))*dt)  + xsol1(7) * cos((Gamma-xsol1(15))*dt); 
//  iLm(gamma)   N mode.  beta -- gamma
I_lm_beta_d =  xsol1(11) - Lambda * (Gamma-xsol1(15))*dt;

DCMB2_PON = %F
if  (I_lr_0_d > I_lm_0_d) & (I_lr_beta_d < I_lm_beta_d) & (mm2_alpha<1) & (xsol1(15)< Gamma) &( xsol1(14)< xsol1(15)) & (xsol1(15) >0) & (xsol1(14)>0) & (info==1) then
    DCMB2_PON =%T
end

correctMode = DCMB2_PON;
x1=xsol1;

endfunction
