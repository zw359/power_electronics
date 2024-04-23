// Algorithm from paper [Yu2012]
// DCMB2 mode solver  same as PON mode

function [x1,correctMode] = CCMB_PN_solver_ini2(Cr, Lr, Lm, N, Rload,Fn,x0)   
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
// mc(gamma)    x(3)
// iLr(0)       x(4)
// iLr(alpha)   x(5)
// iLr(gamma)   x(6)
// iLm(0)       x(7)
// iLm(alpha)   x(8)
// iLm(gamma)   x(9)
// M            x(10)
// alpha        x(11)

//  mc(alpha)   P mode  0-- alpha
func1 = 'res(1) = -x(2) + (x(1)-1/x(10)+1).*cos(x(11)) + x(4).*sin(x(11)) +1/x(10) -1';
//  iLr(alpha)  P mode  0-- alpha
func2 = 'res(2) = -x(5)  +  (-x(1) +1/x(10) -1 ) .* sin(x(11)) + x(4) .*cos(x(11))';
//  iLm(alpha)  P mode  0-- alpha
func3 = 'res(3) = -x(8) + x(7) + Lambda * x(11)';
//  mc(gamma)    N mode alpha -- gamma
func4 = 'res(4) =  -x(3) + (x(2) -1/x(10) -1) .* cos(Gamma-x(11))  + x(5).* sin(Gamma-x(11)) + 1/x(10) +1' ;
//  iLr(gamma)  N mode alpha -- gamma
func5 = 'res(5) = -x(6) +  (-x(2) +1/x(10) +1) .* sin(Gamma-x(11))  + x(5) .* cos(Gamma-x(11))'; 
//  iLm(gamma)   N mode alpha -- gamma
func6 = 'res(6) = -x(9) + x(8) - Lambda * (Gamma-x(11))';
//  mc(0) + mc(gamma) =0
func7 = 'res(7) = x(1)+x(3)';
// iLr(0) + iLr(gamma) =0
func8 = 'res(8) = x(4) +x(6)';
// iLm(0) + iLm(gamma) =0
func9 = 'res(9) = x(7) + x(9)';
// iLr(alpha) - iLm(alpha) =0
func10 ='res(10) =  x(5) -x(8)' ;
// iout  voltage balance balance
func11 ='res(11) =  (( -x(1)+1/x(10)-1 ) .* (1-cos(x(11)))  + x(4) .*sin(x(11)) -x(7).*x(11) -0.5*Lambda*x(11).^2  + x(8) .*(Gamma-x(11)) - 0.5*Lambda*( Gamma -x(11)).^2 + ( -x(2) + 1/x(10) +1 ) .* ( cos(Gamma - x(11))-1) - x(5) .* sin(Gamma -x(11))) * rL - Gamma'; 
deff('res=CCMB_mode(x)',[func1; func2; func3; func4; func5; func6; func7; func8; func9; func10; func11]);

x0_11 = [x0(1); x0(2); x0(4); x0(5); x0(6);x0(8);x0(9); x0(10); x0(12); x0(13); x0(14)];
[xsol1, res1, info] =fsolve(x0_11, CCMB_mode, tol=1D-20); 


Vbase = 0.5*Vin*xsol1(10);
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
Vc_alpha = xsol1(2) * Vbase;
I_lr_0= xsol1(4) * Ibase;
I_lr_alpha= xsol1(5) * Ibase;

mm2_0 = abs((-xsol1(1) + 1/xsol1(10) )/(1+Lambda));
mm2_alpha = abs( (-xsol1(2) +1/xsol1(10))/(1+Lambda) );
mm2_gamma = abs( (-xsol1(3) + 1/xsol1(10))/(1+Lambda) ); // after the Gamma, Vin = 1/M


dt = 0.5
//  iLr(alpha)  P mode  0-- alpha
I_lr_0_d = (-xsol1(1) +1/xsol1(10) -1 ) .* sin(xsol1(11)*dt) + xsol1(4) .*cos(xsol1(11)*dt);
//  iLm(alpha)  P mode  0-- alpha
I_lm_0_d = xsol1(7) + Lambda * xsol1(11)*dt;
//  iLr(gamma)  N mode alpha -- gamma
I_lr_alpha_d = (-xsol1(2) +1/xsol1(10) +1) .* sin((Gamma-xsol1(11))*dt)  + xsol1(5) * cos((Gamma-xsol1(11))*dt); 
//  iLm(gamma)   N mode alpha -- gamma
I_lm_alpha_d = xsol1(8) - Lambda * ((Gamma-xsol1(11))*dt)';

CCMB_PN =%F;
if (I_lr_alpha_d < I_lm_alpha_d ) & ( I_lr_0_d > I_lm_0_d) &(mm2_alpha>1)  & (info ==1) then
   CCMB_PN = %T;
end
correctMode = CCMB_PN;
//  mc(0)     mc(alpha)  mc(beta)   mc(gamma)  iLr(0)    iLr(alpha)   iLr(betta)   iLr(gamma)   iLm(0)     iLm(alpah)   iLm(beta)   iLm(gamma)  M           alpha     beta
x1=[xsol1(1);  xsol1(2);  xsol1(3); xsol1(3);  xsol1(4);  xsol1(5);    xsol1(6);    xsol1(6);  xsol1(7);  xsol1(8);    xsol1(9);    xsol1(9) ; xsol1(10); xsol1(11); Gamma];
endfunction
