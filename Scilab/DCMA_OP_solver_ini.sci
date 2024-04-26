// Algorithm ref [Yu2012] time base analysis
// DCMBA_OP mode solver
// Fs > Fr,  lightload
// Set Rload =5, Fs=180k
// mmc_alpha =1
function [x1,correctMode]=  DCMA_OP_solver_ini(Cr, Lr, Lm, N, Rload,Fn,x0)
Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 
Fr = 1/(2*%pi*sqrt(Lr*Cr));
Lambda = Lr/Lm;
Gamma = %pi/Fn;
Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;

// general variables  TWO sections
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

//  mc(alpha)   O mode  0-- alpha
func1 = 'res(1) = -x(2) + (x(1)-1/x(10)) .* cos(Kx*x(11))  + (x(4)/Kx) .* sin(Kx*x(11)) + 1/x(10)';
//  iLr(alpha)  O mode  0-- alpha
func2 = 'res(2) = -x(5) + (-x(1) +1/x(10)) *Kx .* sin(Kx*x(11))  + x(4) .* cos(Kx*x(11))';
//  iLm(alpha)  O mode  0-- alpha
func3 = 'res(3) = -x(8) + x(5)';
//  mc(gamma)    P mode alpha -- gamma
func4 = 'res(4) =  -x(3) + (x(2)-1/x(10)+1) .* cos(Gamma-x(11)) + x(5).*sin(Gamma-x(11)) +1/x(10) -1';
// iLr(gamma)   P mode alpha -- gamma
func5 = 'res(5) =  -x(6) + (-x(2) +1/x(10) -1 ) .* sin(Gamma-x(11)) + x(5) .*cos(Gamma-x(11))' ;
// iLm(gamma)   P mode alpha -- gamma
func6 = 'res(6) =-x(9) + x(8) + Lambda * (Gamma-x(11))';
//  mc(0) + mc(gamma) =0
func7 = 'res(7) = x(1)+x(3)';
// iLr(0) + iLr(gamma) =0
func8 = 'res(8) = x(4) +x(6)';
// iLm(0) + iLm(gamma) =0
func9 = 'res(9) = x(7) + x(9)';
//  mm(alpha) = 1
//  mc(alpha) = 1/M - Lambda -1
func10 = 'res(10) = - 1/x(10) +x(2) + Lambda + 1';
//equation iout balance
func11 = 'res(11) = ( (-x(2)+1/x(10)-1 ).*(1-cos(Gamma-x(11))) + x(5) .*sin(Gamma-x(11)) - x(8).*(Gamma-x(11)) -0.5*Lambda*(Gamma-x(11)).^2 ) * rL -Gamma'; 
deff('res=DCMA_mode(x)',[func1; func2; func3; func4; func5; func6; func7; func8; func9; func10; func11]);


x0_4 = [x0(1); x0(2); x0(3); x0(5); x0(6); x0(7); x0(9); x0(10); x0(11); x0(13); x0(14)];
[xsol1, res1, info] =fsolve(x0_4, DCMA_mode); 

mc_0= xsol1(1);
mc_alpha = xsol1(2);
mc_gamma = xsol1(3);
I_lr_0 = xsol1(4);
I_lr_alpha = xsol1(5);
I_lr_gamma = xsol1(6);
I_lm_0 =xsol1(7);
I_lm_alpha =xsol1(8)
I_lm_gamma =xsol1(9)
M = xsol1(10);
Alpha = xsol1(11);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = mc_0 * Vbase;
Vc_alpha = mc_alpha * Vbase;


mm2_0 = abs((-mc_0 + 1/M )/(1+Lambda));
mm2_alpha = abs((-mc_alpha +1/M)/(1+Lambda));
mm2_gamma = abs((-mc_gamma + 1/M)/(1+Lambda));

// copy from func5 func6
dt=0.5;

// iLr(gamma)   P mode alpha -- gamma
I_lr_alpha_d  = (-xsol1(2) +1/xsol1(10) -1 ) .* sin((Gamma-xsol1(11))*dt) + xsol1(5) .*cos((Gamma-xsol1(11))*dt);
// iLm(gamma)   P mode alpha -- gamma
I_lm_alpha_d  = xsol1(8) + Lambda * (Gamma-xsol1(11))*dt;


DCMA_OP =%F;
if (I_lr_alpha_d > I_lm_alpha_d ) & (xsol1(6) == xsol1(9)) & (xsol1(11)>0) & (xsol1(11)<Gamma) & (xsol1(10)>0) & (info==1) then
    DCMA_OP = %T;
end
correctMode = DCMA_OP;
//  mc(0)     mc(alpha)  mc(beta)   mc(gamma)  iLr(0)    iLr(alpha)   iLr(betta)   iLr(gamma)   iLm(0)     iLm(alpah)   iLm(beta)   iLm(gamma)  M           alpha     beta
x1=[xsol1(1);  xsol1(2);  xsol1(3); xsol1(3);  xsol1(4);  xsol1(5);    xsol1(6);    xsol1(6);  xsol1(7);  xsol1(8);    xsol1(9);    xsol1(9) ; xsol1(10); xsol1(11); Gamma];
endfunction
