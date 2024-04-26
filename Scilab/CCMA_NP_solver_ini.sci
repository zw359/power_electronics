// Algorithm from paper [Hu2016]
// CCMA  NP mode   most common mode when Fs > Fr
// Norminal load condition  set at Rload =1
// mm2_0 >1
// mm2_alpha >1
function [x1, correctMode] = CCMA_NP_solver_ini(Cr, Lr, Lm, N, Rload,Fn, x0)
Fr = 1/(2*%pi*sqrt(Lr*Cr));

Lambda = Lr/Lm;
Gamma = %pi/Fn;

Im_alpha = -Lambda*Gamma/2;  // i_lm(Alpha)
Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;

// 6 unknow variables
// mc(0)        x(1)
// mc(alpha)    x(2)
// iLr(0)       x(3)
// iLm(0)       x(4)
// M            x(5)
// alpha        x(6)

// equation (25a)
func1 =  'res(1) = x(1) + (x(2) -1/x(5) +1) .* cos(Gamma-x(6)) + Im_alpha .* sin(Gamma-x(6)) + 1/x(5)-1';
// equation (22a)
func2 =  'res(2) = -x(2) + (x(1) -1/x(5) -1).*cos(x(6)) + x(3) .* sin(x(6)) + 1/x(5) +1' ;
// equation (22d)
func3 =  'res(3) = x(4)- Lambda * x(6) - Im_alpha';
// equation (25b)
func4 =  'res(4) = x(3) + (-x(2) +1/x(5) -1) .* sin(Gamma-x(6)) + Im_alpha .* cos(Gamma-x(6))';
// equation (22c)
func5 =  'res(5) = -Im_alpha + (-x(1) +1/x(5) +1) .* sin(x(6)) + x(3) .* cos(x(6))';
// euqation (25c)
func6 =  'res(6) = rL * (1/Gamma) * ( x(6) .* x(4) - 0.5* Lambda * x(6).^2 - (-x(1) +1/x(5) +1) .* (1-cos(x(6)))  - x(3).*sin(x(6)) + (-x(2)+1/x(5)-1) .* (1-cos(Gamma-x(6))) +Im_alpha * sin(Gamma-x(6)) - Im_alpha*(Gamma-x(6)) -0.5*Lambda*(Gamma-x(6)).^2) -1' ; 
deff('res=CCMA_NP_mode(x)',[func1; func2; func3; func4; func5; func6]);

x0_4 =[x0(1); x0(2) ; x0(5);  x0(9); x0(13); x0(14)];
[xsol1, res1, info] = fsolve(x0_4, CCMA_NP_mode); 

mc_0= xsol1(1);
mc_alpha = xsol1(2);
I_lr_0 = xsol1(3);
I_lm_0 = xsol1(4);
M = xsol1(5);
Alpha = xsol1(6);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
Vc_alpha = xsol1(2) * Vbase;

mm2_0 = (-mc_0 + 1/M )/(1+Lambda);
mm2_alpha = (-mc_alpha +1/M)/(1+Lambda);
mm2_gamma = (mc_0 +1/M)/(1+Lambda);

// copy from func5 , func6
dt=0.5;
I_lr_0_d =  (-xsol1(1)+1/xsol1(5)+1) *sin(xsol1(6) * dt) + xsol1(3) * cos(xsol1(6) * dt);
I_lm_0_d =xsol1(4) - Lambda * xsol1(6) * dt;
I_lr_alpha_d  = (-xsol1(1)+1/xsol1(5)-1) * sin((Gamma-xsol1(6))*dt) + Im_alpha * cos((Gamma-xsol1(6))*dt) ;
I_lm_alpha_d  = Im_alpha + Lambda*(Gamma-xsol1(6))*dt;


CCMA_NP =%F;
if (I_lr_0_d < I_lm_0_d ) & (I_lr_alpha_d > I_lm_alpha_d ) & (xsol1(5)>0) & (xsol1(6)>0)  &  (xsol1(6)<Gamma) & (mm2_alpha >1) then
    CCMA_NP = %T;
end
correctMode = CCMA_NP;

//  mc(0)     mc(alpha)  mc(beta)   mc(gamma)  iLr(0)    iLr(alpha)   iLr(betta)   iLr(gamma)   iLm(0)     iLm(alpah)   iLm(beta)   iLm(gamma)      M       alpha     beta
x1=[xsol1(1);  xsol1(2); -xsol1(1); -xsol1(1);  xsol1(3); Im_alpha;   -xsol1(3);    -xsol1(3);  xsol1(4);  Im_alpha;    -xsol1(4);    -xsol1(4) ; xsol1(5); xsol1(6); Gamma];
endfunction
