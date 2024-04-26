// Algorithm from paper [Yu2012]
// DCMB PO mode solver  most common below Fr
// circuit parameters
function [x1, correctMode] = DCMB_PO_solver_ini2(Cr, Lr, Lm, N, Rload,Fn, x0)
Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 
Fr = 1/(2*%pi*sqrt(Lr*Cr));

Lambda = Lr/Lm;
Gamma = %pi/Fn;

Zbase = sqrt(Lr/Cr);
rL= N^2*Rload/Zbase;

// x4 unknow variables
// mc(0)        x(1)
// iLr(0)       x(2)
// M            x(3)
// alpha        x(4)

// equation (18a)
func1 = 'res(1) = x(1) + ( ((x(1)-1/x(3)+1).*cos(x(4)) + x(2).*sin(x(4)) +1/x(3) -1) -1/x(3)) .* cos(Kx *(Gamma-x(4))) + (1/Kx) * (( -x(1) +1/x(3) -1 ) .* sin(x(4)) + x(2) .*cos(x(4)) ) .* sin(Kx*(Gamma-x(4))) + 1/x(3)';
// equation (18b)
func2 = 'res(2) = x(2) +  Kx *(-((x(1)-1/x(3)+1).*cos(x(4)) + x(2).*sin(x(4)) +1/x(3) -1) +1/x(3)).* sin(Kx*(Gamma -x(4))) + ((-x(1) +1/x(3) -1 ) .* sin(x(4)) + x(2) .*cos(x(4))) .* cos(Kx*(Gamma -x(4)))' ;
// equation (18c)
func3 = 'res(3) = ((-x(1) +1/x(3) -1 ) .* sin(x(4)) + x(2) .*cos(x(4))) - x(2) - Lambda*x(4) ';
// equation (18d)
//func4 = 'res(4) =  rL * ((-x(1) + 1/x(3) -1) .* (1- cos(x(4)))  + x(2) .* sin(x(4)) -x(2) .* x(4) - 0.5 * Lambda * x(4).^2) -Gamma';
func4 = 'res(4) =  (rL/Gamma) * ((-x(1) + 1/x(3) -1) .* (1- cos(x(4)))  + x(2) .* sin(x(4)) -x(2) .* x(4) - 0.5 * Lambda * x(4).^2) -1';
deff('res=DCMB_mode(x)',[func1; func2; func3; func4]);

x0_4 = [x0(1); x0(5); x0(13); x0(14)];
[xsol1, res1, info] =fsolve(x0_4, DCMB_mode); 


mc_0= xsol1(1);
I_lr_0 = xsol1(2);
M = xsol1(3);
Alpha = xsol1(4);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
I_lr_0= xsol1(2) * Ibase;

iL_alpha = (-xsol1(1)+1/xsol1(3)-1) * sin(xsol1(4)) + xsol1(2)*cos(xsol1(4));
mc_alpha = (xsol1(1)-1/xsol1(3)+1) * cos(xsol1(4)) + xsol1(2)*sin(xsol1(4)) +1/xsol1(3) -1;

mm2_0 = abs((-mc_0 + 1/M )/(1+Lambda));
mm2_alpha = abs((-mc_alpha +1/M)/(1+Lambda));
mm2_gamma = abs ((mc_0 +1/M)/(1+Lambda));


dt = 0.5;
I_lr_0_d = (-xsol1(1) + 1/xsol1(3) -1).*sin(xsol1(4)*dt) + xsol1(2).*cos(xsol1(4)*dt);
I_lm_0_d = xsol1(2) + Lambda * xsol1(4)*dt;

DCMB_PO = %F
if  (I_lr_0_d > I_lm_0_d) & (mm2_alpha<1) & (mm2_0 >1) & (mm2_gamma <=1) & (xsol1(3)>0) & (xsol1(4)>0) &(xsol1(4)< Gamma) & (info ==1) then
    DCMB_PO =%T
end
iLr_0 = xsol1(2); 
iLm_0 = xsol1(2);
alpha = xsol1(4);

correctMode = DCMB_PO;

//  mc(0)     mc(alpha)  mc(beta)   mc(gamma)  iLr(0)    iLr(alpha)   iLr(betta)   iLr(gamma)   iLm(0)     iLm(alpah)   iLm(beta)   iLm(gamma)      M       alpha     beta
x1=[xsol1(1);  mc_alpha; -xsol1(1); -xsol1(1);  xsol1(2);  iL_alpha;   -xsol1(2);    -xsol1(2);  xsol1(2);  iL_alpha;    -xsol1(2);    -xsol1(2) ; xsol1(3); xsol1(4); Gamma];


endfunction
