// Algorithm from paper [Yu2012]
// DCMB PO mode solver  most common below Fr
// circuit parameters
function [mc_0, iLr_0, iLm_0, alpha, M, correctMode] = DCMB_PO_solver(Cr, Lr, Lm, N, Rload,Fn)
//Vin=280;
//Rload =5.4;
//N = 16;
//Cr =25E-9;
//Lr = 47.0212E-6
//Lm = 175.7023E-6;

Omega0= 1/sqrt(Lr*Cr);
Omega1= 1/sqrt((Lr+Lm)*Cr);
Kx= Omega1/Omega0; 
Fr = 1/(2*%pi*sqrt(Lr*Cr));

//Fn = Fs/Fr;
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

// Initial condition
x1_0=0
x2_0=0
x3_0=1.2  // M >1
x4_0=%pi/2  // alpha 

x0 = [x1_0; x2_0; x3_0; x4_0];
[xsol1, res1, info] =fsolve(x0, DCMB_mode); 


mc_0= xsol1(1);
I_lr_0 = xsol1(2);
M = xsol1(3);
Alpha = xsol1(4);

Vbase = 0.5*Vin*M;
Ibase = Vbase/Zbase;

Vo = Vbase/N;
Vc_0 = xsol1(1) * Vbase;
I_lr_0= xsol1(2) * Ibase;

mc_alpha = (xsol1(1)-1/xsol1(3)+1) * cos(xsol1(4)) + xsol1(2)*sin(xsol1(4)) +1/xsol1(3) -1;

mm2_0 = abs((-mc_0 + 1/M )/(1+Lambda));
mm2_alpha = abs((-mc_alpha +1/M)/(1+Lambda));
mm2_gamma = abs ((mc_0 +1/M)/(1+Lambda));


dt = 0.5;
I_lr_0_d = (-xsol1(1) + 1/xsol1(3) -1).*sin(xsol1(4)*dt) + xsol1(2).*cos(xsol1(4)*dt);
I_lm_0_d = xsol1(2) + Lambda * xsol1(4)*dt;

DCMB_PO = %F
if  (I_lr_0_d > I_lm_0_d) & (mm2_alpha<1) & (mm2_0 >1) & (mm2_gamma <=1)  & (info ==1) then
    DCMB_PO =%T
end
iLr_0 = xsol1(2); 
iLm_0 = xsol1(2);
alpha = xsol1(4);
correctMode = DCMB_PO;
/*
printf('DCM_PO  mode check %s \n', DCMB_PO)
printf('|mm2_0| \t |mm2_alpha| \t |mm2_gamma|\n')
printf('%f \t %f \t %f\n', abs(mm2_0), abs(mm2_alpha), abs(mm2_gamma));

printf('x1_0=%f\n', xsol1(1));
printf('x2_0=%f\n', xsol1(2));
printf('x3_0=%f\n', xsol1(3));
printf('x4_0=%f\n', xsol1(4));


printf('res1=%f\n', res1);
printf('Vo = %f\n', Vo);
*/

endfunction
