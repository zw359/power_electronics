// Algorithm from paper [Hu2016]
function [Lr, Lp, Fr, Fp, Lambda, Theta] = PON_calc(Vo, Rload, Vin_min, Fsw_min, N, Cr, Theta0, Lambda0)

// x1 = Lambda
// x2 = theta

K1 = 2*N*Cr*Rload*Vo*Vin_min*Fsw_min;
K2 = Cr*Rload*Vin_min^2*Fsw_min;
K3 = N*Vo*Rload*Fsw_min;
// equation (56)

func1 =  'res(1) = (-K1-K2+Vo^2)*cos(x(1)).*x(1) +(-K1-K2+Vo^2) * cos(x(1)).*x(2) + (K1-K2-Vo^2)*sin(x(2)) + (K1+K2-Vo^2)*sin(x(1))';

// equation (60)
func2= ' res(2) = 4* ((K1-K2 -Vo^2) * sin(x(2)) .*x(1) + (K1-K2 -Vo^2) * sin(x(2)) .* x(2) + (-K1 -K2 +Vo^2) * sin(x(1)) .* x(1) +(-K1 - K2 +Vo^2) *sin(x(1)) .* x(2) +(2*K1-2*K2 -2*Vo^2)*cos(x(2)) + (-2*K1-2*K2 +2*Vo^2)*cos(x(1)) + 4*K2) * K3* (x(1) +x(2))  - 4* ( (K1+K2-Vo^2)*cos(x(1)).*x(1) +  (K1+K2-Vo^2) * cos(x(1)) .* x(2) + (-K1+K2+Vo^2) *sin(x(2)) + (-K1-K2 +Vo^2) *sin(x(1))  + 2* Vo^2* x(1) +2*Vo^2*x(2)  ) *Rload *Fsw_min* Vin_min ';

deff('res=pon_calc(x)',[func1; func2]);

x0 = [Lambda0;Theta0];
xsol1 =fsolve(x0,pon_calc) 
res1 = pon_calc(xsol1) 

//printf('check the convergence of pon_calc \n')
//printf('res1=%f\n', res1);

Lambda = xsol1(1);
Theta = xsol1(2);
//res1 = pon_calc(x0); 

den1= K1 * sin(Theta) + K1 * sin(Lambda) -K2 * sin(Theta) +K2* sin(Lambda) - sin(Theta) *Vo^2 -sin(Lambda) *Vo^2;
num1 = -K1 * (Theta + Lambda);
K = num1/den1;


// equation (61) 
K4= ( K1 * cos(Theta) -K1 - K2*cos(Theta) -cos(Theta)*Vo^2 ) / (2*Rload*Fsw_min*Cr*Vin_min);
K5 =-sqrt(Lr*(K+1)*Cr) * sin(Theta) *(K1-K2-Vo^2) / (2*Cr*Rload*Fsw_min*Vin_min *sqrt(Lr*Cr)) ;
K6 = +Vin_min- (K*N*Vo +K *Vin_min + N*Vo)/K;

func3 = 'res = K4*cos(psi) + K5* sin(psi) +K6';
deff('res=psi_calc(psi)',func3);
x0=0;
xsol2 = fsolve(x0, psi_calc);
res2 = psi_calc(xsol2);

//printf('check the convergence of psi_calc \n')
//printf('res2=%f\n', res2);

psi = xsol2;

// try the closed from
psi1= atan(K5/K4) -acos(-K6/sqrt(K4^2 +K5^2));

Omega_r = 2*psi*sqrt(K+1)*Fsw_min + 2*Fsw_min*Lambda + 2*Fsw_min*Theta;
Omega_p = Omega_r/sqrt(K+1);
Fr = Omega_r/1000/2/%pi;
Fp = Omega_p/1000/2/%pi;
Lr = 1/(Cr*Omega_r^2);
Lp = K*Lr;
endfunction


