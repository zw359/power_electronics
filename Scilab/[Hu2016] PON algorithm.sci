// Algorithm from paper [Hu2016]
Vo=12;
Po=600;
Rload = Vo^2/Po;
Io = Vo/Rload;
Cr = 16E-9;
Vin_min = 280;
Fsw_min = 100E3;
N = 16;

// x1 = alpha
// x2 = theta

K1 = 2*N*Cr*Rload*Vo*Vin_min*Fsw_min;
K2 = Cr*Rload*Vin_min^2*Fsw_min;
K3 = N*Vo*Rload*Fsw_min;
// equation (56)

func1 =  'res(1) = (-K1-K2+Vo^2)*cos(x(1)).*x(1) +(-K1-K2+Vo^2) * cos(x(1)).*x(2) + (K1-K2-Vo^2)*sin(x(2)) + (K1+K2-Vo^2)*sin(x(1))';

// equation (60)
func2= ' res(2) = 4* ((K1-K2 -Vo^2) * sin(x(2)) .*x(1) + (K1-K2 -Vo^2) * sin(x(2)) .* x(2) + (-K1 -K2 +Vo^2) * sin(x(1)) .* x(1) +(-K1 - K2 +Vo^2) *sin(x(1)) .* x(2) +(2*K1-2*K2 -2*Vo^2)*cos(x(2)) + (-2*K1-2*K2 +2*Vo^2)*cos(x(1)) + 4*K2) * K3* (x(1) +x(2))  - 4* ( (K1+K2-Vo^2)*cos(x(1)).*x(1) +  (K1+K2-Vo^2) * cos(x(1)) .* x(2) + (-K1+K2+Vo^2) *sin(x(2)) + (-K1-K2 +Vo^2) *sin(x(1))  + 2* Vo^2* x(1) +2*Vo^2*x(2)  ) *Rload *Fsw_min* Vin_min ';

//func1 = 'res(1)=x(2)-(x(1).^2+1)';
//func2 = 'res(2)=x(1)-(2*x(2)-x(2).^2)/3';

deff('res=pon_calc(x)',[func1; func2]);

Theta=2.635710
Alpha=1.089578

x0 = [Alpha;Theta];
xsol1 =fsolve(x0,pon_calc) 
res1 = pon_calc(xsol1) 

Theta = xsol1(2);
Alpha = xsol1(1);
//res1 = pon_calc(x0); 

den1= K1 * sin(Theta) + K1 * sin(Alpha) -K2 * sin(Theta) +K2* sin(Alpha) - sin(Theta) *Vo^2 -sin(Alpha) *Vo^2;

num1 = -K1 * (Theta + Alpha);

K = num1/den1;

printf('res1=%f\n', res1);
printf('Theta=%f\n', Theta);
printf('Alpha=%f\n', Alpha);
printf('K=%f\n', K);

