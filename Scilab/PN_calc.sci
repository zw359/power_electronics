// Algorithm from paper [Hu2016]
function [Lr, Lp, Fr, Fp, Vmargin, Lambda, Theta, K] = PN_calc(Vo, Rload, Vin_min, Fsw_min, N, Cr)
    

num1 = Vin_min*(4*N^2*Cr*Rload*Vo*Fsw_min - 2*N*Cr*Rload*Vin_min*Fsw_min + Vo);
temp1 = (2*N*Cr*Rload*Vo*Vin_min*Fsw_min - Cr*Rload*Vin_min^2 * Fsw_min -Vo^2)
den1 = 2*N*temp1;
Theta = acos(num1/den1);

// Lambda equation(27)
num2 = (2*N*Cr*Rload*Vo*Vin_min*Fsw_min - Cr*Rload*Vin_min^2*Fsw_min - Vo^2)* sin(Theta);
den2 = 2*N*Cr*Rload*Vo*Vin_min*Fsw_min + Cr*Rload*Vin_min^2*Fsw_min - Vo^2;
Lambda = asin(num2/den2);

// K euqation(31)
num3 = -N*Cr*Rload*Vo*Vin_min*Fsw_min*(Lambda + Theta);
den3 = num2;
K = num3/den3;

// Omega_r   equation (33)
Omega_r = 2*(Lambda + Theta)*Fsw_min;
Fr = Omega_r/1000/2/%pi;
// Lr   equation(34)
Lr = 1/(Cr*Omega_r^2);

// Lp  equation (35)
Lp = K*Lr;
Fp =1/(2*%pi*sqrt(Cr*(Lr+Lp)))/1000; //kHz
// margin euqation(38)
num4= temp1 *cos(Theta);
den4= 2*Rload*Fsw_min*Cr*Vin_min;

Vmargin = num4/den4 - N*Vo - N*Vo*(K+1)/K;



endfunction


