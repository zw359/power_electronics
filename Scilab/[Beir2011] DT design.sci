Cp = 550E-12;
Cr =9.9E-9;
Lr = 165E-6;
a = 1.44;
Lm = Lr/ a;
Vin =  370;
Fr =  1/ (2* %pi * sqrt(Lr* Cr));
Omega_0 = 1/sqrt((Lm + Lr)*Cr);
DTn = (sqrt(2)/(4*%pi)) * sqrt(sqrt(1+ 16*(Cp/Cr)^2) -1);
Fn_max = %pi/ ( sqrt(1+1/a) * atan(1/(2*%pi*DTn)));
DT = DTn * 2 * %pi / Omega_0;
disp(DT);

//  a simple calculation assuming the current is like a triangule waveform

Fs = 205E3;
Fs = 221E3;
Ip_pk = Vin/(Lr+Lm)/(4*Fs);
DT2 = Cp *Vin / Ip_pk;
disp(DT2);
