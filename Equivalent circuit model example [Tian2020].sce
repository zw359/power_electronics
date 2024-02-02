Vin=50; // half bridge
Rload= 1;
Lr=360E-9;
Cr=110E-9;
Lm=2.1E-6;
fr=1/(2*3.14159*(Lr*Cr)^0.5);
Omega_o = 1/sqrt(Lr * Cr);
Co = 50E-6
n=5;
Kf= -8*Vin*Lm/(%pi*n*Lr*fr);
Rac=(8/%pi^2)*(n^2)*Rload;

//~~~~~~~~~~~~~~~~~~~~
//  Equivalent circuit model Ref.Paper [Tian2020]
Ln = Lm/Lr;  // Normalized inductance
Q= sqrt(Lr/Cr)/(n^2*Rload); // Load coefficient
fs=1.1E6;
Omega_n = fs/fr;  // Normalized frequency

//  DC gain calculation.  NO resistance considered.
if Omega_n >= 1 then
    alpha = %pi
else
    alpha = Omega_n * %pi
end
M1 =(1/sin(alpha/2));  // DC gain First part
M2 = %i * Omega_n * Ln / (%i * Omega_n * (Ln + 1 - 1/Omega_n^2 ) + (%pi^2/8) * (Q/sin(alpha/2)) * (1- Omega_n^2)*Ln * Omega_n );
M = abs( M1 * M2 ); 

//=============== Control to Vo Transfer function =====
Gdc1 = (Vin/(2*n))* (Ln/(Omega_o * Omega_n));
Gdc2 = (1/Omega_n^2 - Omega_n^2)*(%pi^2*Q*Ln/8)^2 - (Ln+1 - 1/Omega_n^2)*(2/Omega_n^2)
Gdc3 = 1/ (sqrt( (Ln+1-1/Omega_n^2)^2 +  ((1/Omega_n -Omega_n)* (%pi^2/8)*Q *Ln)^2 ))^3;
Gdc = Gdc1 * Gdc2 * Gdc3;

Leq=Lr*(3.14159^2/4);
Weq=n/(Leq*Co)^0.5;
Feq=Weq/(2*3.14159);
Gif_Wz1 = 1/(Co*Rload);
Gif_Kdc = Kf/Rload;
s=%s;
num = -Kf*(1/Rload + s*Co); // use positive Kf here
den = 1 + s* Leq/(n^2*Rload) + s*s*Leq*Co/n^2;
Gif= syslin('c', num,den); // Gid small-signal model
Q_Gif = 0.5 * (1/n/Rload)*sqrt(Leq/Co);
// ---- Set PI controller Gi
Kp1 = -10 * 1/Gif_Kdc;
Kp1 = 1;
T1= 1/(0.1*Weq*Kp1);
Gi= syslin('c', Kp1 * T1 *s +1,  T1*s);
//bode([Gif; Gif*Gc1] ,1E-3, 1E8, ["Gif", "open loop"]); // Gid bode plot
T2= 1/(2*3.14*10E3);
//T2=Rload*Co;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   get the esr resistor in as power loss in circuit
Z_Leq = syslin('c', s*Leq, n^2);
Z_Co=syslin('c', 1, s*Co);
Z_Rload = syslin('c', Rload, 1);
esr = 10e-3;
Z_esr = syslin('c', esr,1);
Z_Vin = syslin('c', Kf/n, 1);

Gif_esr = Z_Vin/(Z_Leq + Z_esr+ Z_Rload*Z_Co/(Z_Rload + Z_Co));
//-----set the LPF for current sense
Fp2 = 1000 // 1kHz pole
Wp2 = 2*3.14*Fp2;
//Wp2=63;
T2 =1/Wp2;
G_current_sense = syslin('c', 1, 1+ s * T2);  // LPF for current sense 1kHz Wp1
//--------------------------------

Gopen = Gi * Gif * G_current_sense;
Gi= syslin('c', 1+s/(2*3.14*50), s*(1+s/(2*3.14*50E3)));
clf; 
bode(Gif,1E-3, 1E8); // Gid bode plot

//---------Get the value from Analog control card 1900-7920 opamp N4B
//----------current controller
R53=47E3; R64=1E6; C42=2.2E-9; C45=2.2E-9;
Wz1_N4B= (R53+R64)/(2*R53*R64*C42); Fz1_N4B = Wz1_N4B/2/3.14;
Wp1_N4B = 1/(R53*C42); Wp2_N4B = 1/(R64*C45); 
Fp1_N4B = Wp1_N4B/2/3.14; Fp2_N4B= Wp2_N4B/2/3.14;
Z1_N4B = syslin('c', 2*R53*R64*C42*s + R53 +R64, (R53*C42*s +1)*(R64*C45*s+1));

//  Get the value from current sense filter
C43= 100E-9;  C69 = 100E-9;
R115=1000;
Wp1_CT_sense = 1/(R115*(C43+C69));
G_CT_sense = syslin('c', 1, 1+ s* R115*(C43+C69));

Gopen_analog = Gif * Z1_N4B * G_CT_sense;

printf('Fz1_N4B=%f \n',Fz1_N4B);
printf('Fp1_N4B=%f \n',Fp1_N4B);
printf('Fp2_N4B=%f \n',Fp2_N4B);
