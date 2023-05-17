Vin=600;
Rload=18;
Lr=30E-6;
Cr=300E-9;
Lm=220E-6;
fr=1/(2*3.14159*(Lr*Cr)^0.5);
Co = 4*220E-6;
Np=10;
Ns=11;
n=Np/Ns;
Kf= -8*Vin*Lm/(3.14159*n*Lr*fr);
Rac=(8/%pi^2)*(Np^2/Ns^2)*Rload;
Q=sqrt(Lr/(Cr * Rac));
printf('Q= %f',Q);

Leq=Lr*(3.14159^2/4);
Weq=n/(Leq*Co)^0.5;
Feq=Weq/(2*3.14159);
Gif_Wz1 = 1/(Co*Rload);
Gif_Kdc = Kf/Rload;
s=%s;
num = -Kf*(1/Rload + s*Co); // use positive Kf here
den = 1 + s* Leq/(n^2*Rload) + s*s*Leq*Co/n^2;
Gif= syslin('c', num,den); // Gid small-signal model

// ---- Set PI controller Gi
Kp1 = -10 * 1/Gif_Kdc;
Kp1 = 1E-4;
T1= 1/(1*Weq*Kp1);
printf('Gi T1=%f \n',T1);
//T1= 1/(Gif_Wz1*Kp1);
Gi= syslin('c', Kp1 * T1 *s +1,  T1*s);
Gif = Gif * 83E3;
//bode([Gif; Gif*Gc1] ,1E-3, 1E8, ["Gif", "open loop"]); // Gid bode plot
T2= 1/(2*3.14*10E3);
//T2=Rload*Co;


//-----set the LPF for current sense
Fp2 = 1000 // 1kHz pole
Wp2 = 2*3.14*Fp2;
//Wp2=63;
T2 =1/Wp2;
G_current_sense = syslin('c', 1, 1+ s * T2);  // LPF for current sense 1kHz Wp1
//--------------------------------

Giopen = Gi * Gif * G_current_sense;
//Gi= syslin('c', 1+s/(2*3.14*50), s*(1+s/(2*3.14*50E3)));
clf;
//bode(Gopen,1E-3, 1E8); // Gid bode plot
Giclose = Giopen/(1+Giopen);
bode(Giclose, 1E-3, 1E8);//


//------------  Voltage loop Control----
Fz2=1E3;
Wz2=2*3.14*Fz2;
Kp2=0.01;
T2=1/(Wz2*Kp2);
Gv= syslin('c', Kp2 * T2 *s +1,  T2*s);
clf;
bode([Gv*Giclose], 1E-3, 1E8);
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
