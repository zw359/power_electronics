Vin=50; // Full bridge
Ro=1;
Lr=360E-9;
Cr=110E-9;
Lm=2.1E-6;
fr=1/(2*3.14159*(Lr*Cr)^0.5);
Fo = fr;
Omega_o = 1/sqrt(Lr * Cr);
Wo=1/sqrt(Lr * Cr);
Co = 50E-6;
Np=10;
Ns=2;
n=Np/Ns;
Kf= -8*Vin*Lm/(%pi*n*Lr*fr);

Ln=Lm/Lr; // important resonant tank parameter
Zo=sqrt(Lr/Cr);
In= Vin/Zo;
Q=Zo/(n^2 * Ro);
Kvco = Fo;
Vm = Fs/Kvco;

Fs = 1.1E6; // swithing freq.
//Fs=1.01*Fo;
Fn= Fs/Fo;



//-------------- Equivalent circuit model [Hsien2020]
//  Fs > Fo  switching freq higher than resonant freq

Re=(Zo/n^2) * (Fn-1)^2 * (4+ 2*%pi*Q) / (1+ 4*(1 + 1/Ln)*(Fn-1));
Le= (%pi^2 /4) *(Lr/n^2) /(1 + (2 +5*Ln)*(Fn-1));
Ce=1/(Le* (2*%pi*(Fs-Fo))^2);
Yc=-(1/(Fn-1)^2) * (n*In/Vm) * (1/(2*Ln*(1+Q)) + Q*(Fn-1));
a1=-(3/(Ln*Q*(1+2*Q)))*(1+6*Q*(5+Ln)*(Fn-1)) * Wo;
a2=(1+(0.4-0.04*Ln)*Q + (0.24+0.13*Ln)*Q^2) *(1 +4*Q*Ln*(Fn-1)) * Wo^2 / (0.5 * (Ln-1)*Q -1.4);

//~~~~~~~~~~~~~~~~~~~~ control to vo function

s= %s;
Hz = syslin('c', 1+s/a1+s^2/a2, 1);
num = Yc *Re * Hz;
den = syslin('c', (1+Re/Ro) + s*(Le/Ro +(Co+Ce)*Re) + s^2 *( Le*(Co+Ce*Re/Ro)) +s^3 *Co*Ce*Le*Re, 1 ); 
den2= syslin('c', (1 + s*(Le/Ro + Co*Re) + s*s*Le*Co)*(1+s*Ce*Re), 1 )
Gfo = num/den;
Gfo1 = num/den2;
clf; 
bode([Gfo; Gfo1], 100, 1E6);
