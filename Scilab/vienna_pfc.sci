Tline = 20E-3; // line frequency
Wline = 2* %pi / Tline;
t= [0: 10E-6 :Tline];
Vin = 230;
Po = 17.5E3;
Iin = Po/3/Vin;

// Modulation Index 
M=1.1;
M=2/sqrt(3);

Vdc = Vin*sqrt(2)*2/M;
V1 = Vdc/2;
V2 = Vdc/2;
C1 = 1880E-6;
C2 = 1880E-6;

// AC source voltage
va=sqrt(2)*Vin*sin(Wline * t);
vb=sqrt(2)*Vin*sin(Wline * t + %pi*2/3);
vc=sqrt(2)*Vin*sin(Wline * t - %pi*2/3);

// AC source current
//phase difference in degrees
phi = -0;
ia=sqrt(2)*Iin*sin(Wline * t + phi*%pi/180);
ib=sqrt(2)*Iin*sin(Wline * t + %pi*2/3 + phi*%pi/180);
ic=sqrt(2)*Iin*sin(Wline * t - %pi*2/3 + phi*%pi/180);

// open loop dutycyle
// derived as boost converter CCM Vin/Vo relation
// (2.14)
da=1 -abs(va)/V1;
db=1 -abs(vb)/V1;
dc=1 -abs(vc)/V1;

// middle point current from each phase leg
ina = ia .* da;
inb = ib .* db;
inc = ic .* dc;
inp = ina +  inb + inc;

// commom mode ref voltage
// control strategy
v0=0;

// currents into the centre point M
vam = va + v0;
vbm = vb + v0;
vcm = vc + v0;

ima = -(1/V1) * abs(vam) .* ia;
imb = -(1/V1) * abs(vbm) .* ib;
imc = -(1/V1) * abs(vcm) .* ic;

// ia + ib + ic = 0
// (2.15)
imp = ima +  imb + imc;
clf;
//plot(t, [va; vb; vc]);
//plot(t, [ia; ib; ic]);
//plot(t, [ imp; inp]);

vnp =  14/(3*Wline*2*C1);

// (2.22)
v0max_a = (Vdc/4) * (sign(ia) +1) -va;
v0max_b = (Vdc/4) * (sign(ib) +1) -vb;
v0max_c = (Vdc/4) * (sign(ic) +1) -vc;
v0max = min(v0max_a, v0max_b, v0max_c);


v0min_a = (Vdc/4) * (sign(ia) -1) -va;
v0min_b = (Vdc/4) * (sign(ib) -1) -vb;
v0min_c = (Vdc/4) * (sign(ic) -1) -vc;
v0min = max(v0min_a, v0min_b, v0min_c);
//plot(t, [v0max_a; v0max_b; v0max_c; v0max]);
clf;plot(t,[v0max; v0min; va; vb ; vc]);
