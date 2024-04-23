funcprot(0);
exec('LLC_M_calc_ini.sci', -1)

Vin=280;

Rload =0.1; // PN
//Rload =0.2; // PON
//Rload =1;   // PO
Rload =0.6;   //OPO

N = 16;
Cr =25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;
Gamma = %pi/0.5;
// Initial condition
// this is important to get the proper iniital condition
// after that the previous calculation result (x1) will be used as (x0) for the next calculation
x1_0=0
x2_0=0
x3_0=0
x4_0=0
x5_0=0
x6_0=0
x7_0=0
x8_0=0
x9_0=0
x10_0=0
x11_0=0
x12_0=0
//x13_0=1  
//x14_0=%pi * 0.5
//x15_0=%pi/2
x13_0=2
x14_0=%pi*0.3
x15_0=%pi*0.6 
//x13_0=1.2
//x14_0=Gamma*0.3
//x15_0=Gamma*0.6 

x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0; x7_0; x8_0; x9_0; x10_0; x11_0; x12_0; x13_0; x14_0; x15_0];


M_0=1.2
alpha_0=0.5
beta_0=3

x=[]; y=[];
for i = 0.5 : 0.01:0.9
    Fn = i;
    //Fn = 0.499
    //pause
    [x1, correctMode]= LLC_M_calc_ini(Cr, Lr, Lm, N, Rload,Fn,x0 )
    Vbase = 0.5*Vin*x1(13);
    Vo = Vbase/N;
    printf('Vo = %f Fn = %f \n', Vo, Fn);
    x=[x i];
    y=[y Vo];
    x0=x1;
//    pause    
end

plot(x,y);

