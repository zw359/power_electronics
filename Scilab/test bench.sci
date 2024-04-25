funcprot(0);
exec('LLC_M_calc.sci', -1)

/*
Vbase = 0.5*Vin*M;
Vo = Vbase/N;
printf('Vo = %f\n', Vo);

*/

Vin=280;

Rload =0.1; // PN
//Rload =0.2; // PON
//Rload =1;   // PO
Rload =2;   //OPO

N = 16;
Cr =25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;

/*
N=5;
Cr=300E-9;
Lr = 30E-6;
Lm =226E-6;
Rload =1;
*/

M_0=2.2;
alpha_0=0;
beta1_0=1;


M_0=2.931317
alpha_0=3.434248
beta1_0=6.042374

x=[]; y=[];
for i = 0.62 : 0.01:0.9
    Fn = i;
//    Fn = 0.54
    [M_final, alpha_final, beta1_final]= LLC_M_calc(Cr, Lr, Lm, N, Rload,Fn);
    Vbase = 0.5*Vin*M_final;
    Vo = Vbase/N;
    printf('Vo = %f Fn = %f \n', Vo, Fn);
    x=[x i];
    y=[y Vo];
    alpha_0= alpha_final;
    beta1_0= beta1_final;
    M_0= M_final;
    
    printf('%f\n', M_final);
    printf('%f\n', alpha_0);
    printf('%f\n', beta1_0);
    pause    
end


