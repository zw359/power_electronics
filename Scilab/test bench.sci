funcprot(0);
exec('LLC_M_calc.sci', -1)

/*
Vbase = 0.5*Vin*M;
Vo = Vbase/N;
printf('Vo = %f\n', Vo);

*/

Vin=600;

Rload =0.1; // PN
//Rload =0.2; // PON
//Rload =1;   // PO
Rload =10;   //OPO

N = 16;
Cr =25E-9;
Lr = 47.0212E-6
Lm = 175.7023E-6;

N=5;
Cr=300E-9;
Lr = 30E-6;
Lm =226E-6;
Rload =1;

for i = 0.2 : 0.05:0.99
    Fn = i;
    M_final= LLC_M_calc(Cr, Lr, Lm, N, Rload,Fn );
    Vbase = 0.5*Vin*M_final;
    Vo = Vbase/N;
    printf('Vo = %f Fn = %f \n', Vo, Fn);
    //printf('%f\n', M_final);
end


