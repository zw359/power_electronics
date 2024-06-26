exec('PN_calc.sci', -1);
exec('PON_calc.sci', -1);

Vo=12;
Po=600;
Rload = Vo^2/Po;
Io = Vo/Rload;

Cr_step = 1E-9; // incremenal values.
Cr_init = 6E-9; // Inicial values.
Cr= Cr_init;

i=0;//  counter starts from zero

Vin_min = 280;
Fsw_min = 100E3;
N = 16;

PN=1
Vmargin=1
Fp=0;
K=1;
fd =mopen('result.txt', 'wt');
mfprintf(fd,'Cr(nF) \t Lr(uH) \t Lp(uH) \t  \Fr(kHz) \t Fp(kHz) \t Mode \t K \t Z0 \n');


while i<25
    Cr = Cr_init + i * Cr_step;
    if PN==1 then 
        [Lr, Lp, Fr, Fp, Vmargin, Lambda, Theta, K] = PN_calc(Vo, Rload, Vin_min, Fsw_min, N, Cr); 
        Lambda0=Lambda;
        Theta0=Theta;
        if Vmargin <=0 then 
            PN=0;
        else
            mfprintf(fd,'%.2f \t %.4f \t %.4f \t %.4f \t%.4f \t %s \t \%0.2f \t %f \n' , Cr*1E9, Lr*1E6, Lp*1E6, Fr, Fp, 'PN', K, sqrt(Lr/Cr)); 
            i = i+1;
        end
    else
        
        break
    end
end



while i<25
    Cr = Cr_init + i * Cr_step;
    if PN==0 then
        [Lr, Lp, Fr, Fp,Lambda, Theta, K] = PON_calc(Vo, Rload, Vin_min, Fsw_min, N, Cr, Theta0, Lambda0);
        Lambda0=Lambda;
        Theta0=Theta;
        mfprintf(fd,'%.2f \t %.4f \t %.4f \t %.4f \t%.4f \t %s \t \%0.2f \t %f \n' , Cr*1E9, Lr*1E6, Lp*1E6, Fr, Fp, 'PON', K, sqrt(Lr/Cr)); 
        i = i+1;
    end
end

            
mclose(fd); // close the result.txt file


