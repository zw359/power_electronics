// DC gain calculation base FHA( fundamental hamonic analysis)
// NOT accurate compare to time-based analysis

function [M1, M2]=LLC_DC_gain(Cr, Lr, Lm, N, Rload,Fn)
    h= Lm/Lr;
    Rac = (8/%pi^2) * N^2 * Rload;
    Q=sqrt(Lr/Cr)/Rac;
    r=1/Q;
    mx=1+h;
    //M= (Fn^2*h)/sqrt((mx*Fn^2-1)^2 + Fn^2 *(Fn^2-1)^2 * h^2 *Q^2 );
    M1= (Fn.^2*h)./sqrt((mx*Fn.^2-1).^2 + Fn.^2 .*(Fn.^2-1).^2 * h^2 *Q^2 );
    M2=1/sqrt( (1+ 1/h - 1/(Fn.^2 *h)).^2 + (1/r^2)*(1/Fn -Fn).^2 );
endfunction
