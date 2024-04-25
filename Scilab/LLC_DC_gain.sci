function [M]=LLC_DC_gain(Cr, Lr, Lm, N, Rload,Fn)
    h= Lm/Lr;
    Rac = (8/%pi^2) * N^2 * Rload;
    Q=sqrt(Lr/Cr)/Rac;
    mx=1+h;
    //M= (Fn^2*h)/sqrt((mx*Fn^2-1)^2 + Fn^2 *(Fn^2-1)^2 * h^2 *Q^2 );
    M= (Fn.^2*h)./sqrt((mx*Fn.^2-1).^2 + Fn.^2 .*(Fn.^2-1).^2 * h^2 *Q^2 );
endfunction
