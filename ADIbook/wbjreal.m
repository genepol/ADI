%wbjreal aligns the spectra for the two sweeps.
    if a+c <= 0
        disp('spectrum not in positive real plane')
        stop
    end
md = (a+c)*(b+d); mn = 2*(b-a)*(d-c);
m = mn/md;
rtm  = sqrt(m*(2+m)); 
kp = 1/(1+m+rtm);
sig = 2*(a+d)/(b+d);
alp = b*sig - a*(1 + kp);
bet = a*(1+kp) - b*sig*kp;
gam = sig - 1 - kp;
del = 1 + kp - sig*kp;
return   