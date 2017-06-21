%cagm is the AGM algoirthm for complex elliptic regions
%Given ph0 and kp.
k = sqrt(1 - kp^2);
Kp = ellipke(kp^2);
aa = 1; bb = k;cc = 1; 
ncag = 1;
phi = ph0;
            while cc > .0001
adp = atan(bb*tan(phi)/aa);
phig0 = phi + adp;
        if phig0 < 2*phi
    while abs(phig0 - 2*phi) > abs(phig0 + pi - 2*phi)
        phig0 = phig0 + pi;
    end
        else
    while abs(phig0 - 2*phi) > abs(phig0 - pi - 2*phi)
        phig0 = phig0 - pi;
    end
        end
phi = phig0;
        aa = (aa+bb)/2; bb = sqrt(aa*bb);
        cc = .5*(aa - bb);
        ncag = ncag + 1;
            end
dd = 2^(ncag-2)*pi;    
vv = phi/dd;
return
        