%vagm is the AGM algoirthm for error reduction v
%Given ph0 and kpdn, K and KP.
bb = sqrt(1 - kpdn^2);
aa = 1; cc = 1; 
ncag = 1;
phi = asin(sqrt(adbdkp));
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
        