%embed1 generates an elliptic function region containing
%a given N-stable spectrum Z.
reZ = real(Z);
aZ = min(reZ); bZ = max(reZ);
N = length(Z);
m = 1;
    for nn = 1:N
        if imag(Z(nn)) >= 0
            Z(m) = Z(nn);
            m = m+1;
        end
    end
    M = m - 1;
Z = Z(1:M);
rtabZ = sqrt(aZ*bZ);
Z1 = Z/rtabZ;
sZ = sqrt(bZ/aZ);
%The spectrum is normalized to Z1 with real interval (1/sZ,sZ).
hold off
logZ = log(Z1);
beta = max(imag(logZ));
rc = max(real(logZ)); 
mull = 1;
                        for jc = 1:5
        minreg = exp(-rc);
        maxreg = 1/minreg;
        cj(jc) = .1*jc;
        angA = beta + cj(jc)*(pi/2 - beta);
        angA = double(angA);
        angA = min(angA,(1+.5*cj(jc))*beta);
        angAj(jc) = angA;
%angA is always greater than the bound on the actual spectrum.
%The maximum angle at real(logZ) = 0 is angA.
%We now compute log(dn)-boundary values for the eigenvalues
rcj(jc) = rc;
dnbook
%We now find a new rc to enclose the eigenvalues.
ssp = rc;
       for m = 1:M
        xx = real(logZ(m));
        yy = imag(logZ(m));
    if yy > angA*sqrt(1 - xx^2/ssp^2)
        ssp = abs(xx)/sqrt(1 - (yy/angA)^2);
    end
       end
        rcj(jc) = ssp;
minreg = exp(-ssp);
maxreg = 1/minreg;
mul = 1;
dnbook
mulj(jc) = mul;
vagm
nits(jc) = mulj(jc)*K/(vv*KP);
qpjc(jc) = exp(-pi*K/KP);
kpjc(jc) = kpdn;
minregjc(jc) = minreg;
                        end
[nitbest ,jcbest] = min(nits);
kpbest = kpjc(jcbest);
bp = exp(rcj(jcbest));
ap = 1/bp;
angbest = angAj(jcbest);
return