%Given minreg = min(real(eig(A))), maxreg = max(real(eig(A))),  
%and angA = max(angle(eig(A))) at sqrt(minreg*maxreg)of matrix A, 
%dnbook computes the dn boundary of the specrum 'eig'.  8/24/13
rtab = sqrt(minreg*maxreg);
adb = minreg/maxreg;
rtadb = sqrt(adb);
cosB = 2/(rtadb + 1/rtadb);
cosA = cos(angA);
mdn = 2*(cosA/cosB)^2 - 1;
    if mdn >= 1
        mul = 1;
    else
        mul = 2;
        mdn = 2*(cosB/cosA)^2 - 1;
        rtab = 1;
        minreg = tan(pi/4 - angA/2);
        maxreg = 1/minreg;
        adb = minreg^2;
    end
kpdn = 1/(mdn + sqrt(mdn^2 - 1));    
mp = kpdn^2;
mm = 1 - mp;
%The dn region has modulus = sqrt(mm).
%We now compute log(dn)-boundary values for the eigenvalues
KP = ellipke(mp);
    if kpdn < exp(-5)
        K = 2*KP*log(4/kpdn)/pi;
    else
        K = ellipke(mm);
    end
adbdkp = adb/kpdn;
dk = 1 - kpdn*adb;
si = sqrt((1-adbdkp)/dk);
ci = sqrt(adbdkp*mm/dk);
di = sqrt(mm/dk);
Z2 = zeros(21,1);
dnfac = rtab/sqrt(kpdn);
        for j = 1:30
    u = (j-1)/60;
[sr,cr,dr] = ellipj(K*u,mm);
%We evaluate dn(K(u-i*r),mm).
Z2(j) = dnfac*(ci*dr*di + i*mm*sr*si*cr)/(ci^2 + mm*si^2*sr^2);
        end
Z2(32:61) = Z2(30:-1:1);
Z2(31) = 1;
logZ2 = log(Z2);
logZ2(1:30) = -real(logZ2(61:-1:32))+i*imag(logZ2(61:-1:32));
logZ2(31) = i*angA;
    if mul == 2
        W = 1./Z2;
        W = 2./(Z2 + W);
        W1 = sqrt(W.*W - 1);
        U = W + W1;
        V = W - W1;
        F = ones(61,1);
F(61) = exp(rcj(jc)); F(1) = -F(61);
F(2:30) = log(V(32:60));
F(32:60) = log(conj(U(2:30)));
F(31) = i*angA;
logZ2 = F;
    end 
rZ = real(logZ2);
iZ = imag(logZ2);
return