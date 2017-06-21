%parsyl computes selected eigenvalues of S and SB.  
%ADI iteration parameters to attain a desired error reduction
%are then generated.
T = sparse(S);
n = length(S);
nlamin = ceil(2*sqrt(n));
nlamax = ceil(sqrt(nlamin));
opts.tol = 1e-6;
opts.disp = 0;
lamin = eigs(T,nlamin,'sr',opts);
lamax = eigs(T,nlamax,'lm',opts);
    if min(real(lamin)) <= 0
        error('Spectrum of A not in positive real half plane.');
    end
%We choose a conservative bound on the spectral radius:
b = 3*abs(lamax(1));
a = .8*min(real(lamin));
%We repeat the above for B:
TB = sparse(SB);
m = length(SB);
mlamin = ceil(2*sqrt(m));
mlamax = ceil(sqrt(mlamin));
Blamin = eigs(TB,mlamin,'sr',opts);
Blamax = eigs(TB,mlamax,'lm',opts);
    if min(real(Blamin)) <= 0
        error('Spectrum of B not in positive real half plane.');
    end
%We choose a conservative bound on the spectral radius:
d = 3*abs(Blamax(1));
c = .8*min(real(Blamin));
%We align the real eigenvalue intervals.
wbjreal
lamin = (del*lamin - bet)./(alp - gam*lamin);
Blamin = (del*Blamin + bet)./(alp + gam*Blamin);
lamax = (del*lamax - bet)./(alp - gam*lamax);
Blamax = (del*Blamax + bet)./(alp + gam*Blamax);
Aang = angle(lamin);
angA = max(Aang);
Bang = angle(Blamin);
angB = max(Bang);
%The real intercepts of both spectra are [kp,1].
%Eigenvalues at angles >= pi/3  are chosen
%as iteration parameters and excluded from 
%eigenvalues enclosed by the elliptic function region.
%This guarantees residual spectra with mm > 1.
j1 = 0;
J = 0;
ka = 1;
j1s = 0; j1t = 0;
                if angA >= pi/3
	for j = 1:nlamin
		if Aang(j) >= pi/3
			wr(ka) = lamin(j);
            wr(ka+1) = conj(lamin(j));
            ka = ka+2;
            lamin(j) = conj(lamin(j));
		end
    end
j1 = ka - 1;
    if j1 ~= 0
        wr = wr(1:j1);
%We back transform to obtain:
        ws = (alp*wr + bet)./(del + gam*wr);
        wt = ws;
    end
                end
Aang = angle(lamin);                
angA = max(Aang);
	for j = 1:mlamin
        if Bang(j) >= pi/3
			wr(ka) = Blamin(j);
            wr(ka+1) = conj(Blamin(j));
            ka = ka+2;
            Blamin(j) = conj(Blamin(j));
		end
    end
        if ka - 1 ~= j1
            wr = wr(1:ka-1);
%We back transform to obtain            
        wt(1+j1:ka-1) = (alp*wr(1+j1:ka-1) - bet)./(del - gam*wr(1+j1:ka-1));
        ws(1+j1:ka-1) = wt(1+j1:ka-1);
        end
j1 = ka -1;
    if j1 ~= 0
ws = ws(1:j1); wt = wt(1:j1);
    end
Bang = angle(Blamin);
angB = max(Bang);
angC = max(angA,angB);
ap = kp; bp = 1;
cp = kp; dp = 1;
                                    if angC < .1
%The spectra are treated as real.
q2= errY^2*(1+errY^2/4)^2/16;
    if kp^2 < .5
        qp = (kp/4)*(1+kp^2/4);
        qp = qp^2;
        q = exp(pi^2/log(qp));
    else
        ksq = 1 - kp^2;
        q = ksq*(1+ksq/4)^2/16;
        qp = exp(pi^2/log(q));
    end
J = ceil(.025*log(q2)/log(q));
Jli = J;
rtkp = sqrt(kp);
    for j = 1:J
        r= (2*j-1)/(2*J);
        nw = 1 + qp^(1-r) + qp^(1+r);
        dw = 1 + qp^r + qp^(2-r);
        xp = (2*r-1)/4;
        qpr = qp^xp;
        wJ(j) = qpr*nw/dw;
    end
    wJ = wJ(1:J);
%We now remove the WBJ transformation:
	wws = (alp*wJ - bet)./(del - gam*wJ);	                
	wwt = (alp*wJ + bet)./(del + gam*wJ);
        if j1 ~= 0
    ws = [wws ws]; wt = [wwt wt];
        end
                                    else
%The spectra are complex and must be aligned.                                        
csqa = cos(angA)^2;
csqb = 4/(2+kp+1/kp);
mmA = 2*csqa/csqb - 1;
csqa = cos(angB)^2;
mmB = 2*csqa/csqb - 1;
%By restricting to angC < pi/3 we assure mm > 1.
c11 = cos(angA);
c22 = cos(angB);
    if c11 > c22
        c1 = c11; c2 = c22;
    else
        c1 = c22; c2 = c11;
    end
e = sqrt(kp);    
f0 = (e+1/e)/2;
mu = (f0 - 1/f0)/2;
a3 = (mu*(c1 + c2) - (1 - c1*c2))/(c1-c2);
tau = a3 + sign(a3)*sqrt(a3^2-1);
betc = 1/(tau + sign(tau)*sqrt(tau^2-1));
c1p = (c1*tau + 1)/(tau + c1);
c2p = (c2*tau - 1)/(tau - c2);
A1 = acos(c1p);
A2 = acos(c2p);
f1 = (f0*tau + 1)/(tau + f0);
f2 = (f0*tau - 1)/(tau - f0);
zets = f1*c1p;
zett = f2*c2p;
    if abs(zets-zett) > .01
        error('alignment incorrect')
    end
ali = (f1 - sqrt(f1^2 - 1))/2; bli = 1/ali;
cli = (f2 - sqrt(f2^2 - 1))/2; dli = 1/cli;
mm = 2*zets^2 - 1;
%The spectra are aligned.
%The error reduction errY now determines J.
%Since MATLAB ellipj computes elliptic integrals
%we use sn(z) = sin(phi) and F(phi) to compute vv.
kp = mm + sqrt(mm^2-1);
kp = 1/kp; 
kp2 = kp^2;
sig = ap/(bp*kp);
fl = 1 - ap*kp/bp;
el = 1 - ap/(bp*kp);
snrK = sqrt(el/fl);
ph0 = asin(snrK);
cagm
y = .5*vv/Kp;
vJ = 1 - 2*y;
zz = (1 - sqrt(k))/(2*(1+sqrt(k)));
    if zz^2 < .5
        qp = zz*(1+zz^4);
        q = exp(pi^2/log(qp));
    else
        zp = sqrt(1-zz^2);
        q = zp*(1+zp^4);
        qp = exp(pi^2/log(q));
    end
Jli = ceil(log(errY/40)/(2*vJ*log(q)));
wli = zeros(1,Jli);    
for j = 1:Jli
	rr = (2*j-1)/(2*Jli);
	qexp = (2*rr-1)/4;
	numer = 1 + qp^(1-rr) + qp^(1+rr);
	denom = 1 + qp^rr + qp^(2-rr);
	wli(j) = qp^qexp*numer/denom;
end
wli = wli(1:Jli);
%We now transform back to the spectrum before complex
%alignment:
wws = e*(wli - betc)./(1 - betc*wli);
wwt = e*(wli + betc)./(1 + betc*wli);
%We now remove the WBJ transformation:
wws = (alp*wws + bet)./(del + gam*wws);
wwt = (alp*wwt - bet)./(del - gam*wwt);
    if j1 == 0
        ws = wws;  wt = wwt;
    else        
        ws = [wws ws]; wt = [wwt wt];
    end
                                    end
disp('Now the total number of parameters is:')
J = length(ws) 
disp('The number of complex parameters is:')
j1                                    
usol = zeros(nA,mB);                                    
return