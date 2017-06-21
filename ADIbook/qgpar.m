%qgpar computes the selected ceil(sqrt(n)) small eigenvalues   
%and ceil(n^1/4)) large eigenvalues of S.
%ADI iteration parameters to attain a desired error reduction
%are then generated.
T = sparse(S);
nlamin = ceil(3*sqrt(n));
nlamax = ceil(sqrt(nlamin));
opts.tol = 1e-6;
opts.disp = 0;
lamin = eigs(T,nlamin,'sr',opts);
lamax = eigs(T,nlamax,'lm',opts);
j1 = 0;
wmag = 0;
if min(real(lamin)) <= 0
disp('Spectrum not in positive real half plane.');
disp('Hit [enter] to continue anyway and [ctrl-c] to end.')
pause
end
%We choose a conservative bound on the spectral radius:
b = 3*abs(lamax(1));
ang = angle(lamin);
angmax = max(ang);
ka = 1;
            if angmax >= pi/3
%Eigenvalues at angles greater than 60 deg are chosen
%as iteration parameters and excluded from 
%eigenvalues enclosed by the elliptic function region.
%This assures only real parameters for the residual regions.
	for j = 1:nlamin
		if ang(j) > pi/3
			ww(ka) = lamin(j);
            ww(ka+1) = conj(lamin(j));
            ka = ka+2;
            lamin(j) = b;
		end
	end
angmax = max(angle(lamin));
            end
rat = real(lamin(1))/b;
if rat < 1e-6
	disp('Poorly conditioned problem with eig ratio=')
rat
	disp('Hit [enter] to continue and [ctrl-c] to stop.')
pause
end
j = 1;
	while rat < 1e-6
		ww(ka) = lamin(j);
		lamin(j) = b;
		ka = ka + 1;
		j = j + 1;
		rat = abs(lamin(j))/b;
    end
%We choose a conservative lower bound on the real components:
a = 0.8*min(real(lamin));
ap = a;
bp = b;
rtab = sqrt(a*b);
    j1 = ka-1;
    %We normalize ww:    
	for j = 1:j1
	   ww(j) = ww(j)/rtab;
    end
csqa = cos(angmax)^2;
csqb = 4/(2+a/b+b/a);
mm = 2*csqa/csqb - 1;
%By choosing angmax < pi/3, we ensure mm > 1.
%For larger angles admitted.**********************************
    if mm < 1
%Complex parameters must be computed
%from a dual spectrum.
x1 = lamin/rtab;
x2 = rtab./lamin;
x12 = [x1;x2];
x3 = lamax/rtab;
x4 = rtab./lamax;
x34 = [x3;x4];
x13 = [x1;x3];
x24 = [x2;x4];
xx = 2./(x13 + x24);
z1 = xx + sqrt(xx.*xx - 1);
z2 = xx - sqrt(xx.*xx - 1);
z = [z1;z2];
z = real(z) + i*abs(imag(z));
mm = 2*csqb/csqa - 1;
ap = tan(pi/4 - angmax/2);
bp = 1/ap;
wmag = 1
    end
%Resume with mm >= 1*****************************************    
%The error reduction errY now determines J.
kp = mm + sqrt(mm^2-1);
kp = 1/kp;
kp2 = kp^2;
k = sqrt(1-kp2);
zz = (1 - sqrt(k))/(2*(1+sqrt(k)));
    if zz^2 < .5
        qp = zz*(1+zz^4);
        q = exp(pi^2/log(qp));
    else
        zp = sqrt(1-zz^2);
        q = zp*(1+zp^4);
        qp = exp(pi^2/log(q));
    end
        if angmax < .01
%When the angle is < .01 we compute the
%ADI error reduction for a real spectrum.
	Jl = ceil(0.5*log(errY/40)/log(q));
        else
%Since MATLAB ellipj computes elliptic integrals
%we use sn(z) = sin(phi) and F(phi) to compute vv.
rtkp = sqrt(kp);
sig = ap/(bp*kp);
fl = 1 - ap*kp/bp;
el = 1 - ap/(bp*kp);
snrK = sqrt(el/fl);
ph0 = asin(snrK);
Kp = ellipke(kp2);
cagm
y = .5*vv/Kp;
vJ = 1 - 2*y;
J1 = ceil(log(errY/40)/(2*vJ*log(q)));
        end
    if j1 ~= 0
        ww = [zeros(1,J1) ww(1:j1)];
    else
        ww = zeros(1,J1);
    end
for j = 1:J1
	rr = (2*j-1)/(2*J1);
	qexp = (2*rr-1)/4;
	numer = 1 + qp^(1-rr) + qp^(1+rr);
	denom = 1 + qp^rr + qp^(2-rr);
	ww(j) = qp^qexp*numer/denom;
end
        if wmag == 0 %Always for angmax < pi/3.
%We now denormalize ww:
            ww = rtab*ww;
        else
	JJ = floor(Jl/2);
	for j = 1:JJ
	 par = 2/(ww(j+j1) + 1/ww(j+j1));
	 angj = acos(par);
	 wwp(j) = rtab*exp(i*angj);
	end
	if j1 ~= 0
		ww(1:j1) = rtab*ww(1:j1);
	end
			for j = 1:JJ
                ww(2*j-1) = wwp(j);
                ww(2*j) = wwp(j)';
            end
          if Jl/2 ~= JJ
		    ww(Jl) = rtab;	
          end
        end
J = J1 + j1;
%The iteration parameters for ADI error
%reduction [terror] are the J values in ww(1:J)
%with the j1 discrete complex values at the end.
disp('qgpar ended with j1 and J = ')
j1
J
return