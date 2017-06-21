%precadi generates the precondition ADI pj and qj.
%H,F,V, and G are given as sparse matrices.
nH = length(H); nV = length(V);
opts.tol = 1e-6;
opts.disp = 0;
a = eigs(H,F,1,'sa',opts);
b = eigs(H,F,1,'la',opts);
c = eigs(V,G,1,'sa',opts);
d = eigs(V,G,1,'la',opts);
    if a + c <= 0
        error('Spectra not in positive real plane')
    end 
wbjreal
%wbjreal returns the parameters kp,alp,bet,gam,del
%for spectral alignment
    eps = input('Desired error bound is:')
	if isempty(eps)
		eps = .01
	end
q2= eps^2*(1+eps^2/4)^2/16; 
qp= kp^2*(1+kp^2/4)^2/16;
J = ceil(.25*log(q2)*log(qp)/pi^2);
rtkp = sqrt(kp);
ww = zeros(1,J);
    for j = 1:J
        r= (2*j-1)/(2*J);
        nw = 1 + qp^(1-r) + qp^(1+r);
        dw = 1 + qp^r + qp^(2-r);
        xp = (2*r-1)/4;
        qpr = qp^xp;
        ww(j) = rtkp*qpr*nw/dw;
    end
%ww(1:J) are the ADI parameters for the aligned spectra.
pj = (alp*ww - bet)./(del - gam*ww);
qj = (alp*ww + bet)./(del + gam*ww);
return