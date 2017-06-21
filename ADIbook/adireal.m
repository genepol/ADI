%ADIREAL yields iteration parameters for
%the real two-variable ADI iteration 
H = input('Symmetric matrix H is:');
	if isempty(H)
        nH = input('order of H is:')
            if isempty(nH)
                nH = 100;
            end
        H = triu(ones(nH,nH),-1);
        H = -tril(H,1)+ 3*diag(diag(H));
%Default H is the tridiagonal matrix [-1,2,-1] of order 100. 
    end
F = input('Positive diagonal matrix F is:');
    if isempty(F)
        F = diag(diag(ones(nH)));
    else
        if length(F) ~= nH
            error('Inconsistent F input.')
        end
    end
%The default F is the identity matrix.    
F = sparse(F); H = sparse(H);
a = eigs(H,F,1,0);
b = eigs(H,F,1);
        nV = input('order of V is:')
            if isempty(nV)
                nV = 50;
            end
V = input('Symmetric matrix V is:');
	if isempty(V)
        V = triu(ones(nV,nV),-1);
        V = -tril(V,1)+ 3*diag(diag(V));
%Default V is the tridiagonal matrix [-1,2,-1] of order 50. 
    else
        if length(V) ~= nV
            error('Inconsistent V input.')
        end
    end            
G = input('Positive diagonal matrix G is:');
    if isempty(G)
        G = diag(diag(ones(nV)));        
    else
        if length(G) ~= nV
            error('Inconsistent G input.')
        end
    end
%The default G is the identity matrix.    
V = sparse(V); G = sparse(G);
c = eigs(V,G,1,0);
d = eigs(V,G,1);
    if a + c <= 0
        error('Spectra not in positive real plane')
    end      
wbjreal
%wbjreal returns the parameters kp,alp,bet,gam,del
%for spectral alignment
    eps = input('Desired error bound is:')
	if isempty(eps)
		eps = 1e-4
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
aditer
return