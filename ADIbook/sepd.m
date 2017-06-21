%sepd computes E5 and F5 to approximate D5 with F5'*E5.
%The number of columns is nH and rows is nV.
%The tridiagonal Li and Lj matrices are generated. 
%Program generated on 1/2/2012.
E5 = ones(1,nH);
Li = zeros(nH);
Lj = zeros(nV);
Flast = ones(1,nV);
r1 = ones(1,nV);
c1 = ones(1,nH);
Drat = ones(nV,nH);
F5 = Flast;
k = 1;
    for j = 1:nV %row j
        cj1 = max(D5(j,:)); cj2 = min(D5(j,:));
        F5(j) = sqrt(cj1*cj2);       
    end
frat = 1;    
            while frat > .001
    for i = 1:nH %col i
        for j = 1:nV %row j
            r1(j) = D5(j,i)/F5(j);
        end
        ri1 = max(r1); ri2 = min(r1);
        E5(i) = sqrt(ri1*ri2);
    end
    for j = 1:nV %row j
        for i = 1:nH %col i
            c1(i) = D5(j,i)/E5(i);
        end
        cj1 = max(c1); cj2 = min(c1);
        F5(j) = sqrt(cj1*cj2);
    end
    frat = norm(Flast - F5)/norm(F5);
    Flast = F5;
            end
                for i = 1:nH %col i
                    for j = 1:nV %row j
                        Drat(j,i) = D5(j,i)/(E5(i)*F5(j));
                    end
                end
maxrat = max(max(Drat));
minrat = min(min(Drat));
pprec = maxrat/minrat;
mamD5 = max(max(D5)); minD5 = min(min(D5));
pD5 = mamD5/minD5;
benef = pD5/pprec;
    for ic = 2:nH-1
        Li(ic,ic+1) = -E5(ic);
        Li(ic,ic-1) = -E5(ic-1);
        Li(ic,ic) = E5(ic) + E5(ic-1);
    end
    Li(1,2) = -E5(1); Li(1,1) = 2*E5(1);
    Li(nH,nH) = E5(nH) + E5(nH-1);
    Li(nH,nH-1) = -E5(nH-1);
H = sparse(Li);
F = sparse(diag(.5*diag(H)));
    for jc = 2:nV-1
        Lj(jc,jc+1) = -F5(jc);
        Lj(jc,jc-1) = -F5(jc-1);
        Lj(jc,jc) = F5(jc) + F5(jc-1);
    end
    Lj(1,2) = -F5(1); Lj(1,1) = 2*F5(1);
    Lj(nV,nV) = F5(nV) + F5(nV-1);
    Lj(nV,nV-1) = -F5(nV-1);
V = sparse(Lj);
G = sparse(diag(.5*diag(V)));
return