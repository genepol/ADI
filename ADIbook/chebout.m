%chebout is the Chebyshev outer iteration
%for ADI preconditioned heat diffusion
dsam
%The sample problem is the heat diffusion equation over 
%a 16x31 grid with  x and y increments of unity snd 
%diffusion coefficient as on P. 55 in 5x10 blocks.
%The boundary condition is zero value at increments of
%unity from the boundary and symmetric diffusion coefficients 
%along the grid boundary. 
%dsam calls sepd which computes the iteration matrices
%Li= H,Gi=F,Lj=V,Gj=G for LiGj + LjGi = RS.
%The separable preconditioner yields an outer 
%iteration condition number of ~30 for exact ADI solution 
bout = max(max(D5./(F5'*E5)));
aout = min(min(D5./(F5'*E5)));
%sam5pt generates the five-point equations for the dsam problem.
sam5pt
precadi
%precadi computes the eigenvalue bounds a,b for
%(H - pF) and c,d, for (V - qG).
%Iteration parameters p(s) and q(s) for ADI iterations 
%are computed with use of wbjreal for spectral alignment.
%The error is reduced by a factor of eps
%so the Chebyshev outer iteration is with:
bout = (1+eps)*bout;
aout = (1-eps)*aout;
z = (bout+aout)/(bout-aout);
%The spectral radius for iteration without extrapolation is 1/z.
%The initial estimate is zero,
%The initial right-hand side is set to unity at all grid points
RS = ones(nV,nH);
%The ADI iterations are performed for the first outer iteration.
        vsol = zeros(nV,nH);
    for j = 1:J
        M = pj(j)*F + H;
        rhsj = (pj(j)*G - V)*vsol + RS;
        vsol = M\rhsj';
        rhsj = ((qj(j)*F - H)*vsol)' + RS;
        M = qj(j)*G + V;
        vsol = M\rhsj;
    end
usol = F\vsol';
usol = usol';
alph = 2/(aout+bout);
usol = alph*usol;
delold = usol;
uold = usol;
errsol = input('The solution accuracy boundis:')
    if isempty(errsol)
        errsol = 1e-4
    end
Jout = ceil(acosh(1/errsol)/acosh(z));   
%The error reduction estimate is precise so Jout 
%outer iterations suffice.
rout = z;
        for jout = 2:Jout
%The right hand side must be generated for each outer iteration
Auold = zeros(nV,nH);
for m = 1:nV %row m
    for n = 1:nH %column n
        Auold(m,n) = cP(m,n)*uold(m,n);
        if n < nH
            Auold(m,n) = Auold(m,n) - cE(m,n)*uold(m,n+1);
        end
        if n > 1
            Auold(m,n) = Auold(m,n) - cW(m,n)*uold(m,n-1);
        end
        if m < nV
            Auold(m,n) = Auold(m,n) - cN(m,n)*uold(m+1,n);
        end
        if m > 1
            Auold(m,n) = Auold(m,n) - cS(m,n)*uold(m-1,n);
        end
    end
end
%The inner iterations are performed:
        vsol = zeros(nV,nH);
        RHS = RS - Auold;
    for j = 1:J
        M = pj(j)*F + H;
        rhsj = (pj(j)*G - V)*vsol+ RHS;
        vsol = M\rhsj';
        rhsj = ((qj(j)*F - H)*vsol)' + RHS;
        M = qj(j)*G + V;
        vsol = M\rhsj;
    end
usol = F\vsol';
delu = usol';
%The Chebyshev outer iteration is performed
rout1 = 2*z - 1/rout;  %rout1 = T_{j+1}/T_j}.
alphout = 4/(rout1*(bout-aout));
betout = 1/(rout1*rout);
rout = rout1;
delu = alphout*delu + betout*delold;
uold = uold + delu;
delold = delu;
        end
disp('The solution is uold:')
disp('The Value of norm(RS -  Auold,2)/norm(RS,2) is:')
norm(RS -  Auold,2)/norm(RS,2)
disp('[J = ADI inners per outer  Jout = outers]')
[J Jout]
return