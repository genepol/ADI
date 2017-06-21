%ran5pt generates coefficients for 5-point 
%difference equations with Diff coef = 2*(minel + rand).
nH = input('number of columns nH =');
if isempty(nH)
    nH = 31;
end
nV = input('number of rows nV=');
if isempty(nV)
    nV = 16;
end
minel = input('minel is lower bound on elements');
if isempty(minel)
    minel = .001;
end
D5 = minel + rand(nV,nH);
%D5(n,m)is ONE HALF the diffusion coefficient in quadrant 1 at point (n,m).
cN = zeros(nV,nH);cE = zeros(nV,nH);cW = zeros(nV,nH);cS = zeros(nV,nH);cP = zeros(nV,nH);
%The diffusion equation coefficients are
%cN, cE, cW,cS, cP where
%for the interioe points:
    for n = 2:nV-1 %row n
        for m = 2:nH-1 %column m
            cN(n,m) = D5(n,m-1) + D5(n,m);
            cE(n,m) = D5(n,m) + D5(n-1,m);
            cW(n,m) = D5(n,m-1) + D5(n-1,m-1);
            cS(n,m) = D5(n-1,m-1) + D5(n-1,m);            
        end        
    end
%for the left and right edges:    
        for n = 2:nV
            cN(n,1) = 2*D5(n,1);
            cE(n,1) = D5(n,1) + D5(n-1,1);
            cW(n,1) = cE(n,1);
            cS(n,1) = 2*D5(n-1,1);
            cN(n,nH) = 2*D5(n,nH-1);
            cE(n,nH) = D5(n,nH-1) + D5(n-1,nH-1);
            cW(n,nH) = cE(n,nH);
            cS(n,nH) = 2*D5(n-1,nH);
        end
%for the top and bottome edges:        
        for m = 2:nH
            cE(1,m) = 2*D5(1,m);
            cN(1,m) = D5(1,m)+ D5(1,m-1);
            cW(1,m) = 2*D5(1,m-1);
            cS(1,m) = cN(1,m);
            cE(nV,m) = 2*D5(nV-1,m);
            cW(nV,m) = cE(nV,m-1);
            cS(nV,m) = D5(nV-1,m-1) + D5(nV-1,m);
            cN(nV,m) = cS(nV,m);
        end
%for the four corner points:
cN(1,1) = cS(2,1); cS(1,1) = cN(1,1); cN(1,nH) = cS(2,nH); 
cS(1,nH) = cN(1,nH); cS(nV,1) = cN(nV-1,1); cN(nV,1) = cS(nV,1);
cE(1,1) = cW(1,2); cW(1,1) = cE(1,1);
cE(nV,1) = cW(nV,2); 
cW(nV,1) = cE(nV,1); 
%The diagonal coefficient is:
cP = cE + cN + cS + cW;
%All five coefficients have been set for the nV*nH grid points.
D5 = 2*D5; %Half the D5 was used in the coefficient calculation.
%Now D5 is the fiffusion coefficient in quadrant 1 at column n and row m.
D6 = D5;
return                        