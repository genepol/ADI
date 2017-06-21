%sam5pt generates 5-point difference equations 
%for the sample problem dsam.
%D5(n,m)is ONE HALF the diffusion coefficient in quadrant 1 at point (n,m).
D5 = D5/2;
cN = zeros(nV,nH);cE = zeros(nV,nH);cW = zeros(nV,nH);cS = zeros(nV,nH);cP = zeros(nV,nH);
%The diffusion equation coefficients are
%cN, cE, cW,cS, cP where (N,E,W,S = north,east,west south and p = center.
%For the interior points:
    for n = 2:nV-1 % row index
        for m = 2:nH-1 % column index          
            cN(n,m) = D5(n,m-1) + D5(n,m);
            cE(n,m) = D5(n,m) + D5(n-1,m);
            cW(n,m) = D5(n,m-1) + D5(n-1,m-1);
            cS(n,m) = D5(n-1,m-1) + D5(n-1,m);            
        end        
    end
%For the left and right edges:    
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
%For the top and bottom edges:        
        for m = 2:nH %col m
            cE(1,m) = 2*D5(1,m);
            cN(1,m) = D5(1,m)+ D5(1,m-1);
            cW(1,m) = 2*D5(1,m-1);
            cS(1,m) = cN(1,m);
            cE(nV,m) = 2*D5(nV-1,m);
            cW(nV,m) = cE(nV,m-1);
            cS(nV,m) = D5(nV-1,m-1) + D5(nV-1,m);
            cN(nV,m) = cS(nV,m);
        end
%For the corner points (1,1), and (nV,1):        
cN(1,1) = cS(2,1); cS(1,1) = cN(1,1); 
cS(nV,1) = cN(nV-1,1); cN(nV,1) = cS(nV,1);
cE(1,1) = cW(1,2); cW(1,1) = cE(1,1); 
cE(nV,1) = cW(nV,2); cW(nV,1) = cE(nV,1); 
%The diagonal coefficient is:
cP = cE + cN + cS + cW;
%All five coefficients have been set for the nV*nH grid points.
D6 = 2*D5; %Half the D5 was used in the coefficient calculation.
%Now D6 is the fiffusion coefficient in quadrant 1 at column m and row n.
D5 = .5*F5'*E5;
%Coefficients are now computed for the separable problem.
    for n = 2:nV-1 % row index
        for m = 2:nH-1 % ccolumn index          
            ccN(n,m) = D5(n,m-1) + D5(n,m);
            ccE(n,m) = D5(n,m) + D5(n-1,m);
            ccW(n,m) = D5(n,m-1) + D5(n-1,m-1);
            ccS(n,m) = D5(n-1,m-1) + D5(n-1,m);            
        end        
    end
%for the left and right edges:    
        for n = 2:nV
            ccN(n,1) = 2*D5(n,1);
            ccE(n,1) = D5(n,1) + D5(n-1,1);
            ccW(n,1) = ccE(n,1);
            ccS(n,1) = 2*D5(n-1,1);
            ccN(n,nH) = 2*D5(n,nH-1);
            ccE(n,nH) = D5(n,nH-1) + D5(n-1,nH-1);
            ccW(n,nH) = ccE(n,nH);
            ccS(n,nH) = 2*D5(n-1,nH);
        end
%for the top and bottom edges:        
        for m = 2:nH %ccol m
            ccE(1,m) = 2*D5(1,m);
            ccN(1,m) = D5(1,m)+ D5(1,m-1);
            ccW(1,m) = 2*D5(1,m-1);
            ccS(1,m) = ccN(1,m);
            ccE(nV,m) = 2*D5(nV-1,m);
            ccW(nV,m) = ccE(nV,m-1);
            ccS(nV,m) = D5(nV-1,m-1) + D5(nV-1,m);
            ccN(nV,m) = ccS(nV,m);
        end
%for the corner points (1,1) and (nV,1):        
ccN(1,1) = ccS(2,1); ccS(1,1) = ccN(1,1); 
ccS(nV,1) = ccN(nV-1,1); ccN(nV,1) = ccS(nV,1);
ccE(1,1) = ccW(1,2); ccW(1,1) = ccE(1,1); 
ccE(nV,1) = ccW(nV,2); ccW(nV,1) = ccE(nV,1); 
%The diagonal ccoefficient is:
ccP = ccE + ccN + ccS + ccW;
%All five separable coefficients have been set for the nV*nH grid points.
D5 = 2*D5;
%D5 are the diffusion coefficients for the separable problem.
%D6 are the diffusion coefficients for the actual problem.
return                        