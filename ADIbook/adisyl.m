%adisyl performs the ADI iteration to solve
%the reduced Sylvester equation.
n = length(T); m = length(TB);
    if isempty(usol)
        usol = zeros(n,m);
    end
%usol is the initial estimate.
    if isempty(CS) 
        error('no source term')
    end
    for j = 1:J
        rhsj = ws(j)*usol - usol*TB + CS;
        M = ws(j)*speye(n) + T;
        usol = M\rhsj;
        usp = usol';
        rhsj = wt(j)*usp - usp*T' + CS';
        M = wt(j)*speye(m) + TB';
        usol = (M\rhsj)';
    end
usol = real(usol);    
Cest = T*usol + usol*TB;
disp('Estimated ||error in Y||/||Y||')
Yerr = norm(Cest - CS,1)/norm(CS,1)
Y = usol;
return            