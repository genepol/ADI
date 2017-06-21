%adilap performs the ADI iteration to solve
%the reduced Lyapunov equation.
    if isempty(usol)
        usol = zeros(n);
    end
%usol is the initial estimate.
    if isempty(CS) 
        error('no source term')
    end
        Tp = T';
    for j = 1:J
        rhsj = ww(j)*usol - usol*Tp+ CS;
        M = ww(j)*speye(n) + T;
        usol = M\rhsj;
        usp = usol';
        rhsj = ww(j)*usp - usp*Tp + CS;
        usol = (M\rhsj)';
    end
%usol is now Y.
Cest = T*usol + usol*T';
disp('Estimated ||error in Y||/||Y||')
Yerr = norm(Cest - CS,1)/norm(CS,1)
return            