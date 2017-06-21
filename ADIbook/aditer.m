%aditer performs the ADI iteration
    if isempty(usol)
        vsol = zeros(nV,nH);
    else
        vsol = F*usol;
    end
%usol is the initial estimate.
    if isempty(RS) 
        RS = ones(nV,nH);
    end
%RS is the given right hand side.
    for j = 1:J
        rhsj = ((pj(j)*G - V)*vsol+ RS)';
        M = pj(j)*F + H;
        vsol = M\rhsj;
        rhsj = ((qj(j)*F - H)*vsol)' + RS;
        M = qj(j)*G + V;
        vsol = M\rhsj;
    end
usol = F\vsol';    
usol = usol';
return            