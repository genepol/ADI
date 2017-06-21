%adilow is preceded by qglow and qgpar.  
%A low rank Lyapunov equation is solved with Li-White ADI iteration
[n,m] = size(CS);
Z = zeros(n,m);
Y = zeros(n);
w = ww(1);
rtw = sqrt(w+conj(w));
RS = rtw*CS; %RS is the right-hand side for j = 1.
M = T + w*speye(n);
Z = M\RS;
Y = real(Z*Z');
%Y is the result of the first ADI iteration.
            for j = 2:J
wj = ww(j);
M = T + wj*speye(n); 
rtwj = sqrt(wj+conj(wj));
rt = rtwj/rtw;
rtw = rtwj;
RS = rt*Z;   %This is the right-hand side for iteration j. 
M = T + wj*speye(n);
wpwj = conj(wj) + wj;
wp = conj(w) + wj;
w = wj;
Z = M\RS;
Z = RS - wp*Z;
%Z is the result of iteration j.
Y = Y + real(Z*Z');
            end
disp('Estimated ||error in Y||/||Y||')
errinY = norm(FCS - T*Y - Y*T',1)/norm(FCS,1)
return