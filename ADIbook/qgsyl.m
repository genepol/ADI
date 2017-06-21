%qgsyl solves the Sylvester matrix equation with
%a full right-hand side.
A = input('matrix A is:');
    if isempty(A)
        randn('seed',0)
        A = randn(40) + 6.5*eye(40);
    end
nA = length(A);
B = input('matrix B = ');
    if isempty(B)
        randn('seed',0)
        B = randn(20) + 4.2*eye(20);
    end
mB = length(B);
Xtr = zeros(nA,mB);
C = input('Matrix C = ');
    if isempty(C)
        Xtr = rand(nA,mB);
        C = A*Xtr + Xtr*B;
%This sets the solution to Xtr.
    end    
  if size(C) ~= [nA mB]
        error('C wrong dimension.')
        return
  end
%We now compite the banded upper-Hessenberg matrix 
%S similar to A:
myzero = 1e-14;
n = nA; 
snorm = 1/sqrt(norm(A,1)*norm(A,inf)); 
tnorm = 1/sqrt(norm(B,1)*norm(B,inf)); 
snorm = sqrt(snorm*tnorm);
k1 = 1; 
S = snorm*A;
CS = snorm*C;
GI = eye(n);
tol = input('Bound on gaussian multipliers is:');
if isempty(tol)
    tol = 1e3
end
                    for k = 1:n-2
kp1 = k+1; kp2 = k+2;kp3 = min(k+3,n);
colvec = S(kp1:n,k); Y = max(abs(colvec));
%The column is reduced with a HH transformation:
vk = S(kp1:n,k); 
    if vk(1) == 0
        sgn = 1
    else
        sgn = sign(vk(1));
    end
u1 = sgn*norm(vk,2);
uk = [u1;zeros(n-kp1,1)]; wk = uk + vk; nmk = wk'*wk;
HHk = eye(n-k) - 2*(wk*wk')/nmk;
S(kp1:n,k:n) = HHk*S(kp1:n,k:n);
CS(kp1:n,:) = HHk*CS(kp1:n,:);
S(k1:n,kp1:n) = S(k1:n,kp1:n)*HHk;
S(kp2:n,k) = zeros(n-kp1,1); %To eliminate roundoff error.
%The column has been reduced.
%The matrix is banded above k.
GI(:,kp1:n) = GI(:,kp1:n)*HHk;
rowvec = S(k1,kp1:n); Z = max(abs(rowvec));
                            if k < n-2
nm = norm(S(k1,kp2:n),inf);
    if nm < myzero
%The row does not have to be reduced.        
        S(k1,kp2:n) = zeros(1,n-kp1);
        k1 = k1+1; 
    else
%Row k1 is reduced beyond kp2 with a HH transformation .
vk = S(k1,kp2:n); 
    if vk(1) == 0
        sgn = 1
    else
        sgn = sign(vk(1));
    end
u1 = sgn*norm(vk,2);
uk = [u1 zeros(1,n-kp2)]; wk = uk + vk; nmk = wk*wk';
HHk = eye(n-kp1) - 2*(wk'*wk)/nmk;
S(k1:n,kp2:n) = S(k1:n,kp2:n)*HHk;
GI(:,kp2:n) = GI(:,kp2:n)*HHk;
S(kp2:n,kp1:n) = HHk*S(kp2:n,kp1:n);
CS(kp2:n,:) = HHk*CS(kp2:n,:);
S(k1,kp3:n) = zeros(1,n-kp3+1); %To eliminate roundoff error.
dm = S(k1,kp1);
           if abs(dm) > sqrt(myzero)
tolk = abs(S(k1,kp2)/dm);
			if tolk > tol
%Element S(k1,kp2) cannot be reduced with a gaussian transformation.
            else
%We reduce S(k1,kp2) to zero with a gaussian transformation.
w = S(k1,kp2)/dm;
S(k1+1:n,kp2) = S(k1+1:n,kp2) - w*S(k1+1:n,kp1);
S(k1,kp2) = 0;
S(kp1,kp1:n) = S(kp1,kp1:n) + w*S(kp2,kp1:n);
CS(kp1,:) = CS(kp1,:) + w*CS(kp2,:);
GI(:,kp2) = GI(:,kp2) - w*GI(:,kp1);
k1 = k1+1;
            end	
           end        
    end      
                            end %on k < n-2
                    end %on k =1:n-2
        if k1 == n-2
dm = S(k1,kp1);
           if abs(dm) > sqrt(myzero)
tolk = abs(S(k1,kp2)/dm);
			if tolk > tol
%Element S(k1,kp2) cannot be reduced with a gaussian transformation.
            else
%We reduce S(k1,kp2) to zero with a gaussian transformation.
w = S(k1,kp2)/dm;
S(k1+1:n,kp2) = S(k1+1:n,kp2) - w*S(k1+1:n,kp1);
S(k1,kp2) = 0;
S(kp1,kp1:n) = S(kp1,kp1:n) + w*S(kp2,kp1:n);
CS(kp1,:) = CS(kp1,:) + w*CS(kp2,:);
GI(:,kp2) = GI(:,kp2) - w*GI(:,kp1);
k1 = k1+1;
            end	
           end        
        end
k = n-1;
%S is snorm*G*A*G^-1, CS = snorm*G*C, GI is G^-1.
disp('1-norm of S is')
norm(S,1)
disp('number of elements beyond tridiagonal')
xtra = nnz(triu(S,2))
%We now reduce B to banded upper Hessenberg form SB.
m = mB; 
k1 = 1; 
H = eye(m);
SB = snorm*B;
                    for k = 1:m-2
kp1 = k+1; kp2 = k+2;kp3 = min(k+3,m);
colvec = SB(kp1:m,k); Y = max(abs(colvec));
%The column is reduced with a HH transformation:
vk = SB(kp1:m,k); 
    if vk(1) == 0
        sgn = 1
    else
        sgn = sign(vk(1));
    end
u1 = sgn*norm(vk,2);
uk = [u1;zeros(m-kp1,1)]; wk = uk + vk; nmk = wk'*wk;
HHk = eye(m-k) - 2*(wk*wk')/nmk;
SB(kp1:m,k:m) = HHk*SB(kp1:m,k:m);
SB(k1:m,kp1:m) = SB(k1:m,kp1:m)*HHk;
SB(kp2:m,k) = zeros(m-kp1,1); %To eliminate roundoff error.
CS(:,kp1:m) = CS(:,kp1:m)*HHk;
H(kp1:m,:) = HHk*H(kp1:m,:);
%The column has been reduced.
%The matrix is banded above k.
rowvec = SB(k1,kp1:m); Z = max(abs(rowvec));
                            if k < m-2
nm = norm(SB(k1,kp2:m),inf);
    if nm < myzero
%The row does not have to be reduced.        
        SB(k1,kp2:m) = zeros(1,m-kp1);
        k1 = k1+1; 
    else
%Row k1 is reduced beyond kp2 with a HH transformation .
vk = SB(k1,kp2:m); 
    if vk(1) == 0
        sgn = 1
    else
        sgn = sign(vk(1));
    end
u1 = sgn*norm(vk,2);
uk = [u1 zeros(1,m-kp2)]; wk = uk + vk; nmk = wk*wk';
HHk = eye(m-kp1) - 2*(wk'*wk)/nmk;
SB(k1:m,kp2:m) = SB(k1:m,kp2:m)*HHk;
CS(:,kp2:m) = CS(:,kp2:m)*HHk;
H(kp2:m,:) = HHk*H(kp2:m,:);
SB(kp2:m,kp1:m) = HHk*SB(kp2:m,kp1:m);
SB(k1,kp3:m) = zeros(1,m-kp3+1); %To eliminate roundoff error.
dm = SB(k1,kp1);
           if abs(dm) > sqrt(myzero)
tolk = abs(SB(k1,kp2)/dm);
			if tolk > tol
%Element SB(k1,kp2) cannot be reduced with a gaussian transformation.
            else
%We reduce SB(k1,kp2) to zero with a gaussian transformation.
w = SB(k1,kp2)/dm;
SB(k1+1:m,kp2) = SB(k1+1:m,kp2) - w*SB(k1+1:m,kp1);
SB(k1,kp2) = 0;
SB(kp1,kp1:m) = SB(kp1,kp1:m) + w*SB(kp2,kp1:m);
CS(:,kp2) = CS(:,kp2) - w*CS(:,kp1);
H(kp1,:) = H(kp1,:) + w*H(kp2,:);
k1 = k1+1;
            end	
           end        
    end      
                            end %on k < m-2
                    end %on k =1:m-2
        if k1 == m-2
dm = SB(k1,kp1);
           if abs(dm) > sqrt(myzero)
tolk = abs(SB(k1,kp2)/dm);
			if tolk > tol
%Element SB(k1,kp2) cannot be reduced with a gaussian transformation.
            else
%We reduce SB(k1,kp2) to zero with a gaussian transformation.
w = SB(k1,kp2)/dm;
SB(k1+1:m,kp2) = SB(k1+1:m,kp2) - w*SB(k1+1:m,kp1);
SB(k1,kp2) = 0;
SB(kp1,kp1:m) = SB(kp1,kp1:m) + w*SB(kp2,kp1:m);
CS(:,kp2) = CS(:,kp2) - w*CS(:,kp1);
H(kp1,:) = H(kp1,:) + w*H(kp2,:);
k1 = k1+1;
            end	
           end        
        end
k = m-1; 
%SB is snorm*H*B*H^-1, and CS = snorm*G*C*H^-1
disp('1-norm of SB is')
norm(SB,1)
disp('number of elements beyond tridiagonal')
xtra = nnz(triu(SB,2))
errY = input('errY =')
    if isempty(errY)
        errY = 1e-6
    end
%We now compute the ADI iteration parameters.    
parsyl
%We now iterate
adisyl
    if Yerr > 5*errY
%We cycle through the iterations one more time.        
        adisyl
    end        
	X = GI*Y*H;
	Z = A*X;
	W = X*B;
	Xerr = norm(C - Z - W,1)/norm(C,1);
        if Xerr > 10*errY
%We cycle through the iterations one more time.
            adisyl
	X = GI*Y*H;
	Z = A*X;
	W = X*B;
        end
disp('The estimated error is:')    
	Xerr = norm(C - Z - W,1)/norm(C,1)    
        if Xtr ~= zeros(nA,mB)
            disp('The true error is:')
            norm(X-Xtr,1)/norm(Xtr,1)
        end
return