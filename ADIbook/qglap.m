%QGLAP computes a banded upper Hessenberg matrix S similar to %snorm*A.
%The last nonzero element in row k is element k of sband.
%Columns are reduced successively with Householder transformations.
%After column k is reduced, row k1 (the last unreduced row beyond column 
%k+1)is reduced with a householder transformation beyond column k+2.  
%Element S(k1,k+2) then has magnitude 
%equal to the 2-norm of S(k1,k+2:n) prior to reduction. If 
%|S(k1,k+2)/S(k1,k+1)| < tol, S(k1,k+2) is reduced to zero with a gaussian
%transformation and k1 is increased to k1+1.
%A file GI is also generated for recovery of X = GI*Y*GI^T.
%AX +XA' = C is solved for X by calling qgpar and adilap.
%Eigenvalues of A must lie in the right half plane bounded away
%from the imaginary axis and C must be SPD.
A = input('Matrix A is:');
if isempty(A)
    randn('seed',0)
    A = randn(30) + 5.5*eye(30);;
end
nA = length(A);
Xtr = zeros(nA);
C = input('matrix C is:');
if isempty(C)   
    C = rand(nA); C = C+C';
    c1 = min(eig(C));
    C = C + (.1+abs(c1))*eye(nA);
%This yields an spd C.    
end
myzero = 1e-14;
n = nA; snorm = 1/sqrt(norm(A,1)*norm(A,inf)); 
k1 = 1; 
S = snorm*A;
%This normalizes S so its eigenvalue magnitudes are bounded by unity.
CS = snorm*C;
GI = eye(n);
band = linspace(n,n,n); 
tol = input('Bound on gaussian multipliers is:');
if isempty(tol)
    tol = 1e3
end
%The 1-norm of S is initially around 1.  It has been observed that
%The 1-norm of the banded upper Hessenberg matrix similar to the initial S
%is of order magnitude less than tol.  Eigenvalues of the reduced matrix are
%reasonably close to those of S even when tol is large.  
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
CS(:,kp1:n) = CS(:,kp1:n)*HHk;
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
        k1 = k1+1; band(k1) = kp1;
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
CS(:,kp2:n) = CS(:,kp2:n)*HHk;
GI(:,kp2:n) = GI(:,kp2:n)*HHk;
S(kp2:n,kp1:n) = HHk*S(kp2:n,kp1:n);
CS(kp2:n,:) = HHk*CS(kp2:n,:);
S(k1,kp3:n) = zeros(1,n-kp3+1); %To eliminate roundoff error.
band(k1) = kp2;
dm = S(k1,kp1);
           if abs(dm) > sqrt(myzero)
tolk = abs(S(k1,kp2)/dm);
			if tolk > tol
%Element S(k1,kp2) cannot be reduced with a gaussian transformation.
            else
%We reduce S(k1,kp2) to zero with a gaussian transformation.
w = S(k1,kp2)/dm;
S(k1+1:n,kp2) = S(k1+1:n,kp2) - w*S(k1+1:n,kp1);
CS(:,kp1) = CS(:,kp1) + w*CS(:,kp2);
S(k1,kp2) = 0;
S(kp1,kp1:n) = S(kp1,kp1:n) + w*S(kp2,kp1:n);
CS(kp1,:) = CS(kp1,:) + w*CS(kp2,:);
GI(:,kp2) = GI(:,kp2) - w*GI(:,kp1);
band(k1) = kp1; k1 = k1+1;
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
CS(:,kp1) = CS(:,kp1) + w*CS(:,kp2);
S(k1,kp2) = 0;
S(kp1,kp1:n) = S(kp1,kp1:n) + w*S(kp2,kp1:n);
CS(kp1,:) = CS(kp1,:) + w*CS(kp2,:);
GI(:,kp2) = GI(:,kp2) - w*GI(:,kp1);
band(k1) = kp1; k1 = k1+1;
            end	
           end        
        end
k = n-1; 
sband = band;
%*******************S is banded upper Hessenberg****************
disp('1-norm of S is')
norm(S,1)
figure
spy(S);
hold on
title('Figure 1:  Profile of matrix S')
ylabel('matrix row'); xlabel('nonzero elements')
hold off
disp('number of elements beyond tridiagonal')
xtra = nnz(triu(S,2))
%ADI iteration parameters are determined:
errY = input('allowed error in Y')
    if isempty(errY)
        errY = 1e-6
    end
qgpar
usol = zeros(n);
adilap
    if Yerr > 10*errY
%Another set of ADI iterations is performed.        
        adilap
    end
X = GI*usol*GI';
disp('The estimated value for ||error in X||/||X|| is')
norm(C-A*X-X*A',1)/norm(C,1)
    if Xtr ~= zeros(nA)
disp('The true value for ||error in X||/||X|| is')
norm(X - Xtr,1)/norm(Xtr,1)
    end
return