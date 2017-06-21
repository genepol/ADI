%qglow solves for X the Lyapunov equation AX + XA' = C where eigenvalues 
%of A (of order n) are in the right half plane bounded away from the imaginary axis. 
%C = cl*cl' with cl of order nxm with m<n.
%A is reduced with a similarity transformation to a banded upper 
%Hessenberg matrix S = snorm*G*A*inv(G). A matrix CS = snorm*G*cl is also generated.
%The value of snorm is the reciprocal of the square root 
%of the product of the 1-norm and infinity-norm of matrix A. This bounds magnitudes 
%of eigenvalues of S by unity. The reduced Lyapunov equation is SY + YS' = B where
%B = CS*CS'. Eigenvalues of S are determined with eigs. ADI iteration parameters are 
%computed with qgpar  The Li-White ADI iteration is performed with adilow
%to achieve a prescribed accuracy.  The inverse GI of G is also generated for 
%recovery of X = GI*Y*GI'/snorm.  
%This m-file was generated on 12/31/12.
randn('seed',0)
A = input('Matrix A is ');
if isempty(A)
disp('A = randn(30) + s*eye(30) to assure N-stability')
A = randn(30);
A = A + (.3 + abs(min(real(eig(A)))))*eye(30);
%We require an N-stable matrix A.
end
n = length(A);
cl = input('Matrix cl is ');
if isempty(cl)
disp('cl = rand(n,2)')
cl = rand(n,2);
end
C = cl*cl';
if length(C) ~= n
disp('size of cl and A not consistent')
return
end
snorm = 1/sqrt(norm(A,1)*norm(A,inf)); 
k1 = 1; 
S = snorm*A;
%This normalizes S so its eigenvalue magnitudes are bounded by unity.
CS = sqrt(snorm)*cl;
GI = eye(n);
band = linspace(n,n,n); 
tol = input('Bound on gaussian multipliers is:');
if isempty(tol)
    tol = 1e3
end
%The 1-norm of S is initially around 1.  Eigenvalues of the reduced 
%matrix are reasonably close to those of S even when tol is large.  
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
S(k1:n,kp1:n) = S(k1:n,kp1:n)*HHk;
CS(kp1:n,:) = HHk*CS(kp1:n,:);
S(kp2:n,k) = zeros(n-kp1,1); %To eliminate roundoff error.
%The column has been reduced.  The matrix is banded above k.
GI(:,kp1:n) = GI(:,kp1:n)*HHk;
rowvec = S(k1,kp1:n); Z = max(abs(rowvec));
                            if k < n-2
nm = norm(S(k1,kp2:n),inf);
    if nm < 1e-10
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
S(kp2:n,kp1:n) = HHk*S(kp2:n,kp1:n);
S(k1,kp3:n) = zeros(1,n-kp3+1); %To eliminate roundoff error.
GI(:,kp2:n) = GI(:,kp2:n)*HHk;
CS(kp2:n,:) = HHk*CS(kp2:n,:);
band(k1) = kp2;
    end
                            end %on k < n-2 
band(k1) = kp2;
dm = S(k1,kp1);
           if abs(dm) > 1e-5
tolk = abs(S(k1,kp2)/dm);
			if tolk > tol
%Element S(k1,kp2) cannot be reduced with a gaussian transformation.
            else
%We reduce S(k1,kp2) to zero with a gaussian transformation.
w = S(k1,kp2)/dm;
S(k1+1:n,kp2) = S(k1+1:n,kp2) - w*S(k1+1:n,kp1);
S(k1,kp2) = 0;
GI(:,kp2) = GI(:,kp2) - w*GI(:,kp1);
S(kp1,kp1:n) = S(kp1,kp1:n) + w*S(kp2,kp1:n);
CS(kp1,:) = CS(kp1,:) + w*CS(kp2,:);
band(k1) = kp1; k1 = k1+1;
            end	
           end 
T = sparse(S);           
                    end %on k =1:n-2
k = n-1;                    
sband = band;
disp('1-norm of the reduced matrix S is')
norm(T,1)
figure
spy(S);
hold on
title('Figure 1:  Profile of matrix S')
ylabel('matrix row'); xlabel('nonzero elements')
hold off
disp('number of elements beyond tridiagonal')
xtra = nnz(triu(S,2))
%ADI iteration parameters are determined.
errY = input('Bound on error in Y is ')
    if isempty(errY)
        errY = 1e-6;
    end
FCS = CS*CS';
qgpar
adilow
X = GI*Y*GI';
disp('Estimated ||error in X||/||X||')
errX = norm(C - A*X - X*A',1)/norm(C,1)
return