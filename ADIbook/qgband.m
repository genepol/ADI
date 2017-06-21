%QGBAND computes a banded upper Hessenberg matrix S similar to %snorm*A.
%The last nonzero element in row k is element k of sband.
%Columns are reduced successively with Householder transformations.
%After column k is reduced, row k1 (the last unreduced row beyond column 
%k+1is reduced with a householder transformation beyond column k+2.  
%Element S(k1,k+2) then has magnitude 
%equal to the 2-norm of S(k1,k+2:n) prior to reduction. If 
%|S(k1,k+2)/S(k1,k+1)| < tol, S(k1,k+2) is reduced to zero with a gaussian
%transformation and k1 is increased to k1+1. 10/31/11
A = input('Matrix A is:');
if isempty(A)
    A = randn(30);
end
myzero = 1e-10;
n=length(A); nA = n; snorm = 1/sqrt(norm(A,1)*norm(A,inf)); 
k1 = 1; S = snorm*A;
%This normalizes S so its eigenvalue magnitudes are bounded by unity.
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
%MATLAB sgn(0) = 0 and we need sgn(0)=1 here so:
    if vk(1) == 0
        sgn = 1;
    else
        sgn = sign(vk(1));
    end
u1 = sgn*norm(vk,2);
uk = [u1;zeros(n-kp1,1)]; wk = uk + vk; nmk = wk'*wk;
HHk = eye(n-k) - 2*(wk*wk')/nmk;
S(kp1:n,k:n) = HHk*S(kp1:n,k:n);
S(k1:n,kp1:n) = S(k1:n,kp1:n)*HHk;
S(kp2:n,k) = zeros(n-kp1,1); %To eliminate roundoff error.
%The column has been reduced.
%The matrix is banded above k.
rowvec = S(k1,kp2:n); Z = max(abs(rowvec));
                            if k < n-2
    if Z < myzero
%The row does not have to be reduced.        
        S(k1,kp2:n) = zeros(1,n-kp1);
        band(k1) = kp1; 
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
S(kp2:n,kp1:n) = HHk*S(kp2:n,kp1:n);
S(k1,kp3:n) = zeros(1,n-kp3+1); %To eliminate roundoff error.
band(k1) = kp2;
    end
                            end %on k < n-2 
      if Z > myzero            
           dm = S(k1,kp1);
           if abs(dm) > sqrt(myzero)
%We need a significant dm to attempt row reduction.               
tolk = abs(S(k1,kp2)/dm);
                if tolk < tol
%We reduce S(k1,kp2) to zero with a gaussian transformation.
%Else element S(k1,kp2) cannot be reduced with a gaussian transformation.
w = S(k1,kp2)/dm;
S(k1+1:n,kp2) = S(k1+1:n,kp2) - w*S(k1+1:n,kp1);
S(k1,kp2) = 0;
S(kp1,kp1:n) = S(kp1,kp1:n) + w*S(kp2,kp1:n);
band(k1) = kp1; k1 = k1+1;
                end
           end
      end %on Z > myzero
        end %on k =1:n-2        
k = n-1;                
if abs(S(n-2,n))< myzero
    S(n-2,n) = 0;
    band(n-2) = n-1;
end
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
eigon = input('To compute eigenvalues enter 1, else enter 0')
    if eigon == 1
%nA-2 eigenvalues of S are found with MATLAB sparse matrix routine eigs.
%My MATLAB eigs does not permit all nA values with one command:
T = sparse(S);
opts.tol = 1e-6;
opts.disp = 0;
lamS = eigs(T,nA-2,0,opts);
eigsofA = lamS/snorm;
disp('eigenvalues of A are in eigsofA')
    end
return