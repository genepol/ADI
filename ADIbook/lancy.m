%LANCY is the Lanczos algorithm with matrix G
%approximated by low rank H = VT*T*VT'. 10/28/12.
G = input('matrix G is:');
    if isempty(G)
        G1 = 10*rand(30); G = G1 + G1';
    end
if norm(G - G',1) > 1e-8
    error('G not symmetric')
end
n = length(G);
gg = norm(G,1);ig = n/gg;
G = ig*G;
Gin = G;
    lowbnd = input('low bound on elements of T:')
        if isempty(lowbnd)
            lowbnd = .001;
        end
mm = input('maximum rank of W');
    if isempty(mm)
        mm = ceil(sqrt(n));
    end
Vk = zeros(n,mm);
Tk = zeros(mm);
cold = 0;
korder = 1;
H = zeros(n);    
delEG = 1.0;
errG = input('Desired accuracy is:')
    if isempty(errG)
        errG = .001
    end
k = 1;
% k is the number of low rank matrices.
krank = 0;
%krank is the rank of matrix k.
EG = G;
V1 = ones(n,1)/sqrt(n);
%V1 is the initial Lanczos vector.
                    while delEG > errG
V = zeros(n,mm);
V(:,1) = V1;
j = 1;
bett = 1;
bet(1) = 1;
G = EG;
while bett > lowbnd & j <= mm
        wv = G*V(:,j);
            if j ~=1 
                wv = wv - bet(j)*V(:,j-1);
            end
        alp(j) = wv'*V(:,j);
        wv = wv - alp(j)*V(:,j);
        bet(j+1) = norm(wv);
        V(:,j+1) = wv/bet(j+1);
        bett = bet(j+1)/bet(1);
        j = j+1;
    end
    r = j - 1;
    V1 = V(:,j);    
    T = zeros(r);
    T(1,1) = alp(1);    
        if r > 1
            T(1,2) = bet(2);
                for j = 2:r-1
        T(j,j) = alp(j);
        T(j-1,j) = bet(j);    
        T(j,j-1) = bet(j);
                end
            T(r-1,r) = bet(r);
            T(r,r-1) = bet(r);
            T(r,r) = alp(r);
        end
VT = V(:,1:r);
H = H + VT*T*VT';
EG = Gin - H;
korder(k) = r;
Vk(:,cold+1:cold+r) = VT;
[rold cold] = size(Vk);
Tk(1:r,1:r,k) = T;
k = k+1;
delEG = norm(EG,1)/n;
                    end
gi = 1/ig;
G = gi*Gin;
Tk = gi*Tk;
H = gi*H;
disp('norm(G - H,1)/norm(G,1) =')
delEG
return