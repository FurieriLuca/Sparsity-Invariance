
CYs = sym('CY',[n n*(N+1)]);
CYv = sdpvar(n,n*(N+1));          % decision variables for Y

CUs = sym('CU',[m n*(N+1)]);
CUv = sdpvar(m,n*(N+1));          % decision variables for U


CWs = sym('CW',[n m*(N+1)]);
CWv = sdpvar(n,m*(N+1));          % decision variables for W


CZs = sym('CZ',[m m*(N+1)]);
CZv = sdpvar(m,m*(N+1));          % decision variables for Z


%%Express Y,U,W,Z as FIR transfer matrices of order N
Y = zeros(n,n);
U = zeros(m,n);
W = zeros(n,m);
Z = zeros(m,m);
for(t=1:N+1)
    Y = Y + CYs(:,[(t-1)*n+1:t*n])/z^(t-1);
    U = U + CUs(:,[(t-1)*n+1:t*n])/z^(t-1);
    W= W + CWs(:,[(t-1)*m+1:t*m])/z^(t-1);
    Z= Z + CZs(:,[(t-1)*m+1:t*m])/z^(t-1);
end


%%%achievability constraints
ach1 = Y - G*U-eye(n);
ach2 = W - G*Z;
ach3 = -Y*G+W;
ach4 = -U*G + Z-eye(m);
Constraints=[];


for(i = 1:n)       %ach1
    for(j = 1:n)
        fprintf(' ach1:  Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
        [num,~] = numden(ach1(i,j));
        cc = coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CYs);vec(CUs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CYv);vec(CUv)]== b_eqs];
    end
end

for(i = 1:n)       %ach2
    for(j = 1:m)
        fprintf(' ach2:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
        [num,~] = numden(ach2(i,j));
        cc = coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CWs);vec(CZs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CWv);vec(CZv)]== b_eqs];
    end
end
for(i = 1:n)       %ach3
    for(j = 1:m)
        fprintf(' ach3:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
        [num,~] = numden(ach3(i,j));
        cc = coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CYs);vec(CWs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CYv);vec(CWv)]== b_eqs];
    end
end
for(i = 1:m)       %ach4
    for(j = 1:m)
        fprintf(' ach4:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/m );
        [num,~] = numden(ach4(i,j));
        cc = coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CUs);vec(CZs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CUv);vec(CZv)]== b_eqs];
    end
end
