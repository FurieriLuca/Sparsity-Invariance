
CXs = sym('CX',[n n*(N+1)]);
CXv = sdpvar(n,n*(N+1));          % decision variables for X

CYs = sym('CY',[m n*(N+1)]);
CYv = sdpvar(m,n*(N+1));          % decision variables for Y


CWs = sym('CW',[n m*(N+1)]);
CWv = sdpvar(n,m*(N+1));          % decision variables for W


CZs = sym('CZ',[m m*(N+1)]);
CZv = sdpvar(m,m*(N+1));          % decision variables for Z


%%Express X,Y,W,Z as FIR transfer matrices of order N
X=zeros(n,n);
Y=zeros(m,n);
W=zeros(n,m);
Z=zeros(m,m);
for(t=1:N+1)
    X=X+CXs(:,[(t-1)*n+1:t*n])/z^(t-1);
    Y=Y+CYs(:,[(t-1)*n+1:t*n])/z^(t-1);
    W=W+CWs(:,[(t-1)*m+1:t*m])/z^(t-1);
    Z=Z+CZs(:,[(t-1)*m+1:t*m])/z^(t-1);
end


%%%achievability constraints
ach1 = X-G*Y-eye(n);
ach2 = W-G*Z;
ach3 = -X*G+W;
ach4 = -Y*G+Z-eye(m);
Constraints=[];


for(i=1:n)       %ach1
    for(j=1:n)
        fprintf(' ach1:  Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
        [num,~]=numden(ach1(i,j));
        cc=coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CXs);vec(CYs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CXv);vec(CYv)]== b_eqs];
    end
end

for(i=1:n)       %ach2
    for(j=1:m)
        fprintf(' ach2:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
        [num,~]=numden(ach2(i,j));
        cc=coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CWs);vec(CZs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CWv);vec(CZv)]== b_eqs];
    end
end
for(i=1:n)       %ach3
    for(j=1:m)
        fprintf(' ach3:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
        [num,~]=numden(ach3(i,j));
        cc=coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CXs);vec(CWs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CXv);vec(CWv)]== b_eqs];
    end
end
for(i=1:m)       %ach4
    for(j=1:m)
        fprintf(' ach4:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/m );
        [num,~]=numden(ach4(i,j));
        cc=coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CYs);vec(CZs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CYv);vec(CZv)]== b_eqs];
    end
end
