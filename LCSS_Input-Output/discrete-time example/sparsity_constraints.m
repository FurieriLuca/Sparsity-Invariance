fprintf('Step 1: Encoding the constraint Y(s) in T ...\n')

for i = 1:m
    for j = 1:n
        if Tbin(i,j) == 0 % main cycle
            for(k=1:N+1)
                Constraints = [Constraints, CYv(i,j+(k-1)*n) == 0];
                
            end
        end
    end
end

if(QI==0)
    for i = 1:n
        for j = 1:n
            if Rbin(i,j) == 0 % main cycle
                for(k=1:N+1)
                    Constraints = [Constraints, CXv(i,j+(k-1)*n) == 0];
                end
            end
        end
    end
end
fprintf('Done\n')


%{
[num,~]=numden(Y);
for(i=1:m)       %sparsity
    for(j=1:n)
        if(Tbin(i,j)==0)
            cc=coeffs(num(i,j),z,'All');
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CYs)]);
            A_eqs   = double(A_eq);
            b_eqs    = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(CYv)]== b_eqs];
        end
    end
%}
