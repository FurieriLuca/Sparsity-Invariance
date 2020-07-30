
   
Gnom = G*inv(eye(n)-Knom*G);

I = eye(n);
% Defines CQ and DQ both as symbolic variables and sdpvar variables
CQs = sym('CQ',[m n*N]);                            
DQs = sym('DQ',[m n]);
CQv = sdpvar(m,n*N);                                                                        % decision variables
DQv = sdpvar(m,n);

Constraints = [];

%% Constraint Type 1: Y(s) \in Sparse(T) in terms of CQv and DQv (same as [2])
fprintf('Step 1: Encoding the constraint Y(s) in T ...')
for i = 1:m                                                                                  
    for j = 1:n
        if Tbin(i,j) == 0 % main cycle
            Constraints = [Constraints, CQv(i,[(j-1)*N+1:j*N]) == 0, DQv(i,j) == 0];
        end
    end
end
fprintf('Done \n')

%% Constraint Type 2: G(s)Y(s) \in Sparse(R) in terms of CQv and DQv 
fprintf('Step 2: Encoding the constraint GY(s) in R ...\n')
if QI == 0        %This cycle is useless if QI (redundant constraints). Hence, we skip it in this case.
     Gi = (s*eye(N)-AiQ)\BiQ;
    for i = 1:n                                                                                                                
        for j = 1:n
            fprintf('   Encoding constraints progress %6.4f \n', 100*(n*(i-1)+j)/n/n );                                                        
            if Rbin(i,j) == 0                                                                       % Whenever we need GY(i,j) = 0 ....
                Uij = Gnom(i,:)* (CQs(:,(j-1)*N+1:j*N) * Gi + DQs(:,j));        % Performs the symbolic matrix product GY(i,j)=\sum_l G(i,l)Y(l,j)
                [num,~] = numden(Uij);
                cc      = coeffs(num,s);                                                       % All elements of this vector must be 0....
                A_eq    = equationsToMatrix(cc,[vec(CQs);vec(DQs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
                A_eqs   = double(A_eq);                                                     %A_eqs is the same as A_eq, for computation with sdpvars
                Constraints = [Constraints, A_eqs*[vec(CQv);vec(DQv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
            end
        end
    end
end
fprintf('Encoding the constraint GY(s) in R ...Done\n')