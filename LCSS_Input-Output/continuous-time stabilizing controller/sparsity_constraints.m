



%%These are needed to build the transfer function from the variables
shift = 1;
AiQ   = -a*eye(N) + diag(ones(N-abs(shift),1),shift);        
AQ    = kron(eye(n),AiQ); %from size of P21...2n
I     = eye(N);
BiQ   = I(:,N);
BQ    = kron(eye(n),BiQ);


% Defines CQ and DQ both as symbolic variables and sdpvar variables, for
% four paramenters X=I+GY,Y,W=(I+GY)G,Z=I+YG

CQYs = sym('CQY',[m n*N]);
DQYs = sym('DQY',[m n]);
CQYv = sdpvar(m,n*N);          % decision variables
DQYv = sdpvar(m,n);



%Force Y in RHinfinity

for(i=1:m)
    for(j=1:n)
        Y(i,j)=(CQYs(i,[(j-1)*N+1:j*N])*inv(s*eye(N)-AiQ)*BiQ+DQYs(i,j));
    end
end



Constraints = [];

%% Constraint "Achievability": X=I+GY,W=XG,Z=I+YG



fprintf('Step 1-a: Encoding the constraint X=I+GY is stable \n')
constraints_X;
 
fprintf('Step 1-b: Encoding the constraint W=(I+GY)G is stable \n')
constraints_W;

fprintf('Step 1-c: Encoding the constraint Z=I+YG is stable \n')
Constraints_Z;


%%% Constraints "Sparsity"
sparsity_Y

sparsity_X



