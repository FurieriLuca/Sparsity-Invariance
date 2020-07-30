%%%     Code associated with the paper 
%%%     "Sparsity Invariance for Convex Design of Distributed Controllers"
%%%
%%%     Authors: Luca Furieri, Yang Zheng, Antonis Papachristodoulou, Maryam
%%%     Kamgarpour
%%%
%%%     All rights reserved

%%%     To be cited as follows:
%{

@article{furieri2020sparsity,
  title={Sparsity invariance for convex design of distributed controllers},
  author={Furieri, Luca and Zheng, Yang and Papachristodoulou, Antonis and Kamgarpour, Maryam},
  journal={IEEE Transactions on Control of Network Systems},
  year={2020},
  publisher={IEEE}
}

%}

clear all;
clc;





% Sparsity patterns for simulation
S = [1 0 0 0 0;1 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1];     %not QI, feasible with SI
S2 = [0 0 0 0 0;0 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %QI, closest to S, not feasible! (closest QI subset approach)
S3= [1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];     %QI, closest superset

%%  SELECT WHICH SPARSITY TO BE SIMULATED (Write either "S", "S2" or "S3")
Sbin = S;



%% generate plant data
N = 8;           % order of the controller
a = 2;           % defines the basis for RH_infinity as {1/(s+a)^i}
generate_plant_data;         



%Binary matrices for SI approach
Tbin = Sbin;                                    % Matrix "T"
Rbin = generate_SXlessS(Tbin);      %Matrix "R^*" generated according to SI algorithm
QI   = test_QI(Sbin,Delta);            % variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time





fprintf('==================================================\n')
fprintf('              Using the Youla parametrization         \n')
fprintf('              Continuous-time example        \n')
fprintf('==================================================\n')

%% Encode sparsity Y(s) \in Sparse(T) and G(s)Y(s) in Sparse(R)
sparsity_constraints;

%% H2 norm minimization, SDP formulation
fprintf('Step 3: Encoding the other LMI constraint ...')
Qstatic = [CQv DQv];
P       = sdpvar(size(A1_hat,1),size(A1_hat,1));
S       = sdpvar(size(A1_hat,2),size(B2_hat,1));
R       = sdpvar(size(A2_hat,1),size(A2_hat,1));
L       = sdpvar(size(C2_hat,1),size(C2_hat,1));
gamma   = sdpvar(1,1);

Constraints = [Constraints,trace(L)<=gamma, P>=0,R>=0];                                                                   % (27)-(28) in our paper, Section V
Constraints = [Constraints,
    [A1_hat*P+P*A1_hat' A1_hat*S-S*A2_hat+A_hat+B_hat*Qstatic*C_hat B1_hat+B_hat*Qstatic*F_hat-S*B2_hat;
    (A1_hat*S-S*A2_hat+A_hat+B_hat*Qstatic*C_hat)' R*A2_hat+A2_hat'*R R*B2_hat;
    (B1_hat+B_hat*Qstatic*F_hat-S*B2_hat)' B2_hat'*R -eye(size(B2_hat,2))] <= 0];                                         % (29)

Constraints = [Constraints,
    [P zeros(size(P,1),size(R,2)) P*C1_hat';zeros(size(R,1),size(P,2)) R (C2_hat+E_hat*Qstatic*C_hat+C1_hat*S)';
    C1_hat*P C2_hat+E_hat*Qstatic*C_hat+C1_hat*S L] >= 0,DQv == 0];                                                       % (30), DQ=0 to guarantee that \mathcal{D}=0.
fprintf('Done \n')

% options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
fprintf('Step 4: call SDP solver to obtain a solution ... \n')
fprintf('==================================================\n')

options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
sol     = optimize(Constraints,gamma,options);

%CQ   = round(value(CQv),6);  %rounding to avoid false non-zeros
%DQ   = round(value(DQv),6);
CQ  = value(CQv);
DQ  = value(DQv);

Vgamma = sqrt(value(gamma));  % value of the H2 norm!
fprintf('\n H2 norm of the closed loop system is %6.4f \n', Vgamma);

%% RECOVER Q(s) and K(s) to check sparsities
Gi = (s*eye(size(AiQ,1))-AiQ)\BiQ;
for i = 1:n
    for j = 1:m
        Y(i,j) = CQ(i,(j-1)*N+1:j*N)*Gi + DQ(i,j);
    end
end
%K = Y/(eye(n)+Gs*Y);
K = Knom+Y*inv(eye(n)+Gnom*Y);
GQ      = G*Y;
GQsubs  = double(subs(GQ,s,rand));           
GQbin   = bin(GQsubs);
Ksubs   = double(subs(K,s,rand));           
Kbin    = bin(Ksubs);



