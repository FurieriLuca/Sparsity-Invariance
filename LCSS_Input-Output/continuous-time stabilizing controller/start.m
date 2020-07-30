% Solves an LP to find a stabilizing continuous time distributed controller
% with our new parametrizations

clear all;
clc;
%clc;

N = 2;           % order of the controller
a = 3;           % defines the basis for RH_infinity as {1/(s+a)^i}

%% generate plant data

s = sym('s');
Gs = [1/(s+1) 0 0 0 0;
    1/(s+1) 1/(s-1) 0 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 0 0;                                              % plant tf definition
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 1/(s-1)];
Delta = [1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];                   % plant binary structure

n=size(Gs,1);
m=size(Gs,2);

%% QI Sparsity Constraints defined in [1]
Sbin1 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 1];
Sbin2 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 1];
Sbin3 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 0;1 1 0 0 1];
Sbin4 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 0;1 1 1 0 1];
Sbin5 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 1 0 0;1 1 1 0 1];
Sbin6 = [1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];
centralized = ones(m,n);

%% non-QI Sparsity Constraint we study
Sbin7 = [1 0 0 0 0;
    1 1 0 0 0;
    0 1 0 0 0;
    0 1 0 0 0;
    0 1 0 0 1];
Sbin8 = [0 1 0 0 0;0 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %QI-closest to Abin9

Sbin = Delta;
QI   = test_QI(Sbin,Delta);        % variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time
%QI = 0; %%%%%REMOVE


Tbin = Sbin;                       % Matrix "T"
Rbin = generate_SXlessS(Tbin);     % Matrix "R_MSI"
%Rbin=eye(n);

fprintf('==================================================\n')
fprintf('          Compute Stabilizing Distributed Controller         \n')
fprintf('==================================================\n')

%% Encode sparsity Y(s) \in Sparse(T) and G(s)Y(s) in Sparse(R)
sparsity_constraints;


% options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
fprintf('Step 5: call SDP solver to obtain a solution ... \n')
fprintf('=====================Yr=============================\n')

options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
sol     = optimize(Constraints,abs(100*CQYv(2,4)+200*CQYv(4,4)),options); %abs(100*CQYv(2,4)+200*CQYv(4,4))

K0 = K


disp('\n\n\n We will now perform stability checks\n')
checks;
