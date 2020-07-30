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
generate_plant_data;         



%Binary matrices for SI approach
Tbin = Sbin;                                    % Matrix "T"
Rbin = generate_SXlessS(Tbin);      %Matrix "R^*" generated according to SI algorithm
QI   = test_QI(Sbin,Delta);            % variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time


N = 15;                                             %order of finite-dimensional approximation           


fprintf('==================================================\n')
fprintf('              Using the IOP parametrization         \n')
fprintf('              Discrete-time example        \n')
fprintf('==================================================\n')

fprintf('Encoding the achievability constraints ...\n')

achievability_constraints;              % IOP achievability constraints
sparsity_constraints;                    % sparsity constraints

fprintf('Encoding the H2 cost ...\n') %Using the FIR interpretation of H2 norm (frobeinus norm of FIR factors)
cost_matrix = [CWv(:,[(1-1)*m+1:1*m]) CYv(:,[(1-1)*n+1:1*n])-eye(n,n);CZv(:,[(1-1)*m+1:1*m])-eye(m) CUv(:,[(1-1)*n+1:1*n])]'*[CWv(:,[(1-1)*m+1:1*m]) CYv(:,[(1-1)*n+1:1*n])-eye(n,n);CZv(:,[(1-1)*m+1:1*m])-eye(m) CUv(:,[(1-1)*n+1:1*n])];
cost = trace(cost_matrix); 
for(t=2:N+1) 
       cost_matrix = [CWv(:,[(t-1)*m+1:t*m]) CYv(:,[(t-1)*n+1:t*n]);CZv(:,[(t-1)*m+1:t*m]) CUv(:,[(t-1)*n+1:t*n])]'*[CWv(:,[(t-1)*m+1:t*m]) CYv(:,[(t-1)*n+1:t*n]);CZv(:,[(t-1)*m+1:t*m]) CUv(:,[(t-1)*n+1:t*n])];
       cost = cost+trace(cost_matrix);
end


%Convex program
options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);  
sol     = optimize(Constraints,cost,options);
if(sol.problem == 0)
H2norm_squared = value(cost); 
fprintf('\n H2 norm of the closed loop system is %6.4f \n', sqrt(H2norm_squared));
else
    fprintf('\n INFEASIBLE, or numerical problem (check) \n');
end



%Check sparsity of the solution
z = sym('z');
Gz = [0.1/(z-0.5) 0 0 0 0;                             
     0.1/(z-0.5) 1/(z-2) 0 0 0;
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0 0;                                              % plant tf definition
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0.1/(z-0.5) 0;
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0.1/(z-0.5) 1/(z-2)]; 
Yr = zeros(n,n);                                                                                  
Ur = zeros(m,n);
Wr = zeros(n,m);
Zr = zeros(m,m);
for(t=1:N+1)
    Yr = Yr + value(CYv(:,[(t-1)*n+1:t*n]))/z^(t-1);                               % Construct the IOP parameters from the solution                         
    Ur = Ur + value(CUv(:,[(t-1)*n+1:t*n]))/z^(t-1);
    Wr = Wr + value(CWv(:,[(t-1)*m+1:t*m]))/z^(t-1);
    Zr = Zr + value(CZv(:,[(t-1)*m+1:t*m]))/z^(t-1);
end
K = Ur*inv(Yr);                                                                                 % Construct the controller
double(subs(K,z,rand))                                                                       % print its sparsity


