% Solve an optimal  control problem using the IO parametrization of the paper [1]
%         min <K> ||P11 + P12 * K(I - GK)^(-1) * P21||
%         s.t.      (I-GK)^{-1},K(I-GK)^{-1}, (I-GK)^{-1}G and (I-KG)^-1
%                   are stable
%
% where P11, P12, G, P21 are transfer functions (plant dynamics). In our
% latest note , we propose the following convex representation of the
% problem
%
%         min <X,Y> ||P11 + P12 * Y * P21||
%         s.t.      X, Y, W, Z are stable
%                   [I -G][X W| = [I 0]
%                         |Y Z]
%                         [X W|[-G| = [0|
%                         |Y Z]| I]   |I]
%
% Here we consider a simplified discrete time scenario, where the cost
% admits a FIR representation, allowing for simplified computation of the 
% H2 norm.
%
% Authors: L. Furieri, Automatic Control Laboratory, ETH.
%          Y. Zheng, Department of Engineering Science, University of Oxford
%[1] "An Input-Output Parametrization of Stabilizing Controllers: amidst Youla and System Level Synthesis"


clear all;
clc;
%clc;

N =10;


%% generate plant data
generate_plant_data;         
% QI Sparsity Constraints defined in [1]
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
    0 1 0 0 1];  %not QI
Sbin8 = [0 1 0 0 0;0 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %QI-closest to Abin9
Sbin9 = [1 1 0 0 0;1 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %not QI
Sbin10 = [1 0 0 0 0;1 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %not QI, feasible with sparsity invariance
Sbin11 = [0 0 0 0 0;0 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %QI, closest to Sbin10, not feasible!

Sbin = Sbin6;
QI   = test_QI(Sbin,Delta);        % variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time

Tbin = Sbin;                       % Matrix "T"
Rbin = generate_SXlessS(Tbin);     % Matrix "R_MSI"
%Rbin=eye(n);

%}



fprintf('==================================================\n')
fprintf('              New controller parametrization         \n')
fprintf('==================================================\n')
%% Encode sparsity Y(s) \in Sparse(T) and G(s)Y(s) in Sparse(R)
fprintf('Encoding the achievability constraints ...\n')

achievability_constraints;
sparsity_constraints;

%% Hinf norm minimization, SDP formulation
fprintf('Encoding the H2 cost ...\n')



%%H2cost
cost_matrix=[CWv(:,[(1-1)*m+1:1*m]) CXv(:,[(1-1)*n+1:1*n])-eye(n,n);CZv(:,[(1-1)*m+1:1*m])-eye(m) CYv(:,[(1-1)*n+1:1*n])]'*[CWv(:,[(1-1)*m+1:1*m]) CXv(:,[(1-1)*n+1:1*n])-eye(n,n);CZv(:,[(1-1)*m+1:1*m])-eye(m) CYv(:,[(1-1)*n+1:1*n])];
cost=trace(cost_matrix); %matrix J[0]

for(t=2:N+1) %%build(J[t])
       cost_matrix=[CWv(:,[(t-1)*m+1:t*m]) CXv(:,[(t-1)*n+1:t*n]);CZv(:,[(t-1)*m+1:t*m]) CYv(:,[(t-1)*n+1:t*n])]'*[CWv(:,[(t-1)*m+1:t*m]) CXv(:,[(t-1)*n+1:t*n]);CZv(:,[(t-1)*m+1:t*m]) CYv(:,[(t-1)*n+1:t*n])];
       cost=cost+trace(cost_matrix);
end


fprintf('=====================Yr=============================\n')

options = sdpsettings('allownonconvex',0,'solver','sedumi','verbose',1);
sol     = optimize(Constraints,cost,options);

Vgamma = value(cost);  % value of the H2 norm!
fprintf('\n H2 norm of the closed loop system is %6.4f \n', sqrt(Vgamma));



%CHECKS %%% WARNING, THESE BECOME VERY SLOW FOR N>5
z=sym('z');
Gz = [0.1/(z-0.5) 0 0 0 0;                             
     0.1/(z-0.5) 1/(z-2) 0 0 0;
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0 0;                                              % plant tf definition
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0.1/(z-0.5) 0;
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0.1/(z-0.5) 1/(z-2)]; 
Xr=zeros(n,n);
Yr=zeros(m,n);
Wr=zeros(n,m);
Zr=zeros(m,m);
for(t=1:N+1)
    Xr=Xr+value(CXv(:,[(t-1)*n+1:t*n]))/z^(t-1);
    Yr=Yr+value(CYv(:,[(t-1)*n+1:t*n]))/z^(t-1);
    Wr=Wr+value(CWv(:,[(t-1)*m+1:t*m]))/z^(t-1);
    Zr=Zr+value(CZv(:,[(t-1)*m+1:t*m]))/z^(t-1);
end
K=Yr*inv(Xr);

%checks;


