%PLANT DEFINITION
% Example taken from [2]
z = sym('z');
G = [0.1/(z-0.5) 0 0 0 0;                             
     0.1/(z-0.5) 1/(z-2) 0 0 0;
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0 0;                                              % plant tf definition
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0.1/(z-0.5) 0;
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0.1/(z-0.5) 1/(z-2)];                               
Delta = [1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];                   % plant binary structure
n     = size(G,1);
m     = size(G,2);
