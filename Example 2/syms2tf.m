% Convert Symbolic Transfer Function to ZPK Transfer Function
% Crystal Nassouri 2009
% Allows for substitution/manipulation that can only be done with syms
%
% Ex: Gs = syms2tf(G)
% Where G is a symbolic equation and Gs is a zpk transfer function

function[ans] = syms2tf(G)
n=size(G,1);
m=size(G,2);
for(i=1:n)
    for(j=1:m)
        [symNum,symDen] = numden(G(i,j)); %Get num and den of Symbolic TF
        TFnum = sym2poly(symNum);    %Convert Symbolic num to polynomial
        TFden = sym2poly(symDen);    %Convert Symbolic den to polynomial
        ans(i,j) =tf(TFnum,TFden);
    end
end


