function [ G_simplified ] = simplify_tf( G )
m=size(G,1);
n=size(G,2);
[num,den]=numden(G);
s=sym('s');

for(i=1:m)
    for(j=1:n)
        i;
        j;

        num_roots=round(double(roots(coeffs(num(i,j),s,'All'))),6);
        den_roots=round(double(roots(coeffs(den(i,j),s,'All'))),6);
        num_new=1;
        den_new=1;
        for(k=1:size(num_roots))
            num_new=num_new*(s-num_roots(k));
        end
        for(k=1:size(den_roots))
            den_new=den_new*(s-den_roots(k));
        end
            G_simplified(i,j)=simplifyFraction(num_new/den_new*s/s);
            
    end
end

end

