equationX=eye(n)+Gs*Y; %will force this to be stable
[NUM,DEN]=numden(equationX);
maxdegnum=0;
for i=1:n
    for j=1:n
        if(size(coeffs(NUM(i,j),s),2)>maxdegnum)
            maxdegnum=size(coeffs(NUM(i,j),s),2);
        end
    end
end

freeXs = sym('freeX',[n n*maxdegnum]);
freeXv = sdpvar(n,n*maxdegnum);   
for(i=1:maxdegnum)
    powers_s(i,1)=s^(i-1);
end



for i=1:n
    for j=1:n
        freepoly(i,j)=freeXs(i,[(j-1)*maxdegnum+1:j*maxdegnum])*powers_s;
    end
end


for i = 1:n
    for j = 1:n        
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
        poles=roots(coeffs(DEN(i,j),s,'All'));
        U=1;
        for(k=1:size(poles))
            if(real(poles(k))>0)
                U=U*(s-poles(k));
            end
        end
        polyno=NUM(i,j)-U*freepoly(i,j);
        cc=coeffs(polyno,s);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQYs);vec(DQYs);vec(freeXs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
        A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CQYv);vec(DQYv);vec(freeXv)]== b_eqs];        
    end
end
fprintf('Done \n')
clear freepoly
clear powers_s