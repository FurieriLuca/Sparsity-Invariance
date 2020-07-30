
%CHECKS

CQYr  = value(CQYv);
DQYr  = value(DQYv);
for(i=1:m)
        for(j=1:n)
        Yr(i,j)=(CQYr(i,[(j-1)*N+1:j*N])*inv(s*eye(N)-AiQ)*BiQ+DQYr(i,j));
        end
end
eq=eye(n)+Gs*Yr;

K=Yr*inv(eq);
Ksubs=subs(K,s,rand);
for(i=1:m)
    for(j=1:n)
        if(abs(Ksubs(i,j))<0.0001)
            K(i,j)=0;
        end
    end
end
K=simplifyFraction(K);
[num,den]=numden(K);
for(i=1:20)
    powers_vec(i)=s^(i-1);
end
assume(s,'real');
for(i=1:m)
    for(j=1:n)
        if(K(i,j)~=0)
            numc=coeffs(num(i,j),s,'All');
            denc=coeffs(den(i,j),s,'All');
            numc=round(double(numc/denc(1)),5);
            denc=round(double(denc/denc(1)),5);
            K(i,j)=flip(numc)*powers_vec(1:size(numc,2))'/(flip(denc)*powers_vec(1:size(denc,2))');
        end
    end
end
K=simplifyFraction(K)


%Check stability parameters.  A BIT CHALLENGING
tollerance=1e-3;

disp('Check (I-GK)^-1')
eq1 = inv(eye(n)-Gs*K); %%(I-G*K)^-1=(I-G*Y*X^-1)^-1=(I-G*Y*(I+G*Y)^-1)^-1
eq1subbed=subs(eq1,s,rand);
for i = 1:n
        for j = 1:n
                fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n )
                
                if(isstable(syms2tf(eq1(i,j)))==0 && abs(eq1subbed(i,j))>tollerance)
                        fprintf('   Apparently unstable: check close zero/pole cancellations\n')
                        [num,den]=numden(eq1(i,j));
                        poles=double(roots(coeffs(den,s,'All')));
                        zeroes=double(roots(coeffs(num,s,'All')));
                        
                        unstable_poles=[];
                        for(k=1:size(poles))
                                if(real(poles(k))>=0)
                                        unstable_poles=[unstable_poles, poles(k)];
                                end
                        end
                        stable=zeros(1,size(unstable_poles,2));
                        zeroes_comparison=zeroes;
                        for(l=1:size(unstable_poles,2))
                                done=0;
                                for(v=1:size(zeroes,1))
                                        if(abs(zeroes_comparison(v)-unstable_poles(l))<tollerance &&done==0) %%still imprecise... in general, I should also remove the cancelling zero/pole and start again. In this example this does not matter
                                                stable(l)=1;
                                                zeroes_comparison(v)=inf; %not to use it again for the next unstable poles, this zero is already used.
                                                done=1;
                                        end
                                end
                        end
                        for(z=1:size(stable,2))
                                if(stable(z)==0)
                                        disp('******ACTUALLY UNSTABLE :-( ******')
                                end
                        end
                        
                end
        end
end


eq2 = K*inv(eye(5)-Gs*K);
eq2subbed=subs(eq2,s,rand);

disp('Check K(I-GK)^-1')

for i = 1:m
        for j = 1:n
                fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n )
                
                if(isstable(syms2tf(eq2(i,j)))==0 && abs(eq2subbed(i,j))>tollerance)
                        fprintf('   Apparently unstable: check close zero/pole cancellations\n')
                        [num,den]=numden(eq2(i,j));
                        poles=double(roots(coeffs(den,s,'All')));
                        zeroes=double(roots(coeffs(num,s,'All')));
                        
                        unstable_poles=[];
                        for(k=1:size(poles))
                                if(real(poles(k))>=0)
                                        unstable_poles=[unstable_poles, poles(k)];
                                end
                        end
                        stable=zeros(1,size(unstable_poles,2));
                        zeroes_comparison=zeroes;
                        for(l=1:size(unstable_poles,2))
                                done=0;
                                for(v=1:size(zeroes,1))
                                        if(abs(zeroes_comparison(v)-unstable_poles(l))<tollerance && done==0) %%still imprecise... in general, I should also remove the cancelling zero/pole and start again. In this example this does not matter
                                                stable(l)=1;
                                                zeroes_comparison(v)=inf; %not to use it again for the next unstable poles, this zero is already used.
                                                done=1;
                                        end
                                end
                        end
                        for(z=1:size(stable,2))
                                if(stable(z)==0)
                                        disp('******ACTUALLY UNSTABLE :-( ******')
                                end
                        end
                        
                end
        end
end


eq3 = inv(eye(n)-Gs*K)*Gs;
eq3subbed = subs(eq3,s,rand);
disp('Check (I-GK)^-1G')

for i = 1:n
        for j = 1:m
                fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n )
                
                if(isstable(syms2tf(eq3(i,j)))==0 && abs(eq3subbed(i,j))>tollerance)
                        fprintf('   Apparently unstable: check close zero/pole cancellations\n')
                        [num,den]=numden(eq3(i,j));
                        poles=double(roots(coeffs(den,s,'All')));
                        zeroes=double(roots(coeffs(num,s,'All')));
                        
                        unstable_poles=[];
                        for(k=1:size(poles))
                                if(real(poles(k))>=0)
                                        unstable_poles=[unstable_poles, poles(k)];  %%THIS IS SOMEWHAT FLAWED, does not recognize the case double 1.
                                end
                        end
                        stable=zeros(1,size(unstable_poles,2));
                        zeroes_comparison=zeroes;
                        for(l=1:size(unstable_poles,2))
                                done=0;
                                for(v=1:size(zeroes,1))
                                        if(abs(zeroes_comparison(v)-unstable_poles(l))<tollerance && done==0) %%still imprecise... in general, I should also remove the cancelling zero/pole and start again. In this example this does not matter
                                                stable(l)=1;
                                                zeroes_comparison(v)=inf; %not to use it again for the next unstable poles, this zero is already used.
                                                done=1;
                                        end
                                end
                        end
                        for(z=1:size(stable,2))
                                if(stable(z)==0)
                                        i
                                        j
                                        disp('******ACTUALLY UNSTABLE :-( ******')
                                end
                        end
                        
                end
        end
end


eq4 = inv(eye(n)-K*Gs);
eq4subbed=subs(eq4,s,rand);
disp('Check (I-KG)^-1')

for i = 1:m
        for j = 1:m
                fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n )
                if(isstable(syms2tf(eq4(i,j)))==0 && abs(eq4subbed(i,j))>tollerance)
                        fprintf('   Apparently unstable: check close zero/pole cancellations\n')
                        [num,den]=numden(eq4(i,j));
                        poles=double(roots(coeffs(den,s,'All')));
                        zeroes=double(roots(coeffs(num,s,'All')));
                        
                        unstable_poles=[];
                        for(k=1:size(poles))
                                if(real(poles(k))>=0)
                                        unstable_poles=[unstable_poles, poles(k)];
                                end
                        end
                        stable=zeros(1,size(unstable_poles,2));
                        zeroes_comparison=zeroes;
                        for(l=1:size(unstable_poles,2))
                                done=0;
                                for(v=1:size(zeroes,1))
                                        if(abs(zeroes_comparison(v)-unstable_poles(l))<tollerance && done==0) %%still imprecise... in general, I should also remove the cancelling zero/pole and start again. In this example this does not matter
                                                stable(l)=1;
                                                zeroes_comparison(v)=inf; %not to use it again for the next unstable poles, this zero is already used.
                                                done=1;
                                        end
                                end
                        end
                        for(z=1:size(stable,2))
                                if(stable(z)==0)
                                        disp('******ACTUALLY UNSTABLE :-( ******')
                                end
                        end
                        
                end
        end
end


