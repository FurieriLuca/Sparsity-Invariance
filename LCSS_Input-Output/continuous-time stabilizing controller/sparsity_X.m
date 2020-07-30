fprintf('Step 2-b: Encoding the constraint X in Sparse(R) ... only if not QI')

[NUM,DEN]=numden(equationX);

if QI == 0        %This cycle is useless if QI (redundant constraints). Hence, we skip it in this case.
    for i = 1:n
        for j = 1:n
            fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
            if Rbin(i,j) == 0     % Whenever we need GY(i,j) = 0 ....
                cc      = coeffs(NUM(i,j),s,'All');                                  % All elements of this vector must be 0....
                [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQYs);vec(DQYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
                A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars      
                b_eqs =double(b_eq);
                Constraints = [Constraints, A_eqs*[vec(CQYv);vec(DQYv)] == b_eqs]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
            end
        end
    end
end