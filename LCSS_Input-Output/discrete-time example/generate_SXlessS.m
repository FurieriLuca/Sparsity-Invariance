function [ X ] = generate_SXlessS( S )
% Generate a maximally sparsity-wise invariant (MSI) subplace with respect
% to X
% See Section 3.2 of the following paper
% "Minimizing suboptimality in Distritbued Control: a Framework-Independent Approach"

m = size(S,1);
n = size(S,2);
X = ones(n,n);

% Analytical solution with complexity mn^2
for i = 1:m
    for k = 1:n
        if S(i,k)==0
            for j = 1:n
                if S(i,j) == 1
                    X(j,k) = 0;
                end
            end
        end
    end
end

end

