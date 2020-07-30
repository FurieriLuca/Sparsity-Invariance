function [A_bin] = bin(A)
    % return a binary matrix to represent the pattern
    if ~iscell(A)    % not cell matrix
        if isa(A(1,1),'double')         % scalar matrices
            eps   = 1e-5;           % tolarance
            A_bin = zeros(size(A));
            A_bin(abs(A) > eps) = 1;
        else                        % transfer function matrices
            % to do
        end
    else  %% cell matrix
        A_bin = ones(size(A));
        for i = 1:size(A,1)
            for j = 1:size(A,2)
                if length(A{i,j}) == 1 && A{i,j} == 0
                    A_bin(i,j) = 0;
                end
            end
        end
    end
        
end

