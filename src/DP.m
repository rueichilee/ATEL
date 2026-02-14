function F = DP(Y, W, J, use_parallel)
% DP Calculates the F matrix with optional parallel processing.
%
% Inputs:
%   Y            : Input data matrix (N x T)
%   W            : Weight matrix
%   J            : Number of iterations/blocks
%   use_parallel : Boolean (true/false) to enable parfor
%
% Outputs:
%   F            : Result matrix (T x J)

    % 1. Derive T from the data dimensions
    %    (Assuming Y is N x T, and we want the mean over N)
    T = size(Y, 2);
    
    % 2. Pre-allocate F
    F = zeros(T, J);

    % 3. Execute Loop
    if use_parallel
        % Optional: Automatically start pool if needed
        if isempty(gcp('nocreate'))
            parpool; 
        end
        
        % PARALLEL Execution
        parfor j = 1:J
            % Calculate indices for this slice
            idx_start = 1 + (j-1)*T;
            idx_end = j*T;
            
            % Extract the slice of W to minimize data transfer overhead
            w_slice = W(2:end, idx_start:idx_end);
            
            % Compute
            % Note: Y is broadcast to all workers
            F(:, j) = mean(Y .* w_slice, 1)';
        end
        
    else
        % SERIAL Execution (Standard for loop)
        % Useful for debugging or when J is very small
        for j = 1:J
            idx_start = 1 + (j-1)*T;
            idx_end = j*T;
            
            F(:, j) = mean(Y .* W(2:end, idx_start:idx_end), 1)';
        end
    end
end