function W = make_trig_W(x, J)
% MAKE_TRIG_W  Trigonometric polynomial basis (scalar x).
%   x : N-by-1 vector
%   J : desired dimension of W (J >= 1)
%       Basis: 1, cos(2*pi*z), sin(2*pi*z), ..., up to K frequencies.
%
%   W : N-by-J matrix, row i = W(x_i)'

    x = x(:);                        % ensure column
    N = numel(x);

    if J < 1
        error('J must be >= 1.');
    end

    % rescale to [0,1]
    xmin = min(x);
    xmax = max(x);
    if xmax == xmin
        % all x identical: set z = 0.5 (arbitrary)
        z = 0.5 * ones(N,1);
    else
        z = (x - xmin) / (xmax - xmin);
    end

    % number of frequencies
    K = floor((J - 1) / 2);

    % allocate
    W = zeros(N, J);
    col = 1;

    % intercept
    W(:, col) = 1;
    col = col + 1;

    % add cos/sin pairs
    for k = 1:K
        if col <= J
            W(:, col) = cos(2*pi*k*z);
            col = col + 1;
        end
        if col <= J
            W(:, col) = sin(2*pi*k*z);
            col = col + 1;
        end
    end
end