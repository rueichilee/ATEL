function W = make_poly_W(X, J)
%MAKE_POLY_W Construct polynomial basis features from X.
%
%   W = MAKE_POLY_W(X, J) builds W = [X.^1, X.^2, ..., X.^J].
%
% Inputs:
%   X    : N-by-T matrix (or vector) of covariates
%   J    : positive integer, polynomial degree / basis dimension
%
% Output:
%   W : N-by-(J*T) matrix if X is N-by-T
%       N-by-J      matrix if X is N-by-1 (vector)
%
% Notes:
%   - Centering/scaling is strongly recommended for stability.
%   - If you need W as N-by-J-by-T (time-varying basis per period),
%     you can reshape afterwards, or I can provide a 3D version.


    if ~(isscalar(J) && J >= 1 && floor(J) == J)
        error('J must be a positive integer.');
    end

    % Ensure numeric
    if ~isnumeric(X)
        error('X must be numeric.');
    end

    % Convert vector to N-by-1
    if isvector(X)
        X = X(:);
    end

    Xw = X;

    [N, T] = size(Xw);

    % Build polynomial blocks: [X.^1, X.^2, ..., X.^J]
    W = zeros(N, J*T);
    for j = 1:J
        cols = (1:T) + (j-1)*T;
        W(:, cols) = Xw.^j;
    end
end