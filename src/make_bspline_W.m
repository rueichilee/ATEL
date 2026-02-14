function W = make_bspline_W_J(x, J)
% MAKE_BSPLINE_W_J Build a B-spline basis W from x with specified dim J.
%   x  : NJ-by-1 vector of covariates
%   J  : desired basis dimension (2, 3, 4, or 5)
%   W  : N-by-J matrix, row i = W_{it}' (B-spline basis at x(i))
%   breaks, k : the underlying breakpoints and spline order used
%
% This function:
%   - chooses order k and number of breaks from J,
%   - builds a knot sequence with augknt,
%   - uses spcol on unique sorted x to get *function values* (no derivatives),
%   - maps back to original x order, including duplicates.

    x = x(:);  % ensure column
    if J <= 1
        error('J must be larger than 1');
    end

    % ----- choose spline order k and number of breaks -----
    if J == 2
        % Linear splines: order 2, no interior knots
        k        = 2;   % degree 1
        nBreaks  = 2;   % [min, max]
    elseif J == 3
        % Quadratic splines: order 3, no interior knots
        k        = 3;   % degree 2
        nBreaks  = 2;   % [min, max]
    else
        % J > 4 : use cubic splines
        k        = 4;           % degree 3
        nBreaks  = J - 2;       % length(breaks) = J - 2
    end

    % ----- build breakpoints over the range of x -----
    xmin   = min(x);
    xmax   = max(x);
    if nBreaks == 2
        breaks = [xmin, xmax];
    else
        breaks = linspace(xmin, xmax, nBreaks);
    end

    % ----- full knot sequence for order-k splines -----
    knots = augknt(breaks, k);

    % ----- sorted unique x to avoid derivative behavior in spcol -----
    [x_sorted, idx_sort] = sort(x);
    [x_unique, ~, idx_unique] = unique(x_sorted, 'stable');

    % spcol(knots, k, tau) with strictly increasing tau -> 0th derivative (values)
    W_unique = spcol(knots, k, x_unique);

    % Expand back to all sorted points (duplicate rows for duplicate x's)
    W_sorted = W_unique(idx_unique, :);

    % Undo the sort so rows match the original x
    W = zeros(size(W_sorted));
    W(idx_sort, :) = W_sorted;
end