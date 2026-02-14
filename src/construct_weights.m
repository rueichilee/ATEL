function W = construct_weights(X, J, basis_type)
% CONSTRUCT_WEIGHTS Wrapper for user-defined basis functions
%
% Inputs:
%   X          - N x T x P matrix (or vector) of covariates
%   J          - desired basis dimension (2, 3, 4, or 5)
%   basis_type - String: 'bspline', 'polynomial', or 'trigonometric'
%
% Outputs:
%   W          - (N x TJ) array of diversified weights
if ndims(X) == 2
    % X is 2-D (N x T)
    [N, T] = size(X);
    switch lower(basis_type)
        case 'bspline'
            % Call user function: make_bspline_W_J(x, J)
            x_t = X(:);
            W_t = make_bspline_W(x_t, J);  % N-by-J, row i = W_{it}'
            W=[];
            for j =1:J
                W = [W,reshape(W_t(:,j), N, T)];
            end
        case 'polynomial'
            % Call user function: make_trig_W(x, J)
            W = make_poly_W(X, J);
        case 'trigonometric'
            % Call user function: make_poly_W(X, J)
            x_t = X(:);
            W_t = make_trig_W(x_t, J);  % N-by-J, row i = W_{it}'
            W=[];
            for j =1:J
                W = [W,reshape(W_t(:,j), N, T)];
            end            

            
        otherwise
            error('ATEL:UnknownBasis', ...
                  'Basis type must be "bspline", "polynomial", or "trigonometric".');
    end
    
% Check dimensions to prevent hidden bugs
    if size(W, 2) ~= T*J
        error('Basis function output rows (%d) do not match N*TJ (%d).', size(W, 1), T*J);
    end
elseif ndims(X) == 3
    % X is 3-D (N x T x P)
    [N, T, P] = size(X);
    switch lower(basis_type)
        case 'bspline'
            % Call user function: make_bspline_W_J(x, J)
            ws = [];
            for p=1:P
                x = squeeze(X(:,:,p));
                x_t = x(:);
                w_t = make_bspline_W(x_t, J);  % NT-by-J, row i = W_{it}'
                ws=[ws,w_t];
            end
            W=[];
            for j =1:J
                for p=1:P
                    tmp = ws(:,1+(p-1)*J:p*J);
                    W = [W,reshape(tmp(:,j), N, T)];
                end
            end
        case 'polynomial'
            % Call user function: make_trig_W(x, J)
            W = [];
            for p=1:P
                x = squeeze(X(:,:,p));
                x_t = x(:);
                w_t = make_poly_W(x_t, J);
                w=[];
                for j =1:J
                    w = [w,reshape(w_t(:,j), N, T)];
                end
                W = [W,w];
            end
        case 'trigonometric'
            % Call user function: make_poly_W(X, J)
            W = [];
            for p=1:P
                x = squeeze(X(:,:,p));
                x_t = x(:);
                w_t = make_trig_W(x_t, J);  % NT-by-J, row i = W_{it}'
                w=[];
                for j =1:J
                    w = [w,reshape(w_t(:,j), N, T)];
                end
                W = [W,w];
            end
           
        otherwise
            error('ATEL:UnknownBasis', ...
                  'Basis type must be "bspline", "polynomial", or "trigonometric".');
        end
    
% Check dimensions to prevent hidden bugs
    if size(W, 2) ~= T*J*P
        error('Basis function output rows (%d) do not match N*TJP (%d).', size(W, 1), T*J*P);
    end


else
    error('X must be 2-D or 3-D. Got ndims(X) = %d.', ndims(X));
end
end
