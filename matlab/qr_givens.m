function [Q, R] = qr_givens(A)

[m, n] = size(A);
% pre-allocation for performance
G = eye(m, m);

% iteration variables
Q = eye(m, m);
R = A;

% loop over columns of R from left to right
for icol = 1:n

    % make elements below the diagonal zero
    % from bottom row to upper row
    for irow = m:-1:icol+1
        
        % skip if the element to eliminate is already zero
        if ~R(irow, icol)
            continue
        end
        
        % obtain m x m Givens matrix
        G = GivensMatrix(G, irow-1, irow, R(irow-1:irow, icol));
        
        % apply Givens rotation
        R = G' * R;
        Q = Q  * G;
        
        % reset pre-allocated memory to identity matrix
        G(irow-1:irow, irow-1:irow) = eye(2);
    end
end

end

function G = GivensMatrix(G, i, j, x)
% In-place Givens matrix generation

v = abs(x);
if v(1) >= v(2)
    vtan = v(2) / v(1);
    vcos = 1 / sqrt(vtan*vtan + 1);
    vsin = vcos * vtan;
else
    vcot = v(1) / v(2);
    vsin = 1 / sqrt(vcot*vcot + 1);
    vcos = vsin * vcot;
end

G(i,i) = vcos;
G(i,j) = -vsin; % [vcos -vsin]
G(j,i) = vsin;  % [vsin  vcos]
G(j,j) = vcos;

end

