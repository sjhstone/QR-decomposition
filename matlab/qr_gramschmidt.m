function [Q, R] = qr_gramschmidt(A)

[m, n] = size(A);

% iteration variables
Q = zeros(m, m);
R = zeros(m, n);

% loop over columns of A from left to right
for icol = 1:n
    a = A(:, icol);
    q = a;
    
    % subtract ak's projection onto previously found bases from ak
    for ib = 1:icol-1
        r = Q(:, ib).' * a;
        q = q - r * Q(:, ib);
        R(ib, icol) = r;
    end
    
    % normalize
    q = q ./ norm(q);
    
    % write results
    Q(:, icol) = q;
    R(icol, icol) = a.' * q;
end

end
