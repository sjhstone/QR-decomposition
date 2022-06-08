function [Q, R] = qr_householder(A)

[m, n] = size(A);
% pre-allocation for performance
H = eye(m, m);

Q = eye(m, m);
R = A;

for icol = 1:n
    x = R(icol:m,icol);
    
    % skip if no need to find reflected vector
    if ~nnz(x(2:end))
        continue
    end
    
    v = x;
    v(1) = v(1) + sign(x(1)) * norm(x);
    v = v ./ v(1);
    
    h = eye(m-icol+1) - 2 .* (v*v') ./ (v'*v);
    H(icol:m, icol:m) = h;
    
    R = H * R;
    Q = Q * H';
    
    H(icol:m, icol:m) = eye(m-icol+1);
end

end
