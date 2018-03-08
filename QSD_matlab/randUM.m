function U = randUM(d)
% U = randUM(d)
% Returns an Unitary Matrix of size 2^d x 2^d
	n = 2^d;
    M = rand(n)/sqrt(2);
    % M = complex(rand(n),rand(n))/sqrt(2);     % zyz decomposition doesn't work currently for non-SU matrices
    [Q,R] = qr(M);
    Rd = diag(diag(R)./abs(diag(R)));
    U = Q*Rd;
	if det(U)==1	% Only Special Unitary Matrices
		return
	else
		U = randUM(d);
    %U*U';
    %U'*U;
end