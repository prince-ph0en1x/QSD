function U = randUM(d)
% U = randUM(d)
% Returns an Unitary Matrix of size 2^d x 2^d
	n = 2^d;
%	M = rand(n)/sqrt(2);	% If only Real
    M = complex(rand(n),rand(n))/sqrt(2);
    [Q,R] = qr(M);
    Rd = diag(diag(R)./abs(diag(R)));
    U = Q*Rd;
% 	U = U/(det(U))^(1/n);	% If only SU
end