function U = U_CX(posc,post,spc)

	X = [0 1; 1 0];
	H = 1/sqrt(2)*[1 1;1 -1];

	if (posc > post)
% 		U = kron(eye(2^(spc-1-posc)),kron(H,eye(2^(posc))));
% 		U = kron(eye(2^(spc-1-post)),kron(H,eye(2^(post)))) * U;
		I = eye(2^(posc-post));
		O = zeros(2^(posc-post));
		br = kron(eye(2^(posc-post-1)),X);
		CX = [I O;O br];
		U = kron(eye(2^(spc-posc-1)),kron(CX,eye(2^post)));
% 		U = kron(eye(2^(spc-1-posc)),kron(H,eye(2^(posc)))) * U;
% 		U = kron(eye(2^(spc-1-post)),kron(H,eye(2^(post)))) * U;
	else
		U = kron(eye(2^(spc-1-posc)),kron(H,eye(2^(posc))));
		U = kron(eye(2^(spc-1-post)),kron(H,eye(2^(post)))) * U;
		I = eye(2^(post-posc));
		O = zeros(2^(post-posc));
		br = kron(eye(2^(post-posc-1)),X);
		CX = [I O;O br];
		U = kron(eye(2^(spc-post-1)),kron(CX,eye(2^posc))) * U;
		U = kron(eye(2^(spc-1-posc)),kron(H,eye(2^(posc)))) * U;
		U = kron(eye(2^(spc-1-post)),kron(H,eye(2^(post)))) * U;
	end

end