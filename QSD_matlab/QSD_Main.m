function Ud = QSD_Main(U,prnt)
	
	%% Single Elementary Operators and other Global Constants
	
	dim = log2(size(U,1));
	splitpt = 2^(dim-1);
    OU = zeros(splitpt,splitpt);
	O = zeros(2,2);
	I = eye(2);
	X = AP(-pi/2)*Rx(pi);				% [0 1; 1 0] Decomposed 
	H = AP(-pi/2)*Ry(pi/2)*Rx(pi);		% 1/sqrt(2)*[1 1;1 -1] Decomposed
		
	CX = [I O; O X];								% QMUX of {I,X}
	XC = kron(H,H) * CX * kron(H,H);
	SWAP = CX * XC * CX;
	CIX = kron(SWAP,I) * kron(I,CX) * kron(SWAP,I);	
	Tf = [I O O O; O I O O; O O I O; O O O X];		% QMUX of {I,I,I,X}

    %% Decompose U to AB1, CS, AB2
	if (prnt == 1) disp(sprintf('Decompose U using Fat CSD at Level %d\n\n',dim)); end
	[L0,L1,cc,ss,R0,R1] = fatCSD(U);

    %% Decompose AB1 to V, D, W (lower dimension)
 
	AB1 = [L0 zeros(size(L0,1),size(L1,2)); zeros(size(L1,1),size(L0,2)) L1];

	U1 = AB1(1:splitpt,1:splitpt);
	U2 = AB1(splitpt+1:end,splitpt+1:end);
	[v,d,~] = eig(U1*U2');

	V = v;
	D = sqrtm(d);
	W = D*V'*U2;
	decomposedAB1 = [V OU; OU V]*[D OU; OU D']*[W OU; OU W];
	
	if (size(V,1) == 2)
		if (isequal(V,I))
			decomposedV = I;
			if (prnt == 1) disp(sprintf('\t\tDecompose V1 at Level %d: I',dim)); end
		else
			[delta,alpha,theta,beta] = zyz(V);
			decomposedV = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);
			if (prnt == 1) disp(sprintf('\t\tDecompose V1 (using ZYZ) at Level %d: I-Ph(%f) * I-Rz(%f,_) * I-Ry(%f,_) * I-Rz(%f,_)',dim,delta,alpha,theta,beta)); end
		end

		K0 = diag([D(1,1),D(1,1)']);
		K1 = diag([D(2,2),D(2,2)']);
		[~,alpha0,~,beta0] = zyz(K0);
		[~,alpha1,~,beta1] = zyz(K1);
		ab = [alpha0+beta0; alpha1+beta1];
		Minv = inv(genMk(1));
		ar = Minv*ab;
		decomposedD = kron(Rz(ar(1)),I) * XC * kron(Rz(ar(2)),I) * XC;
		if (prnt == 1) disp(sprintf('\t\tDecompose D1 (using ZYZ) at Level %d: Rz(%f,_)-I * XC(_,_) * Rz(%f,_)-I',dim,ar(1),ar(2))); end		

		[delta,alpha,theta,beta] = zyz(W);
		decomposedW = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);
		if (prnt == 1) disp(sprintf('\t\tDecompose W1 (using ZYZ) at Level %d: I-Ph(%f) * I-Rz(%f,_) * I-Ry(%f,_) * I-Rz(%f,_)\n',dim,delta,alpha,theta,beta)); end
			
		decomposedAB1 = kron(I,decomposedV)*decomposedD*kron(I,decomposedW);		
	else
		QSD_Main(V,prnt);
		QSD_Main(D,prnt);
		QSD_Main(W,prnt);
	end    

    %% Decompose CS to Ry, CX
	
	CS = [cc ss; -ss cc];
    % Property Test: cc^2 + ss^2 = eye(size(cc,1))
    
	ta = 2*asin(diag(ss));
	Minv = inv(genMk(dim-1));
	tr = Minv*ta;
	dcs = eye(size(CS,1));
	for i = 1:size(ss,1)
		if (i == size(ss,1))
			posc = 0;
		else
			posc = strfind(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1)),1) - 1;
		end		
		dcs = U_CX(posc,dim-1,dim) * kron([cos(tr(i)/2) sin(tr(i)/2); -sin(tr(i)/2) cos(tr(i)/2)],eye(2^(dim-1))) * dcs;
		if (prnt == 1) disp(sprintf('\t\tDecompose CS (using Grey) at Level %d: CX(%d,%d) * I-Ry(%f,_)\n',dim,posc,dim-1,tr(i))); end
	end
	decomposedCS = dcs;
   
    %% Decompose AB2 to V, D, W (lower dimension)
    
	AB2 = [R0 zeros(size(R0,1),size(R1,2)); zeros(size(R1,1),size(R0,2)) R1];

	U1 = AB2(1:splitpt,1:splitpt);
	U2 = AB2(splitpt+1:end,splitpt+1:end);
	[v,d,~] = eig(U1*U2');

	V = v;
	D = sqrtm(d);
	W = D*V'*U2;
	decomposedAB2 = [V OU; OU V]*[D OU; OU D']*[W OU; OU W];	
	
	if (size(V,1) == 2)
		if (isequal(V,I))
			decomposedV = I;
			if (prnt == 1) disp(sprintf('\t\tDecompose V2 at Level %d: I',dim)); end
		else
			[delta,alpha,theta,beta] = zyz(V);
			decomposedV = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);
			if (prnt == 1) disp(sprintf('\t\tDecompose V2 (using ZYZ) at Level %d: I-Ph(%f) * I-Rz(%f,_) * I-Ry(%f,_) * I-Rz(%f,_)',dim,delta,alpha,theta,beta)); end
		end
		
		K0 = diag([D(1,1),D(1,1)']);
		K1 = diag([D(2,2),D(2,2)']);
		[~,alpha0,~,beta0] = zyz(K0);
		[~,alpha1,~,beta1] = zyz(K1);
		ab = [alpha0+beta0; alpha1+beta1];
		Minv = inv(genMk(1));
		ar = Minv*ab;
		decomposedD = kron(Rz(ar(1)),I) * XC * kron(Rz(ar(2)),I) * XC;
		if (prnt == 1) disp(sprintf('\t\tDecompose D2 (using ZYZ) at Level %d: Rz(%f,_)-I * XC(_,_) * Rz(%f,_)-I',dim,ar(1),ar(2))); end	

		[delta,alpha,theta,beta] = zyz(W);
		decomposedW = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);
		if (prnt == 1) disp(sprintf('\t\tDecompose W2 (using ZYZ) at Level %d: I-Ph(%f) * I-Rz(%f,_) * I-Ry(%f,_) * I-Rz(%f,_)\n',dim,delta,alpha,theta,beta)); end

		decomposedAB2 = kron(I,decomposedV)*decomposedD*kron(I,decomposedW);
	else
		QSD_Main(V,prnt);
		QSD_Main(D,prnt);
		QSD_Main(W,prnt);
	end    

	%% Final Decomposition Testing
	
	% Ud = AB1*CS*AB2;
 	Ud = decomposedAB1 * decomposedCS * decomposedAB2
    
end