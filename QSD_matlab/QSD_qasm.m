function Ud = QSD_qasm(U,prnt,qbsp)
	
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
%	if (prnt == 1) disp(sprintf('Decompose U using Fat CSD at Level %d\n\n',dim)); end
	[L0,L1,cc,ss,R0,R1] = fatCSD(U);
	
	%% Decompose AB2 to V, D, W (lower dimension)

	AB2 = [R0 zeros(size(R0,1),size(R1,2)); zeros(size(R1,1),size(R0,2)) R1];

	U1 = AB2(1:splitpt,1:splitpt);
	U2 = AB2(splitpt+1:end,splitpt+1:end);
		
	if (max(max(abs(U1+U2))) < 1e-10 || max(max(abs(U1-U2))) < 1e-10)
		% TBD: How this section behaves for dim > 2
		if (dim > 2)
			'oop'
		end
		[delta,alpha,theta,beta] = zyz(U1);
		decomposedU1 = Rz(alpha)*Ry(theta)*Rz(beta);
		% disp(sprintf('=====>   Decompose U1 (using ZYZ) at Level %d',dim));
		if (prnt > 1)
			disp(sprintf('rz q%d,%f',qbsp(1),-beta));
			disp(sprintf('ry q%d,%f',qbsp(1),-theta));
			disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
		end
		if (prnt == 2)
			decomposedU1 = AP(delta)*decomposedU1;
			disp(sprintf('ap q%d,%f',qbsp(1),-delta));
		end
		decomposedAB2 = kron(I,decomposedU1);
		if (max(max(U1+U2)) < 1e-10)
			decomposedAB2 = kron(Rz(pi),I)*decomposedAB2;
			if (prnt > 1)
				disp(sprintf('rz q%d,%f',qbsp(2),-pi));
			end
			if (prnt == 2)
				decomposedAB2 = kron(AP(-pi/2),I)*decomposedAB2;
				disp(sprintf('ap q%d,%f',qbsp(2),pi/2));
			end
		end
	else
		[v,d,~] = eig(U1*U2');	% MATLAB
% 		[v,d] = eig(complex(U1),complex(U2'));	% OCTAVE
% 		[v,d] = eig(U1,U2');	% OCTAVE
		V = v;
		if dim == 3
			% TBD: Automate this
			% 1 & 2 Eigenvalues in d are repeated, thus V*V' is not I. Adjustment needed
			% https://nl.mathworks.com/matlabcentral/answers/214557-eigenvectors-are-not-orthogonal-for-some-skew-symmetric-matrices-why
			V(:,[1,2]) = orth(V(:,[1,2]));
		end
		D = sqrtm(d);
		W = D*V'*U2;
		decomposedAB2 = [V OU; OU V]*[D OU; OU D']*[W OU; OU W];
		
% 		if dim == 3 
% 			disp(sprintf('+++++++++++++++++++++++++++++++++')); 
% 		end
		
		if (size(W,1) == 2)
			if (isequal(W,I))
				decomposedW = I;
				% disp(sprintf('=====>   Decompose W2 (Identity) at Level %d',dim));
			else
				% disp(sprintf('=====>   Decompose W2 (using ZYZ) at Level %d',dim));
				[delta,alpha,theta,beta] = zyz(W);
				decomposedW = Rz(alpha)*Ry(theta)*Rz(beta);	
				if (prnt > 1)
					disp(sprintf('rz q%d,%f',qbsp(1),-beta));
					disp(sprintf('ry q%d,%f',qbsp(1),-theta));
					disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
				end
				if (prnt == 2)
					decomposedW = AP(delta)*decomposedW;	
					disp(sprintf('ap q%d,%f',qbsp(1),-delta));
				end
			end
		else
			decomposedW = QSD_qasm(W,prnt,qbsp(1:end-1));
%			maxerrW3 = max(max(W-decomposedW))
		end
	
% 		if dim == 3 
% 			disp(sprintf('+++++++++++++++++++++++++++++++++')); 
% 		end

		ab = 2*log(diag(D))/1i;
		Minv = inv(genMk(dim-1));
		ar = Minv*ab;
		dd = eye(size(AB2,1));
		for i = 1:size(D,1)
			if (i == size(D,1))
				posc = dim-2;
			else
				[~,idx] = find(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1)),1);
				posc = dim-2 - (idx - 1);
%				posc = dim-2 - (strfind(num2str(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1))),num2str(1)) - 1) 
			end
			dd = U_CX(posc,dim-1,dim) * kron(Rz(ar(i)),eye(2^(dim-1))) * dd;
			% disp(sprintf('=====>   Decompose D2 (using Grey) at Level %d',dim));				
			if (prnt > 1)
				disp(sprintf('rz q%d,%f',qbsp(end),-ar(i)));
				disp(sprintf('cnot q%d,q%d',posc,dim-1));
			end
		end
		decomposedD = dd;

% 		if dim == 3 
% 			disp(sprintf('+++++++++++++++++++++++++++++++++')); 
% 		end
		
		if (size(V,1) == 2)
			if (isequal(V,I))
				decomposedV = I;
				% disp(sprintf('=====>   Decompose V2 (Identity) at Level %d',dim));
			else
				% disp(sprintf('=====>   Decompose V2 (using ZYZ) at Level %d',dim));
				[delta,alpha,theta,beta] = zyz(V);
				decomposedV = Rz(alpha)*Ry(theta)*Rz(beta);
				if (prnt > 1)
					disp(sprintf('rz q%d,%f',qbsp(1),-beta));
					disp(sprintf('ry q%d,%f',qbsp(1),-theta));
					disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
				end
				if (prnt == 2)
					decomposedV = AP(delta)*decomposedV;
					disp(sprintf('ap q%d,%f',qbsp(1),-delta));
				end
			end
		else
			decomposedV = QSD_qasm(V,prnt,qbsp(1:end-1));
		end
		
		decomposedAB2 = kron(I,decomposedV)*decomposedD*kron(I,decomposedW);
	end
	decomposedAB2; 
	
% 	if dim == 3
% 		kron(I,decomposedW)
% 	end

    %% Decompose CS to Ry, CX
    if (prnt > 1) disp(sprintf('')); end
% 	if dim == 3 
% 		disp(sprintf('============================================')); 
% 	end
	
	CS = [cc ss; -ss cc];
    % Property Test: cc^2 + ss^2 = eye(size(cc,1))
    
	ta = 2*asin(diag(ss));
	Minv = inv(genMk(dim-1));
	tr = Minv*ta;
	dcs = eye(size(CS,1));
	for i = 1:size(ss,1)
		if (i == size(ss,1))
			posc = dim-2;
		else
			[~,idx] = find(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1)),1);
			posc = dim-2 - (idx - 1);
%			posc = dim-2 - (strfind(num2str(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1))),num2str(1)) - 1);
		end		
		dcs = kron(Ry(tr(i)),eye(2^(dim-1))) * dcs;
		dcs = U_CX(posc,dim-1,dim) * dcs;
		% disp(sprintf('=====>   Decompose CS (using Grey) at Level %d',dim));				
		if (prnt > 1)
			disp(sprintf('ry q%d,%f',qbsp(end),-tr(i)));
			disp(sprintf('cnot q%d,q%d',posc,dim-1));
		end
	end
	decomposedCS = dcs;
	
    %% Decompose AB1 to V, D, W (lower dimension)
    if (prnt > 1) disp(sprintf('')); end
% 	if dim == 3 
% 		disp(sprintf('============================================')); 
% 	end
	
	AB1 = [L0 zeros(size(L0,1),size(L1,2)); zeros(size(L1,1),size(L0,2)) L1];

	U1 = AB1(1:splitpt,1:splitpt);
	U2 = AB1(splitpt+1:end,splitpt+1:end);

	if (max(max(abs(U1+U2))) < 1e-10 || max(max(abs(U1-U2))) < 1e-10)
		% TBD: How this section behaves for dim > 2
		if (dim > 2)
			'oop'
		end
		[delta,alpha,theta,beta] = zyz(U1);
		decomposedU1 = Rz(alpha)*Ry(theta)*Rz(beta);
		% disp(sprintf('=====>   Decompose U1 (using ZYZ) at Level %d',dim));
		if (prnt > 1)
			disp(sprintf('rz q%d,%f',qbsp(1),-beta));
			disp(sprintf('ry q%d,%f',qbsp(1),-theta));
			disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
		end
		if (prnt == 2)
			decomposedU1 = AP(delta)*decomposedU1;
			disp(sprintf('ap q%d,%f',qbsp(1),-delta));
		end
		decomposedAB1 = kron(I,decomposedU1);
		if (max(max(U1+U2)) < 1e-10)
			decomposedAB1 = kron(Rz(pi),I)*decomposedAB1;
			if (prnt > 1)
				disp(sprintf('rz q%d,%f',qbsp(2),-pi));
			end
			if (prnt == 2)
				decomposedAB1 = kron(AP(-pi/2),I)*decomposedAB1;
				disp(sprintf('ap q%d,%f',qbsp(2),pi/2));
			end
		end
	else
		[v,d,~] = eig(U1*U2');	% MATLAB
% 		[v,d] = eig(U1,U2');	% OCTAVE
		V = v;
		if dim == 3
			d;
			% TBD: Automate this
			% 1 & 2 Eigenvalues in d are repeated, thus V*V' is not I. Adjustment needed
			% https://nl.mathworks.com/matlabcentral/answers/214557-eigenvectors-are-not-orthogonal-for-some-skew-symmetric-matrices-why
% 			V(:,[1,2]) = orth(V(:,[1,2]));
		end
		D = sqrtm(d);
		W = D*V'*U2;
		decomposedAB1 = [V OU; OU V]*[D OU; OU D']*[W OU; OU W];

		if (size(W,1) == 2)
			if (isequal(W,I))
				decomposedW = I;
				% disp(sprintf('=====>   Decompose W1 (Identity) at Level %d',dim));
			else
				% disp(sprintf('=====>   Decompose W1 (using ZYZ) at Level %d',dim));
				[delta,alpha,theta,beta] = zyz(W);
				decomposedW = Rz(alpha)*Ry(theta)*Rz(beta);				
				if (prnt > 1)
					disp(sprintf('rz q%d,%f',qbsp(1),-beta));
					disp(sprintf('ry q%d,%f',qbsp(1),-theta));
					disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
				end
				if (prnt == 2)
					decomposedW = AP(delta)*decomposedW;				
					disp(sprintf('ap q%d,%f',qbsp(1),-delta));
				end
			end
		else
			decomposedW = QSD_qasm(W,prnt,qbsp(1:end-1));
		end    
		decomposedW;

		ab = 2*log(diag(D))/1i;
		Minv = inv(genMk(dim-1));
		ar = Minv*ab;
		dd = eye(size(AB1,1));
		for i = 1:size(D,1)
			if (i == size(D,1))
				posc = dim-2;
			else
				[~,idx] = find(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1)),1);
				posc = dim-2 - (idx - 1);
%				posc = dim-2 - (strfind(num2str(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1))),num2str(1)) - 1);
			end
			dd = U_CX(posc,dim-1,dim) * kron(Rz(ar(i)),eye(2^(dim-1))) * dd;
			% disp(sprintf('=====>   Decompose D1 (using Grey) at Level %d',dim));				
			if (prnt > 1)
				disp(sprintf('rz q%d,%f',qbsp(end),-ar(i)));
				disp(sprintf('cnot q%d,q%d',posc,dim-1));
			end
		end
		decomposedD = dd;
	
		if (size(V,1) == 2)
			if (isequal(V,I))
				decomposedV = I;
				% disp(sprintf('=====>   Decompose V1 (Identity) at Level %d',dim));
			else
				% disp(sprintf('=====>   Decompose V1 (using ZYZ) at Level %d',dim));
				[delta,alpha,theta,beta] = zyz(V);
				decomposedV = Rz(alpha)*Ry(theta)*Rz(beta);
				if (prnt > 1)
					disp(sprintf('rz q%d,%f',qbsp(1),-beta));
					disp(sprintf('ry q%d,%f',qbsp(1),-theta));
					disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
				end
				if (prnt == 2)
					decomposedV = AP(delta)*decomposedV;
					disp(sprintf('ap q%d,%f',qbsp(1),-delta));
				end
			end
		else
			decomposedV = QSD_qasm(V,prnt,qbsp(1:end-1));
		end
		decomposedV;
		
		decomposedAB1 = kron(I,decomposedV)*decomposedD*kron(I,decomposedW);
	end
	decomposedAB1;

	%% Final Decomposition Testing
	if (prnt > 1) disp(sprintf('')); end

	Ud = AB1*CS*AB2;
	Ud = decomposedAB1 * decomposedCS * decomposedAB2;
% 	if dim == 3
% 		AB2
% 		decomposedAB2
% 		maxerr2 = max(max(abs(Ud)-abs(U)))
% 	end
	
end