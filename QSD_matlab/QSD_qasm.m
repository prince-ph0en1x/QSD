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
	if (prnt == 1) disp(sprintf('Decompose U using Fat CSD at Level %d\n\n',dim)); end
	[L0,L1,cc,ss,R0,R1] = fatCSD(U);
	
	%% Decompose AB2 to V, D, W (lower dimension)

	AB2 = [R0 zeros(size(R0,1),size(R1,2)); zeros(size(R1,1),size(R0,2)) R1];

	U1 = AB2(1:splitpt,1:splitpt);
	U2 = AB2(splitpt+1:end,splitpt+1:end);

	if (max(max(abs(U1+U2))) < 1e-10 || max(max(abs(U1-U2))) < 1e-10)
		% TBD: How this section behaves for dim > 2
		[delta,alpha,theta,beta] = zyz(U1);
		decomposedU1 = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);
		% disp(sprintf('=====>   Decompose U1 (using ZYZ) at Level %d',dim));
		if (prnt == 1)
			disp(sprintf('rz q%d,%f',qbsp(1),-beta));
			disp(sprintf('ry q%d,%f',qbsp(1),-theta));
			disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
			% disp(sprintf('ap q%d,%f',qbsp(1),-delta));
		end
		decomposedAB2 = kron(I,decomposedU1);
		if (max(max(U1+U2)) < 1e-10)
			if (prnt == 1)
				disp(sprintf('rz q%d,%f',qbsp(2),-pi));
				% disp(sprintf('ap q%d,%f',qbsp(2),-pi/2));
			end
			decomposedAB2 = kron(AP(-pi/2)*Rz(pi),I)*decomposedAB2;
		end
	else
%		[v,d,~] = eig(U1*U2');	% MATLAB
		[v,d] = eig(U1,U2');	% OCTAVE
		V = v;
		D = sqrtm(d);
		W = D*V'*U2;
		decomposedAB2 = [V OU; OU V]*[D OU; OU D']*[W OU; OU W];
		
		if (size(W,1) == 2)
			if (isequal(W,I))
				decomposedW = I;
				% disp(sprintf('=====>   Decompose W2 (Identity) at Level %d',dim));
			else
				% disp(sprintf('=====>   Decompose W2 (using ZYZ) at Level %d',dim));
				[delta,alpha,theta,beta] = zyz(W);
				decomposedW = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);			
				if (prnt == 1)
					disp(sprintf('rz q%d,%f',qbsp(1),-beta));
					disp(sprintf('ry q%d,%f',qbsp(1),-theta));
					disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
					% disp(sprintf('ap q%d,%f',qbsp(1),-delta));
				end
			end
		else
			decomposedW = QSD_Main(W,prnt);
		end 		

		ab = 2*log(diag(D))/1i;
		Minv = inv(genMk(dim-1));
		ar = Minv*ab;
		dd = eye(size(AB2,1));
		for i = 1:size(D,1)
			if (i == size(D,1))
				posc = dim-2;
			else
				posc = dim-2 - (strfind(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1)),1) - 1);
			end		
			dd = U_CX(posc,dim-1,dim) * kron(Rz(ar(i)),eye(2^(dim-1))) * dd;
			% disp(sprintf('=====>   Decompose D2 (using Grey) at Level %d',dim));				
			if (prnt == 1)
				disp(sprintf('rz q%d,%f',qbsp(2),-ar(i)));
				disp(sprintf('cnot q%d,q%d',qbsp(1),qbsp(2)));
			end
		end
		decomposedD = dd;

		if (size(V,1) == 2)
			if (isequal(V,I))
				decomposedV = I;
				% disp(sprintf('=====>   Decompose V2 (Identity) at Level %d',dim));
			else
				% disp(sprintf('=====>   Decompose V2 (using ZYZ) at Level %d',dim));
				[delta,alpha,theta,beta] = zyz(V);
				decomposedV = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);
				if (prnt == 1)
					disp(sprintf('rz q%d,%f',qbsp(1),-beta));
					disp(sprintf('ry q%d,%f',qbsp(1),-theta));
					disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
					% disp(sprintf('ap q%d,%f',qbsp(1),-delta));
				end
			end
		else
			decomposedV = QSD_Main(V,prnt);
		end
		
		decomposedAB2 = kron(I,decomposedV)*decomposedD*kron(I,decomposedW);
	end
	decomposedAB2;    

    %% Decompose CS to Ry, CX
    disp(sprintf(''))
	
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
			posc = dim-2 - (strfind(num2str(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1))),num2str(1)) - 1);
		end		
		dcs = U_CX(posc,dim-1,dim) * kron([cos(tr(i)/2) sin(tr(i)/2); -sin(tr(i)/2) cos(tr(i)/2)],eye(2^(dim-1))) * dcs;
		% disp(sprintf('=====>   Decompose CS (using Grey) at Level %d',dim));				
		if (prnt == 1)
			disp(sprintf('ry q%d,%f',qbsp(2),-tr(i)));
			disp(sprintf('cnot q%d,q%d',qbsp(1),qbsp(2)));
		end
	end
	decomposedCS = dcs;
	
    %% Decompose AB1 to V, D, W (lower dimension)
    disp(sprintf(''))
 
	AB1 = [L0 zeros(size(L0,1),size(L1,2)); zeros(size(L1,1),size(L0,2)) L1];

	U1 = AB1(1:splitpt,1:splitpt);
	U2 = AB1(splitpt+1:end,splitpt+1:end);
	
	if (max(max(abs(U1+U2))) < 1e-10 || max(max(abs(U1-U2))) < 1e-10)
		[delta,alpha,theta,beta] = zyz(U1);
		decomposedU1 = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);
		if(max(max(U1+U2)) < 1e-10)
			if (prnt == 1) disp(sprintf('\tDecompose U1 (using ZYZ) at Level %d:\t Ph-I(%f,_) * Rz-I(%f,_) * I-Ph(%f,_) * I-Rz(%f,_) * I-Ry(%f,_) * I-Rz(%f,_)',dim,-pi/2,pi,delta,alpha,theta,beta)); end
			decomposedAB1 = kron(AP(-pi/2)*Rz(pi),I)*kron(I,decomposedU1);
		else
			if (prnt == 1) disp(sprintf('\tDecompose U1 (using ZYZ) at Level %d:\t I-Ph(%f,_) * I-Rz(%f,_) * I-Ry(%f,_) * I-Rz(%f,_)',dim,delta,alpha,theta,beta)); end
			decomposedAB1 = kron(I,I)*kron(I,decomposedU1);
		end
	else

%		[v,d,~] = eig(U1*U2');	% MATLAB
		[v,d] = eig(U1,U2');	% OCTAVE
		V = v;
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
				decomposedW = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);				
				if (prnt == 1)
					disp(sprintf('rz q%d,%f',qbsp(1),-beta));
					disp(sprintf('ry q%d,%f',qbsp(1),-theta));
					disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
					% disp(sprintf('ap q%d,%f',qbsp(1),-delta));
				end
			end
		else
			decomposedW = QSD_Main(W,prnt);
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
				posc = dim-2 - (strfind(num2str(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1))),num2str(1)) - 1);
			end		
			dd = U_CX(posc,dim-1,dim) * kron(Rz(ar(i)),eye(2^(dim-1))) * dd;
			% disp(sprintf('=====>   Decompose D1 (using Grey) at Level %d',dim));				
			if (prnt == 1)
				disp(sprintf('rz q%d,%f',qbsp(2),-ar(i)));
				disp(sprintf('cnot q%d,q%d',qbsp(1),qbsp(2)));
			end
		end
		[D O; O D'];
		decomposedD = dd;

		if (size(V,1) == 2)
			if (isequal(V,I))
				decomposedV = I;
				% disp(sprintf('=====>   Decompose V1 (Identity) at Level %d',dim));
			else
				% disp(sprintf('=====>   Decompose V1 (using ZYZ) at Level %d',dim));
				[delta,alpha,theta,beta] = zyz(V);
				decomposedV = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta);
				if (prnt == 1)
					disp(sprintf('rz q%d,%f',qbsp(1),-beta));
					disp(sprintf('ry q%d,%f',qbsp(1),-theta));
					disp(sprintf('rz q%d,%f',qbsp(1),-alpha));
					% disp(sprintf('ap q%d,%f',qbsp(1),-delta));
				end
			end
		else
			decomposedV = QSD_Main(V,prnt);
		end
		decomposedV;
		
		decomposedAB1 = kron(I,decomposedV)*decomposedD*kron(I,decomposedW);
	end
	decomposedAB1;

	%% Final Decomposition Testing
	disp(sprintf(''))

	Ud = AB1*CS*AB2;
	Ud = decomposedAB1 * decomposedCS * decomposedAB2;
	
end