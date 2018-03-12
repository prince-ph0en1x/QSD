function Ud = QSD_Main(U)
	
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
	
    disp(sprintf('Decompose U using Fat CSD at Level %d',dim))
	[L0,L1,cc,ss,R0,R1] = fatCSD(U);
    % Property Test: cc^2 + ss^2 = eye(size(cc,1))

    %% Decompose AB1 to V, D, W (lower dimension)
    
	AB1 = [L0 zeros(size(L0,1),size(L1,2)); zeros(size(L1,1),size(L0,2)) L1];
    
	U1 = AB1(1:splitpt,1:splitpt);
    U2 = AB1(splitpt+1:end,splitpt+1:end);
    [v,d,~] = eig(U1*U2');
	
    V = v; 
    D = sqrtm(d);
    W = D*V'*U2;
	
	if (size(V,1) == 2)
		disp(sprintf('\t\tDecompose V1 using ZYZ at Level %d',dim))
		[~,Rza,Ryt,Rzb] = zyz(V);
		if (det(V) == 1)		% TBD: Handle non-SU
% 			[Rza Ryt Rzb]
		end
		disp(sprintf('\t\tDecompose D1 using ZYZ at Level %d',dim))
		[~,Rza,Ryt,Rzb] = zyz(D);
		if (det(D) == 1)		% TBD: Handle non-SU
% 			[Rza Ryt Rzb]
		end
		disp(sprintf('\t\tDecompose W1 using ZYZ at Level %d',dim))
		[~,Rza,Ryt,Rzb] = zyz(W);
		if (det(W) == 1)		% TBD: Handle non-SU
% 			[Rza Ryt Rzb]
		end
	else
		QSD_Main(V);
		QSD_Main(D);
		QSD_Main(W);
	end
	
    decomposedAB1 = [V OU; OU V]*[D OU; OU D']*[W OU; OU W];
    
    %%
	
	CS = [cc ss; -ss cc];
    
	ta = 2*asin(diag(ss));
	Minv = inv(genMk(dim-1));
	tr = Minv*ta;
	% tr = [tr(1); tr(4); tr(3); tr(2)];	% TBD: Why is this required? General unjumble scheme
	dcs = eye(size(CS,1));
	for i = 1:size(ss,1)
		if (i == size(ss,1))
			posc = 0;
		else
			posc = strfind(sprintf(dec2bin(bin2gray(i-1),dim-1)) ~= sprintf(dec2bin(bin2gray(i),dim-1)),1) - 1;
		end		
		dcs = U_CX(posc,dim-1,dim) * kron([cos(tr(i)/2) sin(tr(i)/2); -sin(tr(i)/2) cos(tr(i)/2)],eye(2^(dim-1))) * dcs;
	end
	disp(sprintf('\t\tDecompose CS using Grey at Level %d',dim))
	decomposedCS = dcs;
	    
    %% Decompose AB2 to V, D, W (lower dimension)
    
	AB2 = [R0 zeros(size(R0,1),size(R1,2)); zeros(size(R1,1),size(R0,2)) R1];
    
	U1 = AB2(1:splitpt,1:splitpt);
    U2 = AB2(splitpt+1:end,splitpt+1:end);
    [v,d,~] = eig(U1*U2');
	
    V = v; 
    D = sqrtm(d);
    W = D*V'*U2;
	
	if (size(V,1) == 2)
		disp(sprintf('\t\tDecompose V1 using ZYZ at Level %d',dim))
		[~,Rza,Ryt,Rzb] = zyz(V);
		if (det(V) == 1)		% TBD: Handle non-SU
% 			[Rza Ryt Rzb]
		end
		disp(sprintf('\t\tDecompose D1 using ZYZ at Level %d',dim))
		[~,Rza,Ryt,Rzb] = zyz(D);
		if (det(D) == 1)		% TBD: Handle non-SU
% 			[Rza Ryt Rzb]
		end
		disp(sprintf('\t\tDecompose W1 using ZYZ at Level %d',dim))
		[~,Rza,Ryt,Rzb] = zyz(W);
		if (det(W) == 1)		% TBD: Handle non-SU
% 			[Rza Ryt Rzb]
		end
	else
		QSD_Main(V);
		QSD_Main(D);
		QSD_Main(W);
	end
	
    decomposedAB2 = [V OU; OU V]*[D OU; OU D']*[W OU; OU W];
	
	%% Final Decomposition Testing
	
    % decomposedU = AB1*CS*AB2;
	Ud = decomposedAB1 * decomposedCS * decomposedAB2;
    
end