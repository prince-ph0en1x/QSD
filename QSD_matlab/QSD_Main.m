function QSD_Main(dim)

    splitpt = 2^(dim-1);
    
    U = randUM(dim)
	
	%% Single Elementary Operators

	X = [0 1; 1 0];
	CX = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
	Tf = [eye(6,6) zeros(6,2); zeros(2,6) X];
	H = 1/sqrt(2)*[1 1;1 -1];
	CIX = kron(CX,eye(2)) * kron(kron(H,H),eye(2)) * kron(CX,eye(2)) * kron(kron(H,H),eye(2)) * kron(CX,eye(2)) * kron(eye(2),CX) * kron(CX,eye(2)) * kron(kron(H,H),eye(2)) * kron(CX,eye(2)) * kron(kron(H,H),eye(2)) * kron(CX,eye(2));
	XC = [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0];
% 	XC = [X zeros(2,2); zeros(2,2) eye(2)]

    %% Decompose U to AB1, CS, AB2
    [L0,L1,cc,ss,R0,R1] = fatCSD(U);
    AB1 = [L0 zeros(size(L0,1),size(L1,2)); zeros(size(L1,1),size(L0,2)) L1];
    CS = [cc ss; -ss cc];
    AB2 = [R0 zeros(size(R0,1),size(R1,2)); zeros(size(R1,1),size(R0,2)) R1];
    % cc^2 + ss^2
    decomposedU = AB1*CS*AB2;

    %% Decompose AB1 to V1, D1, W1 (lower dimension)
    U1 = AB1(1:splitpt,1:splitpt);
    U2 = AB1(splitpt+1:end,splitpt+1:end);
    [v,d,~] = eig(U1*U2');
    V1 = v;                     % Recursive breakdown CSD 
    D1 = sqrtm(d);              % Recursive breakdown
    W1 = D1*V1'*U2;             % Recursive breakdown CSD
    %AB1
    decomposedAB1 = [V1 zeros(splitpt,splitpt); zeros(splitpt,splitpt) V1]*[D1 zeros(splitpt,splitpt); zeros(splitpt,splitpt) D1']*[W1 zeros(splitpt,splitpt); zeros(splitpt,splitpt) W1];
    
    %%
	t1 = 2*asin(CS(1,3));
	t2 = 2*asin(CS(2,4));	
	r1q = [cos(t1/2) sin(t1/2); -sin(t1/2) cos(t1/2)];
	r2q = [cos(t2/2) sin(t2/2); -sin(t2/2) cos(t2/2)];
	rt1 = (t1+t2)/2;
    rt2 = (t1-t2)/2;
    r1q = [cos(rt1/2) sin(rt1/2); -sin(rt1/2) cos(rt1/2)];
    r2q = [cos(rt2/2) sin(rt2/2); -sin(rt2/2) cos(rt2/2)];
	decomposedCS = XC * kron(r2q,eye(2)) * XC * kron(r1q,eye(2));
	    
    %% Decompose AB2 to V2, D2, W2 (lower dimension)
    U1 = AB2(1:splitpt,1:splitpt);
    U2 = AB2(splitpt+1:end,splitpt+1:end);
    [v,d,~] = eig(U1*U2');
    V2 = v;                     % Recursive breakdown CSD
    D2 = sqrtm(d);               % Recursive breakdown
    W2 = D2*V2'*U2;             % Recursive breakdown CSD
    decomposedAB2 = [V2 zeros(splitpt,splitpt); zeros(splitpt,splitpt) V2]*[D2 zeros(splitpt,splitpt); zeros(splitpt,splitpt) D2']*[W2 zeros(splitpt,splitpt); zeros(splitpt,splitpt) W2];
	
	%%
	
	decomposedAB1*decomposedCS*decomposedAB2
	
    %% Decompose Dx of size 2x2
%     D2th = log(D2(2,2)/D2(1,1))/i;
%     D2dt = (log(D2(2,2)*D2(1,1))-log(cos(0)^2))/(2*i);   % Arbitrary Phase
%     % Rza = [exp(-i*D2th/2) 0; 0 exp(i*D2th/2)];
%     % exp(i*D2dt)*Rza;
%     disp(sprintf('RZ q0,%0.3f',D2th))
%     disp(sprintf('ArbPhase %0.3f',D2dt))
    
end