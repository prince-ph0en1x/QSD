% Quantum Shannon Decomposition
% 06-03-2018

clear all
clc

%% Single Elementary Operators

X = [0 1; 1 0];
CX = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
Tf = [eye(6,6) zeros(6,2); zeros(2,6) X];
H = 1/sqrt(2)*[1 1;1 -1];
CIX = kron(CX,eye(2)) * kron(kron(H,H),eye(2)) * kron(CX,eye(2)) * kron(kron(H,H),eye(2)) * kron(CX,eye(2)) * kron(eye(2),CX) * kron(CX,eye(2)) * kron(kron(H,H),eye(2)) * kron(CX,eye(2)) * kron(kron(H,H),eye(2)) * kron(CX,eye(2));
% XC = [X zeros(2,2); zeros(2,2) eye(2)];

%% Theorem 1 -- ZYZ Decomposition

	U = randUM(1)
% 	U = [0.0000 - 0.7071i   0.7071 + 0.0000i;   0.7071 + 0.0000i   0.0000 - 0.7071i]
% 	det(U)
% 	U*U'
	[delta,alpha,theta,beta] = zyz(U);
	decomposedU = AP(delta)*Rz(alpha)*Ry(theta)*Rz(beta)
    
%% Theorem 4 -- Demultiplexing A Singly Multiplexed Ry or Rz

% 	t1 = rand*2*pi-pi;
% 	t2 = rand*2*pi-pi;
%	ry1 = [cos(t1/2) sin(t1/2); -sin(t1/2) cos(t1/2)];
% 	ry2 = [cos(t2/2) sin(t2/2); -sin(t2/2) cos(t2/2)];
% 	tgtM = [ry1 zeros(2,2); zeros(2,2) ry2]
% 
% 	decomposed = cRyz(t1,t2)			% works

%% Theorem 5 -- Decompositions of a Two-Qubit Multiplexor

% 	U0 = randUM(1);
% 	U1 = randUM(1);
% 	U = [U0 zeros(2,2); zeros(2,2) U1]
% 
% 	W = U0*U1';
% 	decomposed = [W zeros(2,2); zeros(2,2) eye(2)] * kron(eye(2),U1) 
% 
% 	det(W)
% 	[Phd,alpha,theta,beta] = zyz(W);
% 	WA = [exp(1i*alpha/2) 0; 0 exp(-1i*alpha/2)]*[cos(theta/4) sin(theta/4); -sin(theta/4) cos(theta/4)];
% 	WB = [cos(-theta/4) sin(-theta/4); -sin(-theta/4) cos(-theta/4)]*[exp(1i*(-(alpha+beta)/2)/2) 0; 0 exp(-1i*(-(alpha+beta)/2)/2)];
% 	WC = [exp(1i*((beta-alpha)/2)/2) 0; 0 exp(-1i*((beta-alpha)/2)/2)];
% 	% W0 = WA*WB*WC;
% 	% W1 = WA*X*WB*X*WC;
% 	CX = [X zeros(2,2); zeros(2,2) eye(2)]
% 	decomposed = kron(eye(2),WC) * CX * kron(eye(2),WB) * CX * kron(eye(2),WA) * kron(eye(2),U1)

%% Theorem 6 -- ZYZ Decomposition for Single-Data-Bit Multiplexors | Theorem 8 -- Demultiplexing Multiplexed Rk Gates, k = y,z

% % 1-bit Control
% 
% 	d = 1;
% 	n = 2^d;
% 	U1 = randUM(1);
% 	U2 = randUM(1);
% 	[det(U1) det(U2)]
% 	U = [U1 zeros(n,n); zeros(n,n) U2]
% 
% 	[Phd1,Rza1,Ryt1,Rzb1] = zyz(U1); 
% 	[Phd2,Rza2,Ryt2,Rzb2] = zyz(U2);
% 
% 	t1 = -2*log(Rza1(1,1))/1i;
% 	t2 = -2*log(Rza2(1,1))/1i;
% 	decomposed1 = cRyz(t1,t2);
% 	% tgtM = [Rza1 zeros(2,2); zeros(2,2) Rza2]
% 
% 	t1 = 2*asin(Ryt1(1,2));
% 	t2 = 2*asin(Ryt2(1,2));
% 	decomposed2 = cRyz(t1,t2);
% 	% tgtM = [Ryt1 zeros(2,2); zeros(2,2) Ryt2]
% 
% 	t1 = -2*log(Rzb1(1,1))/1i;
% 	t2 = -2*log(Rzb2(1,1))/1i;
% 	decomposed3 = cRyz(t1,t2);
% 	% tgtM = [Rzb1 zeros(2,2); zeros(2,2) Rzb2]
% 
% 	decomposed = decomposed3 * decomposed2 * decomposed1
    
% %	2-bit Control
% 
% 	d = 1;
% 	n = 2^d;
% 	U1 = randUM(1);
% 	U2 = randUM(1);
% 	U3 = randUM(1);
% 	U4 = randUM(1);
% 	[det(U1) det(U2) det(U3) det(U4)];
% 	U = [U1 zeros(n,3*n); zeros(n,n) U2 zeros(n,2*n); zeros(n,2*n) U3 zeros(n,n); zeros(n,3*n) U4]
% 
% 	[Phd1,Rza1,Ryt1,Rzb1] = zyz(U1); 
% 	[Phd2,Rza2,Ryt2,Rzb2] = zyz(U2);
% 	[Phd3,Rza3,Ryt3,Rzb3] = zyz(U3);
% 	[Phd4,Rza4,Ryt4,Rzb4] = zyz(U4);
% 
% 	% R = [Rza1 zeros(n,3*n); zeros(n,n) Rza2 zeros(n,2*n); zeros(n,2*n) Rza3 zeros(n,n); zeros(n,3*n) Rza4]
% 	t1 = -2*log(Rza1(1,1))/1i;
% 	t2 = -2*log(Rza2(1,1))/1i;
% 	t3 = -2*log(Rza3(1,1))/1i;
% 	t4 = -2*log(Rza4(1,1))/1i;
% 	decomposed1a = cRyz(t1,t2);
% 	decomposed1b = cRyz(t3,t4);    
% 	decomposed1 = Tf * kron(eye(2),decomposed1a) * Tf * kron(eye(2),decomposed1b);
% 
% 	% R = [Ryt1 zeros(n,3*n); zeros(n,n) Ryt2 zeros(n,2*n); zeros(n,2*n) Ryt3 zeros(n,n); zeros(n,3*n) Ryt4]
% 	t1 = 2*asin(Ryt1(1,2));
% 	t2 = 2*asin(Ryt2(1,2));
% 	t3 = 2*asin(Ryt3(1,2));
% 	t4 = 2*asin(Ryt4(1,2));
% 	Minv = inv(genMk(2));
% 	tr = Minv*[t1;t2;t3;t4];
% 	r1q = kron(eye(4),[cos(tr(1)/2) sin(tr(1)/2); -sin(tr(1)/2) cos(tr(1)/2)]);
% 	r2q = kron(eye(4),[cos(tr(2)/2) sin(tr(2)/2); -sin(tr(2)/2) cos(tr(2)/2)]);
% 	r3q = kron(eye(4),[cos(tr(3)/2) sin(tr(3)/2); -sin(tr(3)/2) cos(tr(3)/2)]);
% 	r4q = kron(eye(4),[cos(tr(4)/2) sin(tr(4)/2); -sin(tr(4)/2) cos(tr(4)/2)]);
% 	% 	r1 = [cos(tr(1)/2) sin(tr(1)/2); -sin(tr(1)/2) cos(tr(1)/2)]
% 	%   r2 = [cos(tr(2)/2) sin(tr(2)/2); -sin(tr(2)/2) cos(tr(2)/2)]
% 	%   r3 = [cos(tr(3)/2) sin(tr(3)/2); -sin(tr(3)/2) cos(tr(3)/2)]
% 	%   r4 = [cos(tr(4)/2) sin(tr(4)/2); -sin(tr(4)/2) cos(tr(4)/2)]
% 	% 	A = r4*r3*r2*r1
% 	% 	B = r4*(X*r3*X)*(X*r2*X)*r1
% 	% 	C = (X*r4*X)*(X*r3*X)*r2*r1
% 	% 	D = (X*r4*X)*r3*(X*r2*X)*r1
% 	ctrls = bin2gray([0:3]);
% 	decomposed2 = CIX * r4q * kron(eye(2),CX) * r3q * CIX * r2q * kron(eye(2),CX) * r1q;
% 
% 	% R = [Rzb1 zeros(n,3*n); zeros(n,n) Rzb2 zeros(n,2*n); zeros(n,2*n) Rzb3 zeros(n,n); zeros(n,3*n) Rzb4]
% 	t1 = -2*log(Rzb1(1,1))/1i;
% 	t2 = -2*log(Rzb2(1,1))/1i;
% 	t3 = -2*log(Rzb3(1,1))/1i;
% 	t4 = -2*log(Rzb4(1,1))/1i;
% 	decomposed3a = cRyz(t1,t2);
% 	decomposed3b = cRyz(t3,t4);    
% 	decomposed3 = Tf * kron(eye(2),decomposed3a) * Tf * kron(eye(2),decomposed3b);
% 
% 	decomposed = decomposed3 * decomposed2 * decomposed1

%%



%%
% CUd = [eye(2) zeros(2,2); zeros(2,2) Phd] * [eye(2) zeros(2,2); zeros(2,2) Rza] * [eye(2) zeros(2,2); zeros(2,2) Ryt] * [eye(2) zeros(2,2); zeros(2,2) Rzb]
% 
% 
% 
% WA = [e^(i*alpha/2) 0; 0 e^(-i*alpha/2)]*[cos(theta/4) sin(theta/4); -sin(theta/4) cos(theta/4)];
% WB = [cos(-theta/4) sin(-theta/4); -sin(-theta/4) cos(-theta/4)]*[e^(i*(-(alpha+beta)/2)/2) 0; 0 e^(-i*(-(alpha+beta)/2)/2)];
% WC = [e^(i*((beta-alpha)/2)/2) 0; 0 e^(-i*((beta-alpha)/2)/2)];
% W0 = WA*WB*WC;
% W1 = WA*X*WB*X*WC;
% 
% CW = [eye(2) zeros(2,2); zeros(2,2) W1]
% CWd = kron(eye(2),WA)*CX*kron(eye(2),WB)*CX*kron(eye(2),WC)
% 
% 
% V1 = sqrtm(W1);
% V1dag = conj(V1)';
% 
% nC = 3;
% nCW = W1;
% for i = 1:nC
% 	s = size(nCW,1);
% 	nCW = [eye(2) zeros(2,s); zeros(s,2) nCW];
% end
% nCW
% 
% 


