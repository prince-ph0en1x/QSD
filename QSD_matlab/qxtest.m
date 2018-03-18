clear all
clc

U =[-0.2800   -0.6400   -0.6400   -0.3200
   -0.6400    0.6800   -0.3200   -0.1600
   -0.6400   -0.3200    0.6800   -0.1600
   -0.3200   -0.1600   -0.1600    0.9200];

AB1 = [	0.89443   0.44721   0.00000   0.00000
		0.44721  -0.89443   0.00000   0.00000
		0.00000   0.00000   0.89443  -0.44721
		0.00000   0.00000   0.44721   0.89443];
   
CS = [	0.60000  -0.00000  -0.80000   0.00000
		0.00000   1.00000   0.00000  -0.00000
		0.80000  -0.00000   0.60000  -0.00000
		-0.00000   0.00000   0.00000   1.00000];
			   
AB2 = [	-0.89443  -0.44721   0.00000   0.00000
		0.44721  -0.89443   0.00000   0.00000
		0.00000   0.00000   0.89443   0.44721
		0.00000   0.00000  -0.44721   0.89443]; 
  
AB1*CS*AB2;
   
s = [0.5 0.5 0.5 0.5]';		% 00, 01, 10, 11

%s1 = s;
%s = U*s
%s1 = AB2*s1
%s1.^2
%s1 = CS*s1
%s1.^2
%s1 = AB1*s1
%s = -s + 2*mean(s)

%(-0.158114)^2 + (-0.158114)^2 
%(-0.474341)^2 + (+0.474342)^2
%
%+0.664078^2+0.094868^2
%0.474341^2+0.284605^2
%+0.474342^2+0.474342^2


I = eye(2);
O = zeros(2,2);
X = [0 1; 1 0];
Z = [1 0; 0 -1];
CX = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
XC = [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0];
	
	
% Iteration 1

	AB2;
	dqAB2 = kron(AP(-pi/2),I) * kron(Rz(pi),I) * kron(I,Rz(6.2832)) * kron(I,Ry(0.92730));
	qxAB2 = kron(AP(pi/2),I) * AB2 * s;
	abs(qxAB2).^2;
	
	%ry q0,-0.927295
	%rz q0,-6.2832
	%rz q1,-3.141593

	%(+0.000005,-0.670820) |00> +
	%(-0.000002,-0.223607) |01> +
	%(-0.000005,+0.670820) |10> +
	%(+0.000002,+0.223607) |11> +
	
	CS * AB2;
	dqCSAB2 = XC * kron(Ry(-0.927295),I) * XC * kron(Ry(-0.927295),I) * dqAB2;
	qxCSAB2 = XC * kron(Ry(-0.927295),I) * XC * kron(Ry(-0.927295),I) * qxAB2;
	abs(qxCSAB2).^2;
	
	%ry q1,0.927295
	%cnot q0,q1
	%ry q1,0.927295
	%cnot q0,q1
	
	%(+0.000007,-0.939149) |00> +
	%(-0.000002,-0.223607) |01> +
	%(+0.000001,-0.134164) |10> +
	%(+0.000002,+0.223607) |11> +
	
	U;
	abs(U*s).^2;
	AB1 * CS * AB2;
	dqU =  kron(I,AP(1.570796)) * kron(I,Rz(3.141593)) * kron(I,Ry(0.927295)) * XC * kron(Rz(-1.570796),I) * XC * kron(Rz(1.570796),I) * kron(I,AP(-0.785398)) * kron(I,Rz(-1.570796)) * kron(I,Rz(-3.141593)) * dqCSAB2;
	qxU = kron(I,AP(-1.570796)) * kron(I,AP(0.785398)) * kron(I,AP(1.570796)) * kron(I,Rz(3.141593)) * kron(I,Ry(0.927295)) * XC * kron(Rz(-1.570796),I) * XC * kron(Rz(1.570796),I) * kron(I,AP(-0.785398)) * kron(I,Rz(-1.570796)) * kron(I,Rz(-3.141593)) * qxCSAB2;
	abs(qxU).^2;
	
	%rz q0,3.141593
	%rz q0,1.570796
	%rz q1,-1.570796
	%cnot q0,q1
	%rz q1,1.570796
	%cnot q0,q1
	%ry q0,-0.927295
	%rz q0,-3.141593		
	
	%(-0.664676,-0.664684) |00> +
	%(-0.155560,-0.155567) |01> +
	%(-0.155563,-0.155564) |10> +
	%(+0.098996,+0.098994) |11> +	

	s = qxU - 2*mean(qxU);
	(abs(s).^2)';
	
	%h q0
	%h q1
	%x q0
	%x q1
	%h q0
	%cnot q1,q0
	%h q0
	%x q0
	%x q1
	%h q0
	%h q1
	
	%(-0.226275,-0.226274) |00> +
	%(+0.282841,+0.282844) |01> +
	%(+0.282838,+0.282847) |10> +
	%(+0.537398,+0.537404) |11> +
	
	qxs = abs([-0.226275-0.226274i +0.282841+0.282844i +0.282838+0.282847i +0.537398+0.537404i]).^2
	
% Iteration 2
	
	%(+0.000005,-0.896800) |00> +
	%(+0.000003,-0.004000) |01> +
	%(-0.000011,-0.004000) |10> +
	%(-0.000004,+0.442399) |11> +
	
	qxs = abs([+0.000005-0.896800i +0.000003-0.004000i -0.000011-0.004000i -0.000004+0.442399i]).^2
   
% Iteration 3

	%(+0.560057,-0.560045) |00> +
	%(+0.286359,-0.286341) |01> +
	%(+0.286343,-0.286357) |10> +
	%(+0.149502,-0.149497) |11> +
	
	qxs = abs([+0.560057-0.560045i +0.286359-0.286341i +0.286343-0.286357i +0.149502-0.149497i]).^2

	
%% Test QX, Matlab compatibility

%t = -3.141593;
%Rz = [exp(1i*t/2) 0; 0 exp(-1i*t/2)]
%kron(I,Rz)*s							% rz q0,-t

%t = 0.927295;
%Ry = [cos(t/2) sin(t/2); -sin(t/2) cos(t/2)];
%kron(I,Ry)*s							% ry q0,-t

%s = kron(I,Z)*s						% z q0
%s = kron(Z,I)*s						% z q1
%[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]*s		% cnot q1,q0
%[1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]*s		% cnot q0,q1