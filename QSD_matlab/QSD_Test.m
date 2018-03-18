close all
clear all
clc

%% Generate Unitary for testing

% 	disp('Arbitrary Special Unitary Matrix to be decomposed')
% 	dim = 3;
% 	U = randUM(dim);	% Arbitrary Unitary

	p = '00';			% String to search
	q = 0.25;			% Smoothness factor (1 == Grover)
	d = size(p,2);
	N = 2^d;
	bp = ones(1,N);
	for i = 1:N
		h = pdist([sprintf('%s',dec2bin(i-1,d));p],'hamming')*d;
		bp(i) = sqrt((q^h)*((1-q)^(d-h)));
	end
	Ph = diag([1 1 1 1 1 1 1 exp(1i*pi)]);		% To make det(BO) +1 from -1
	BO = eye(N)-2*bp'*bp;
	dim = size(p,2);
	U = BO;				% Binomial Unitary for QuAM distributed query Oracle
	
%% Test Decomposition	
	
% 	det(U)
	U
	if dim >= 2
		Ud = QSD_Main(U,1);
		maxerr = max(max(Ud-U))
	else
		zyz(U)
	end

%% Generating Binomial Oracle Unitary for testing QuAM with distributed queries
% 

% 
% 	% surf(real(BO))
% 	% BOD = QSD_Main(BO);
% 	% state = ones(N,1)/sqrt(N)
% 	% hold on
% 	% plot([0:N-1],state,'-.r')
% 	% state = BOD*state
% 	% plot([0:N-1],state,'.-b')
% 	% axis([0 N-1 -1 1])