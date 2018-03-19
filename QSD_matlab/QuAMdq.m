close all
clear all
clc

%% Generating Binomial Oracle Unitary for testing QuAM with distributed queries
% 
% 	p = '000';		% String to search
% 	q = 0.25;		% Smoothness factor (1 == Grover)
% 	d = size(p,2);
% 	N = 2^d;
% 	bp = ones(1,N);
% 	for i = 1:N
% 		h = pdist([sprintf('%s',dec2bin(i-1,d));p],'hamming')*d;
% 		bp(i) = sqrt((q^h)*((1-q)^(d-h)));
% 	end
% 	Ph = diag([1 1 1 1 1 1 1 exp(1i*pi)]);		% To make det(BO) +1 from -1
% 	BO = eye(N)-2*bp'*bp
% 
% 	% surf(real(BO))
% 	% BOD = QSD_Main(BO);
% 	% state = ones(N,1)/sqrt(N)
% 	% hold on
% 	% plot([0:N-1],state,'-.r')
% 	% state = BOD*state
% 	% plot([0:N-1],state,'.-b')
% 	% axis([0 N-1 -1 1])

%% QuAM Iteration 1 solution Modified Loop
% 
% 	% Initialization
% 	n = N;		% Number of stored elements in associative memory
% 	r = 1;		% Number of known solutions among stored elements
% 	s = [ones(1,n)/sqrt(n) zeros(1,N-n)];
% 	plot([0:N-1],s,'ob')
% 	plot([0:N-1],s.^2,'-.k')
% 	hold on
% 	s1 = s;
% 	T = 2;
% 	for i = 1:T
% 		s = (BO*s')';			% Distributed Query
% 		s = -s + 2*mean(s);		% Diffuse
% 		s1(1:r) = -s1(1:r);		% q = 0: Mark
% 		s1 = -s1 + 2*mean(s1);	% q = 0: Diffuse
% 	end
% 	plot([0:N-1],s1.^2,'or')	% Good for T = 2
% 	plot([0:N-1],s.^2,'xk')		% Good for T = 3
% 	axis([0 N-1 0 1])
	
%% QuAM Iteration 2 solution Modified Loop
% 
% % searching for 000 when memory has states from 010 to 111.
% % quamdq successfully amplifies 010 and 100, as 000 and 001 is not there. Pmax is 0.61
% % completing query fails to amplify if state not in memory
% 
% 	% Oracle
% 	p = '000';		% String to search
% 	q = 0.25;		% Smoothness factor (1 == Grover)
% 	d = size(p,2);
% 	N = 2^d;
% 	bp = ones(1,N);
% 	for i = 1:N
% 		h = pdist([sprintf('%s',dec2bin(i-1,d));p],'hamming')*d;
% 		bp(i) = sqrt((q^h)*((1-q)^(d-h)));
% 	end
% 	Ph = diag([1 1 1 1 1 1 1 exp(1i*pi)]);		% To make det(BO) +1 from -1
% 	BO = eye(N)-2*bp'*bp;
% 
% 	% Initialization
% 	n = N-2;		% Number of stored elements in associative memory
% 	r = 1;			% Number of known exact solutions possible
% 	sln = 2;		% Number of acceptable suboptimal solutions in memory for Psoln
% 	idx = 3;		% Index of an acceptable solution in memory for Psoln
% 	s = [zeros(1,N-n) ones(1,n)/sqrt(n)];
% 	plot([0:N-1],s.^2,'-k')
% 	hold on
% 	s1 = s;
% 	Psoln = [];
% 	
% 	plot([0:N-1],[ones(1,r) zeros(1,N-r)],'-.m')
% 	plot([0:N-1],bp,'-.c')
% 	
% 	s = (BO*s')';					% Distributed Query
% 	s = -s + 2*mean(s);				% Diffuse
% 	s1(1:r) = -s1(1:r);				% q = 0: Mark
% 	s1 = -s1 + 2*mean(s1);			% q = 0: Diffuse
% 	Psoln = [Psoln sln*s(idx)^2];
% 	
% 	s(N-n+1:end) = -s(N-n+1:end);	% Memorised Oracle
% 	s = -s + 2*mean(s);				% Diffuse
% 	s1(N-n+1:end) = -s1(N-n+1:end);	% q = 0: Mark
% 	s1 = -s1 + 2*mean(s1);			% q = 0: Diffuse
% 	Psoln = [Psoln sln*s(idx)^2];
% 	
% 	T = 420;
% 	for i = 1:T
% 		s = (BO*s')';				% Distributed Query
% 		s = -s + 2*mean(s);			% Diffuse
% 		s1(1:r) = -s1(1:r);			% q = 0: Mark
% 		s1 = -s1 + 2*mean(s1);		% q = 0: Diffuse
% 		Psoln = [Psoln sln*s(idx)^2];
% 	end
% 	plot([0:N-1],s1.^2,'or')
% 	plot([0:N-1],s.^2,'xb')	
% 	sln*s(idx)^2
% 	axis([0 N-1 0 1])
% 	
% 	legend('QuAM','SQ Oracle','DQ Oracle','SQ Psoln','DQ Psoln')
% 	xlabel('State')
% 	ylabel('Probability')
% 	
% % 	figure(2)
% % 	plot([-1:T],Psoln)
% % 	[v,i]=max(Psoln)
% % 	axis([-1 T 0 1])

%% QuAMdq iteration

	% Find how many iterations for max Psoln
	% Continue till abs(<s|bp>) ~= 1

	bp = [4 2 2 1]/5;		% Find nearest to '00'
	N = size(bp,2);
	BO = eye(N)-2*bp'*bp
%	BOD = QSD_Main(BO,1);
	BOD = QSD_qasm(BO,1,[0:1])
	maxerr = max(max(BOD-BO))
	plot([0:N-1],bp,'-.g')
	hold on

	s = [ones(1,N)/sqrt(N)];
	plot([0:N-1],s.^2,'xb')

	% beta0 = s*bp';
	% alpha0 = s(1);
	T = round(2*pi/(2*asin(sum(bp)/sqrt(N))));
	for i = 1:T
		s = (BO*s')';					% Distributed Query
		s = s - 2*mean(s);				% Diffuse
		s.^2;
	end
	
	plot([0:N-1],s.^2,'or')
	axis([0 N-1 0 1])	