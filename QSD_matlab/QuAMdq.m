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
 
 % searching for 000 when memory has states from 010 to 111.
 % quamdq successfully amplifies 010 and 100, as 000 and 001 is not there. Pmax is 0.61
 % completing query fails to amplify if state not in memory
 
 	% Oracle
 	p = '000';			% String to search
 	q = 0.25;			% Smoothness factor (1 == Grover)
 	d = size(p,2);
 	N = 2^d;
 	bp = ones(1,N);
 	for i = 1:N
% 		h = pdist([sprintf('%s',dec2bin(i-1,d));p],'hamming')*d;	% MATLAB
 		h = sum(sprintf('%s',dec2bin(i-1,d)) ~= p);					% OCTAVE
 		bp(i) = sqrt((q^h)*((1-q)^(d-h)));
 	end
 	Ph = diag([1 1 1 1 1 1 1 exp(1i*pi)]);		% To make det(BO) +1 from -1
 	BO = eye(N)-2*bp'*bp

	BOD = QSD_opql(BO,1,[0:d-1])	% Arg2 : 1 - no qasm, no AP; 2 - qasm, AP; 3+ - qasm, no AP
% 	BOD = QSD_qasm(BO,1,[0:d-1])	% Arg2 : 1 - no qasm, no AP; 2 - qasm, AP; 3+ - qasm, no AP
% 	maxerr = max(max(BOD-BO))
	maxerrabs = max(max(abs(BOD)-abs(BO)))

U1 =[
  -0.0000 + 1.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i  -0.0000 - 1.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.0000 + 1.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.0000 - 1.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.0000 + 1.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.0000 - 1.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.0000 + 1.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.0000 - 1.0000i];

U2 =[
    0.6414    0.7672         0         0         0         0         0         0
   -0.7672    0.6414         0         0         0         0         0         0
         0         0    0.6414    0.7672         0         0         0         0
         0         0   -0.7672    0.6414         0         0         0         0
         0         0         0         0    0.6414    0.7672         0         0
         0         0         0         0   -0.7672    0.6414         0         0
         0         0         0         0         0         0    0.6414    0.7672
         0         0         0         0         0         0   -0.7672    0.6414];

U3 =[
  -0.1814 + 0.4379i  -0.0737 + 0.1779i   0.0370 - 0.0893i  -0.2121 + 0.5120i   0.0370 - 0.0893i  -0.2121 + 0.5120i   0.0214 - 0.0515i  -0.1224 + 0.2956i
   0.3306 - 0.7981i   0.0334 - 0.0806i   0.0214 - 0.0515i  -0.1224 + 0.2956i   0.0214 - 0.0515i  -0.1224 + 0.2956i   0.0123 - 0.0298i  -0.0707 + 0.1707i
   0.0370 - 0.0893i  -0.2121 + 0.5120i  -0.2241 + 0.5410i   0.1712 - 0.4133i   0.0214 - 0.0515i  -0.1224 + 0.2956i   0.0123 - 0.0298i  -0.0707 + 0.1707i
   0.0214 - 0.0515i  -0.1224 + 0.2956i   0.3059 - 0.7386i   0.1747 - 0.4219i   0.0123 - 0.0298i  -0.0707 + 0.1707i   0.0071 - 0.0172i  -0.0408 + 0.0985i
   0.0370 - 0.0893i  -0.2121 + 0.5120i   0.0214 - 0.0515i  -0.1224 + 0.2956i  -0.2241 + 0.5410i   0.1712 - 0.4133i   0.0123 - 0.0298i  -0.0707 + 0.1707i
   0.0214 - 0.0515i  -0.1224 + 0.2956i   0.0123 - 0.0298i  -0.0707 + 0.1707i   0.3059 - 0.7386i   0.1747 - 0.4219i   0.0071 - 0.0172i  -0.0408 + 0.0985i
   0.0214 - 0.0515i  -0.1224 + 0.2956i   0.0123 - 0.0298i  -0.0707 + 0.1707i   0.0123 - 0.0298i  -0.0707 + 0.1707i  -0.2383 + 0.5754i   0.2528 - 0.6103i
   0.0123 - 0.0298i  -0.0707 + 0.1707i   0.0071 - 0.0172i  -0.0408 + 0.0985i   0.0071 - 0.0172i  -0.0408 + 0.0985i   0.2977 - 0.7188i   0.2219 - 0.5356i];

U3*U2*U1
	
	% Initialization
	n = N-2;		% Number of stored elements in associative memory
	r = 1;			% Number of known exact solutions possible
	sln = 2;		% Number of acceptable suboptimal solutions in memory for Psoln
	idx = 3;		% Index of an acceptable solution in memory for Psoln
	%s = [zeros(1,N-n) ones(1,n)/sqrt(n)];
	s = sqrt([0 0 1 1 1 1 1 3]/N);
	plot([0:N-1],abs(s).^2,'sr')
	hold on
	s1 = s;
	Psoln = [];
	
% 	plot([0:N-1],[ones(1,r) zeros(1,N-r)],'-.m')
	plot([0:N-1],bp,'-.g')
	
	s
	stry = s;
	s = (BOD*s')';					% Distributed Query
	s
	
	stry = (U1*stry')'
	stry = (U2*stry')'
	stry = (U3*stry')'
	
% 	s = -s + 2*mean(s);				% Diffuse
% 	s
% 	s1(1:r) = -s1(1:r);				% q = 0: Mark
% 	s1 = -s1 + 2*mean(s1);			% q = 0: Diffuse
% 	Psoln = [Psoln sln*s(idx)^2];
% 	
% 	s(N-n+1:end) = -s(N-n+1:end);	% Memorised Oracle
% 	s
% 	s = -s + 2*mean(s);				% Diffuse
% 	s
% 	s1(N-n+1:end) = -s1(N-n+1:end);	% q = 0: Mark
% 	s1 = -s1 + 2*mean(s1);			% q = 0: Diffuse
% 	Psoln = [Psoln sln*s(idx)^2];
% 	
% % 	T = 420;
% 	T = 6;
% 	
% 	for i = 1:T
% 		s = (BOD*s')';				% Distributed Query
% 		s
% 		s = -s + 2*mean(s);			% Diffuse
% 		s
% 		s1(1:r) = -s1(1:r);			% q = 0: Mark
% 		s1 = -s1 + 2*mean(s1);		% q = 0: Diffuse
% 		Psoln = [Psoln sln*s(idx)^2];
% 	end
% % 	plot([0:N-1],s1.^2,'or')
% 	plot([0:N-1],abs(s).^2,'xb')	
% % 	sln*abs(s(idx))^2
% 	axis([0 N-1 0 1])
% % 	
% 	legend('QuAM','DQ Oracle','DQ Psoln')
% % 	legend('QuAM','SQ Oracle','DQ Oracle','SQ Psoln','DQ Psoln')
% 	xlabel('State')
% 	ylabel('Probability')
% % 	
% % 	figure(2)
% % 	plot([-1:T],Psoln)
% % 	[v,i]=max(Psoln)
% % 	axis([-1 T 0 1])

%% QuAMdq iteration (tested, MATLAB+QX)
%
%	% Find how many iterations for max Psoln
%	% Continue till abs(<s|bp>) ~= 1
%
% 	bp = [4 2 2 1]/5;		% Find nearest to '00'
% 	N = size(bp,2);
% 	BO = eye(N)-2*bp'*bp
% 	BOD = QSD_qasm(BO,2,[0:1])
% 	maxerr = max(max(BOD-BO))
% 	maxerrabs = max(max(abs(BOD)-abs(BO)))
% 	plot([0:N-1],bp,'-.g')
% 	hold on
% 	s = [ones(1,N)/sqrt(N)];
% 	plot([0:N-1],abs(s).^2,'xb')
% 	T = round(2*pi/(2*asin(sum(bp)/sqrt(N))));
% 	for i = 1:T
% 		s = (BOD*s')';					% Distributed Query
% 		s = s - 2*mean(s);				% Diffuse
% 		sprob = abs(s).^2
% % 		s
% 	end
% 	plot([0:N-1],abs(s).^2,'or')
% 	axis([0 N-1 0 1])	