close all
clear all
clc

%% Generate Arbitrary Unitary for testing

disp('Arbitrary Special Unitary Matrix to be decomposed')
dim = 3;
% SU = randUM(dim)	% SU matrix for now
% QSD_Main(SU);

%% Generating Binomial Oracle Unitary for testing QuAM

p = '011';		% String to search
q = 0.25;		% Smoothness factor (1 == Grover)
d = size(p,2);
bp = ones(1,2^d);
for i = 1:2^d
	h = pdist([sprintf('%s',dec2bin(i-1,d));p],'hamming')*d;
	bp(i) = sqrt((q^h)*((1-q)^(d-h)));
end
Ph = diag([1 1 1 1 1 1 1 exp(1i*pi)]);		% To make det(BO) +1 from -1
BO = Ph *(eye(2^d)-2*bp'*bp)

% surf(real(BO))
BOD = QSD_Main(BO)
% state = ones(2^d,1)/sqrt(2^d)
% hold on
% plot([0:2^d-1],state,'-.r')
% state = BOD*state
% plot([0:2^d-1],state,'.-b')
% axis([0 2^d-1 -1 1])