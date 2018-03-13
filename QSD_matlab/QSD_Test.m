close all
clear all
clc

%% Generate Arbitrary Unitary for testing

disp('Arbitrary Special Unitary Matrix to be decomposed')
dim = 3;
SU = randUM(dim)	% SU matrix for now
QSD_Main(SU);