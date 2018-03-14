close all
clear all
clc

%% Generate Arbitrary Unitary for testing

disp('Arbitrary Special Unitary Matrix to be decomposed')
dim = 4;
U = randUM(dim);
if dim >= 2
	Ud = QSD_Main(U);
	max(max(Ud-U))
else
	zyz(U)
end