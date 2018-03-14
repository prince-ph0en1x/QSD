function [delta,alpha,theta,beta] = zyz(U)
	
   	%% Qubiter Method
	
	delta = atan2(imag(det(U)),real(det(U)))/size(U,1);
	SU = U/exp(1i*delta);
	A = SU(1,1);
    B = SU(1,2);	
	cw = real(A);
	wx = imag(B);
	wy = real(B);
	wz = imag(A);
	sw = sqrt(wx^2 + wy^2 + wz^2);
	wx = wx/sw;
	wy = wy/sw;
	wz = wz/sw;
	t1 = atan2(wz*sw,cw);
	t2 = atan2(wx,wy);
	alpha = t1 + t2;
	beta = t1 - t2;
	theta = 2*atan2(sw*sqrt(wx^2 + wy^2),sqrt(cw^2 + (wz*sw)^2));
		
	%% My calculations
% 	
% 	A = SU(1,1);
% 	B = SU(1,2);
% 	C = SU(2,1);
% 	D = SU(2,2);
% 	alpha = log(-A*B/(C*D))/2
% 	beta = log(-A*C/(B*D))/2
% 	theta = 2*atan(B*exp(beta)/A)  % +/- n*pi
% 	delta = (log(A*D)-log(cos(theta/2)^2))/(2*1i)
% 	% delta = log(det(U))/1i
	
	%% Tests
% 	
%	Phd = AP(delta);
% 	Rza = Rz(alpha);
% 	Ryt = Ry(theta);
% 	Rzb = Rz(beta);
% 	decomposedSU = Phd*Rza*Ryt*Rzb;
	
end