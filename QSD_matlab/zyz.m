function [Phd,Rza,Ryt,Rzb] = zyz(U)
%function [delta,alpha,theta,beta] = zyz(U)
	
    A = U(1,1);
    B = U(1,2);
    C = U(2,1);
    D = U(2,2);
    
    alpha = log(-A*B/(C*D))/2;
    beta = log(-A*C/(B*D))/2;
    theta = 2*atan(B*exp(beta)/A);  % +/- n*pi
    delta = (log(A*D)-log(cos(theta/2)^2))/(2*1i);
    
    Phd = [exp(1i*delta) 0; 0 exp(1i*delta)];
    Rza = [exp(1i*alpha/2) 0; 0 exp(-1i*alpha/2)];
    Ryt = [cos(theta/2) sin(theta/2); -sin(theta/2) cos(theta/2)];
    Rzb = [exp(1i*beta/2) 0; 0 exp(-1i*beta/2)];
    
    % decomposed = Rzb*Ryt*Rza
    % U == decomposed 
    
end