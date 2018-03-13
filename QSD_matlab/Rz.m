function U = Rz(t)
% Rotation about Z-axis
    U = [exp(1i*t/2) 0; 0 exp(-1i*t/2)];	% TBD: Check sign of matrix for convention
end