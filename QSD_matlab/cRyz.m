function decomposed = cRyz(t1,t2)
% Apply Ry(t1) if control bit 0, and Ry(t2) if control bit 1

    CX = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
    rt1 = (t1+t2)/2;
    rt2 = (t1-t2)/2;
    r1q = [cos(rt1/2) sin(rt1/2); -sin(rt1/2) cos(rt1/2)];
    r2q = [cos(rt2/2) sin(rt2/2); -sin(rt2/2) cos(rt2/2)];
    decomposed = CX * kron(eye(2),r2q) * CX * kron(eye(2),r1q);   
    
end