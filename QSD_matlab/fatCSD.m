% Source: http://www.ar-tiste.com/m-fun/csd_qc.m

function [L0,L1,cc,ss,R0,R1] = fatCSD(U)

% This function performs a special case of fat CSD; 
% namely, the case that has been found useful in quantum computing,
% wherein the matrix U being decomposed 
% is a 2^n dimensional unitary matrix,
% and we partition U into four square matrices of the same size.
% This function calls csd() and is a trivial extension of it.
% csd() performs thin CSD. 

% U = [U00, U01] = [L0    ][ cc  ss][R0    ] 
%     [U10, U11]   [    L1][-ss  cc][    R1]
%
% Thin version of CSD (performed by csd()) gives 
% cc,ss, LO, L1 and R0, but 
% it doesn't give R1.
% This subroutine calls csd() and then calculates R1

    %ns = number of states
    %nb = number of bits
    ns = size(U,1);
    nb = 0;
    k = 1;
    while (k<ns)
        nb = nb+1;	
        k = k*2;
    end
    if (k~=ns)
        error('dimension of input matrix for csd_qc is not power of 2');
    end
    if (k==1)
        error('dimension of input matrix for csd_qc is 1');
    end

    nsh = ns/2; %ns half
    U00 = U(1:nsh, 1:nsh);
    U10 = U(nsh+1:ns, 1:nsh);

    [L0,L1,R0,cc,ss] = thinCSD(U00,U10);
    R0 = R0';
    ss = -ss;

    R1 = zeros(nsh, nsh);
    for j=1:nsh
        if abs(ss(j,j))>abs(cc(j,j))
            U01 = U(1:nsh, nsh+1:ns);
            tmp = (L0'*U01);
            R1(j,:) = tmp(j,:)/ss(j,j);
        else
            U11 = U(nsh+1:ns, nsh+1:ns);
            tmp = (L1'*U11);
            R1(j,:) = tmp(j,:)/cc(j,j);
        end
    end