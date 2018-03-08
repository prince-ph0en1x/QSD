% Source: http://www.ar-tiste.com/m-fun/csd.m

function [u1,u2,v,c,s]=thinCSD(q1,q2)

% Given Q1 and Q2 such that Q1'* Q1 + Q2'* Q2 = I, the
% C-S Decomposition is a joint factorization of the form
%       Q1 = U1*C*V' and Q2=U2*S*V'
% where U1,U2,V are orthogonal matrices and C and S are diagonal
% matrices (not necessarily square) satisfying
%          C'* C + S'* S = I
% The diagonal entries of C and S are nonnegative and the
% diagonal elements of C are in nondecreasing order.
% The matrix Q1 cannot have more columns than rows.
% ( Submitted by S. J. Leon )

    [m,n]=size(q1);
    [p,n]=size(q2);
    [u1,c,v]=svd(q1);
    z=eye(n);z=hankel(z(:,n));
    c(1:n,:)=z*c(1:n,:)*z;u1(:,1:n)=u1(:,1:n)*z;v=v*z;
    q2=q2*v;
    k=1;
    for j=2:n
       if c(j,j)<=1/sqrt(2)
          k=j;
       end
    end
    b=q2(:,1:k);
    [u2,r]=qr(b);
    s=u2'*q2;
    t=min(p,n);tt=min(m,p);
    if k<t
       r2=s(k+1:p,k+1:t);
       [ut,ss,vt]=svd(r2);
       s(k+1:p,k+1:t)=ss;
       c(:,k+1:t)=c(:,k+1:t)*vt;
       u2(:,k+1:p)=u2(:,k+1:p)*ut;
       v(:,k+1:t)=v(:,k+1:t)*vt;
       w=c(k+1:tt,k+1:t);
       [z,r]=qr(w);
       c(k+1:tt,k+1:t)=r;
       u1(:,k+1:tt)=u1(:,k+1:tt)*z;
    end
    for j=1:n
       if c(j,j)<0
          c(j,j)=-c(j,j);
          u1(:,j)=-u1(:,j);
       end
    end
    for j=1:t
       if s(j,j)<0
          s(j,j)=-s(j,j);
          u2(:,j)=-u2(:,j);
       end
    end