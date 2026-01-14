function u=GepOneSide(A,B,c)
[m1,~]=size(A);
[m2,~]=size(B);
e1=ones(m1,1);
e2=ones(m2,1);
con=m1-c*m2;
if con~=0
    tmp=e1'*A-c*e2'*B;
    H=A'*A-c*B'*B-1/con*tmp'*tmp;
    [d,v]=eig(H);
    [tmp2,sign]=min(diag(v));
    w=d(:,sign);
    b=-tmp/con*w;
else
    H=A'*A-c*B'*B;
    [d,v]=eig(H);
    [tmp,sign]=min(diag(v));
    w=d(:,sign);
    b=1/m1*sum(A,1)*w;
end
u=[w;b]';
end