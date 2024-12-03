close all
clear all
m=10;
K=2*diag(ones(2*m+1,1),0)  -1*diag(ones(2*m,1),1) -1*diag(ones(2*m,1),-1);
M=4*diag(ones(2*m+1,1),0) + 1*diag(ones(2*m,1),1) +1*diag(ones(2*m,1),-1);
M(1,1)=2;M(end,end)=2;
K(1,1)=1;K(end,end)=1;
Kred=K(2:end-1,2:end-1);
Mred=M(2:end-1,2:end-1);
% [v,d]=eigs(Kred\Mred);
% plot(v(:,1))
F=zeros(m*2+1,1);
F(10)=10000;
Fred=F(2:end-1);

%Jacobi

err=1;
x_old=zeros(m*2+1-2,1);
D=diag(diag(Kred));
LU=Kred-D;
while err>1e-6
x_new=D\(Fred-LU*x_old);
 err=norm(x_new-x_old);
 x_old=x_new;
end

%direct
x=Kred\Fred;

error=x_new-x;
