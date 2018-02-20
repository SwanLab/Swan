function bulk_interpretation
close all
V = 0.6;
Vini = 0.7985;
p = 1.5; 
nlambda = 10;
lambda = linspace(1,20,nlambda);
colorspec = hsv(nlambda+1);
for ilambda = 1:length(lambda)


x = linspace(0.2,1,1000); 
hx = 2*x.^(-p) + 2*x.^(-p/2);
cx = lambda(ilambda)*(x/V-1).^1;
Jx = hx + cx;

%y = 2*x.^(-p) + lambda*(x/V-1); 
 

[ymin,imin] = min(Jx);
[~,i_ini] = min(abs(x-Vini));
figure(1)
hold on
h(ilambda) = plot(x,Jx,'Color',colorspec(ilambda+1,:));
hold on;
plot(x(imin),Jx(imin),'r+','MarkerSize',10)
plot(x(i_ini),Jx(i_ini),'k+','MarkerSize',10)
hold off
axis([0 1 0 20])
pause(0.5)
legend(h,strcat(repmat(['lambda= '],ilambda,1),num2str(lambda(1:ilambda)')))
end
end