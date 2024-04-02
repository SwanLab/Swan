v  = [1,0.78,0.6,0.52];%0.41];
uy = [7.6e2,8.9e2,1.1e3,1.9e3];%1.9e3];
ux = [2.4e2,2.7e2,3.2e2,4.3e2];%5.6e2];
ct = ones(1,numel(uy));

A  = [(v.^2)',v',ct'];
b1 = uy';
b2 = ux';

theta1 = A'*A\A'*b1;
theta2 = A'*A\A'*b2;

x = 0.4:0.01:1;
f_1 = x;
f_2 = x;
for i = 1:numel(x)
    f_1(i) = theta1'*[x(i)^2;x(i);1];
    f_2(i) = theta2'*[x(i)^2;x(i);1];
end
figure
hold on
plot(v,uy,'o')
plot(x,f_1,'k')
xlabel('Volume fraction')
ylabel('Displacement Y')
grid minor

figure
hold on
plot(v,ux,'o')
plot(x,f_2,'k')
xlabel('Volume fraction')
ylabel('Displacement X')
grid minor

