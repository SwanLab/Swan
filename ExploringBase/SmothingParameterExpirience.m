x = linspace(0,0.925,1000);
alpha = 6;
beta = 20;
q = 4*(1./(1-x.^beta)).^alpha;
figure(1)
hold on
plot(x,q)
