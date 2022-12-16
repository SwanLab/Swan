%% Swan vs Feppon

close all;
clear;
clc;

x = -20:0.25:20;
y = 0:0.01:1;

X = zeros(length(y),length(x));
for i = 1:length(y)
    X(i,:) = x;
end

Y = zeros(length(y),length(x));
for i = 1:length(x)
    Y(:,i) = y;
end

figure
f1 = min(max(-10,X),10);
s1 = surf(X,Y,f1);

figure
f2 = min((X.*Y)./(max(1e-9,Y)),0.9*X.*Y);
s2 = surf(X,Y,f2);