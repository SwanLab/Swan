x = sym('x');
y = sym('x');

theta(:,1) = [x;y];
theta(:,2) = [x-1;y];
theta(:,3) = [x;y-1];

M  = theta'*theta;
Mv = int(int(M,'y',0,1-x),0,1)