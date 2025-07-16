clear all;

% Create Mesh
x1      = linspace(0, 1, 50);
x2      = linspace(0, 1, 50);
[xv,yv] = meshgrid(x1,x2);
[F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
s.coord  = V(:,1:2);
s.connec = F;
mesh = Mesh.create(s);

% f(x,y) = x
linear = LagrangianFunction.create(mesh,1,'P1');
dofs = linear.getDofCoord;
linear.setFValues(dofs(:,1)); % as Lagrangian Function
linear.plot();

s.operation = @(xV) linear.evaluate(xV);
s.mesh      = mesh;
linearDF = DomainFunction(s); % as Domain Function
linearDF.plot();

% g(x,y) = x^2
quadratic = LagrangianFunction.create(mesh,1,'P1');
dofs = quadratic.getDofCoord;
quadratic.setFValues(dofs(:,1).^2);
quadratic.plot();

s.operation = @(xV) quadratic.evaluate(xV)
s.mesh      = mesh;
quadraticDF = DomainFunction(s); % as Domain Function
quadraticDF.plot();

% Product of Two Scalar Lagrangian Functions
% h(x,y) = f(x,y)g(x,y)
cubic = linear.*quadratic       % Scalar Domain Function
cubic.plot();   

% Product of Two Scalar Lagrangian Functions
% h(x,y) = f(x,y)g(x,y)
cubicDF = linearDF.*quadraticDF; % Scalar Domain Function
cubicDF.plot(); 

% Scalar Product of Two Vector Domain Functions
dot = DP(Grad(linear), Grad(linear));
dot.plot();                     % Scalar Domain Function
 
% Scalar Product of a Scalar timea a Vector Domain Function
product = linearDF.*Grad(linear);
% product.plot();               % Should be a Vector Domain Function. Does not work.   
