function rbtesting
file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

%%% mesh natural

ss.mesh=mesh;
ss.fvalues=[0 10 0];
ss.refPoint= [0.02 0.02];

a=RigidBodyFunction(ss);
xV=zeros(2,1,16);
xV(1,:,:)=0.33;
xV(2,:,:)=0.33;
a.evaluate(xV);
p1FUNC = a.project('P1');
p1FUNC.plot
p1FUNC.print('GiD')

% FEM
%solution plot

% eigmodes LHS

% ModalFunction build


% Plot an example with given values


end
