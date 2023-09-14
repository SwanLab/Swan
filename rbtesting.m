clear all
close all
file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

ss.mesh=mesh;
ss.fvalues=[0 0 1];
ss.refPoint= [0.02 0.02];

a=RigidBodyFunction(ss);
xV=zeros(2,1,16);
xV(1,:,:)=0.33;
xV(2,:,:)=0.33;
a.evaluate(xV);
P1FUNC=a.project('P1');
P1FUNC.plot
P1FUNC.print('GiD')

