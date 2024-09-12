function [tXi,volume] = integ_exact(t,p,psi)
% tXi = pdeintrp(p,t,(psi<0) + gamma*(psi>=0)); 
dim.nelem=size(t,2);
dim.nndof=2*size(p,2);
dim.nnode=3;
dim.ndime=2;
dim.npnod=size(p,2);
dim.nunkn=2;% 2 en el codigo de Alex??
dim.nstre=3;% 3 en el codigo de Alex??

%element
element.conectivities=t';
element.boundary_elements=[];
element.type='TRIANGLE';
element.ngaus=1;

problembsc.problemtype='2D';

coordinates=p';

[volume,tXi] = cal_omega(psi,dim,element,problembsc,coordinates);

