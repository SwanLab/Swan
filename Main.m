function Main

%% Preprocess
nunkn = 2;

mesh = Mesh();
bc = BC(nunkn);

%% Process

geometry = Geometry(mesh);
dof = DOF(geometry.nnode,mesh.connec,nunkn,mesh.npnod,bc.displacements);

element = Element_Elastic();
element = element.computeLHS(nunkn,mesh.nelem,geometry);
element = element.computeRHS(nunkn,mesh.nelem,geometry.nnode,bc,dof.idx);

[LHS,RHS] = Assemble.Compute(element,geometry.nnode,nunkn,dof);

d_u = zeros(dof.ndof,1);
d_u = Solver.analytical(d_u,LHS,RHS',dof.vR,dof.vL,bc.displacements);

%% Postprocess
disp(d_u);

end
