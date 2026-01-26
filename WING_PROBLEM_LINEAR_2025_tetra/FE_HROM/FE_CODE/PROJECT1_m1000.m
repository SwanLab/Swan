% Inputs example assigment 2 
% ----------------------------


%%%%%%%%%%%%%%%%%
 % ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB
NameFileMesh = 'malla2.msh'; % Name of the file containing the mesh information (Generated with GID)
[COOR,CN,TypeElement,CONNECTb,TypeElementB]=...
    ReadMeshFile(NameFileMesh)  ;
nnode = size(COOR,1) ;% Number of nodes 

% 2. MATERIAL PROPERTIES: output ConductMglo  
%-----------------------
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements
ConductMglo = zeros(ndim,ndim,nelem) ; 
% Conductivity matrix (isotropic)
kappa =  5  ; %  W/ºC
ConductM = kappa*eye(ndim) ; % eye = IDENTITY MATRIx
for e=1:nelem 
    ConductMglo(:,:,e) = ConductM ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of nodes at which temperature is prescribed
%   Copy here the list of nodes from GID
rnod =[1         0.0 0.0 0.0
2         0.0 0.050000000000000044 0.0
3         0.0 0.099999999999999978 0.0
4         0.10000000000000001 0.025000000000000001 0.0
6         0.0 0.15000000000000002 0.0
9         0.0 0.19999999999999996 0.0
10         0.20000000000000001 0.050000000000000003 0.0
14         0.0 0.25 0.0
17         0.0 0.30000000000000004 0.0
18         0.29999999999999999 0.074999999999999997 0.0
23         0.0 0.34999999999999998 0.0
29         0.0 0.40000000000000002 0.0
30         0.40000000000000002 0.10000000000000001 0.0
36         0.0 0.44999999999999996 0.0
42         0.0 0.5 0.0
44         0.5 0.125 0.0
52         0.0 0.55000000000000004 0.0
59         0.0 0.59999999999999998 0.0
65         0.59999999999999998 0.14999999999999999 0.0
71         0.0 0.65000000000000002 0.0
81         0.0 0.69999999999999996 0.0
86         0.69999999999999996 0.17499999999999999 0.0
92         0.0 0.75 0.0
105         0.0 0.80000000000000004 0.0
111         0.80000000000000004 0.20000000000000001 0.0
118         0.0 0.84999999999999998 0.0
132         0.0 0.90000000000000002 0.0
141         0.90000000000000002 0.22500000000000001 0.0
147         0.0 0.94999999999999996 0.0
164         0.0 1.0 0.0
173         1.0 0.25 0.0
202         1.1000000000000001 0.27500000000000002 0.0
227         1.2 0.29999999999999999 0.0
251         1.3 0.32500000000000001 0.0
275         1.3999999999999999 0.34999999999999998 0.0
299         1.5 0.375 0.0
321         1.6000000000000001 0.40000000000000002 0.0
344         1.7 0.42499999999999999 0.0
367         1.8 0.45000000000000001 0.0
389         1.8999999999999999 0.47499999999999998 0.0
413         2.0 0.5 0.0];  %

rnod =rnod(:,1) ;  
% Prescribed temperature 
temp_DAB = 0 ; % Degrees Celsius 
% Vector of prescribed temperatures
dR = temp_DAB*ones(size(rnod)) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: qFLUXglo, CNb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of boundary nodes whose boundary elements have prescribed heat flux
% (different from zero !)
% Edge DC 
% Copy here the list of nodes from GID
NODESb = [164         0.0 1.0 0.0
166         0.10000000000000009 1.0 0.0
168         0.19999999999999996 1.0 0.0
177         0.30000000000000004 1.0 0.0
187         0.39999999999999991 1.0 0.0
196         0.5 1.0 0.0
211         0.60000000000000009 1.0 0.0
223         0.69999999999999996 1.0 0.0
238         0.80000000000000004 1.0 0.0
252         0.89999999999999991 1.0 0.0
268         1.0 1.0 0.0
286         1.1000000000000001 1.0 0.0
303         1.2 1.0 0.0
320         1.3 1.0 0.0
338         1.3999999999999999 1.0 0.0
358         1.5 1.0 0.0
376         1.6000000000000001 1.0 0.0
394         1.7 1.0 0.0
412         1.8 1.0 0.0
432         1.8999999999999999 1.0 0.0
441         2.0 1.0 0.0] ; 
%
NODESb = NODESb(:,1) ; 
% prescribed flux (constant)
qBAR = 20 ; % W/m
% Initialization  of the global vector of nodal heat flux 
qFLUXglo = zeros(nnode,1) ; 
% Assign prescribed values
qFLUXglo(NODESb) = qBAR ; 
%%%%
% Connectivity matrix for boundary elements (Neumann boundary)
CNb = ElemBnd(CONNECTb,NODESb) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Heat source
%%
f = 6 ; %   W/m2
% Global vector of heat sources (constant)
fNOD = f*ones(nnode,1) ; 
