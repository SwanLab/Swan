function DATAIN = ParamCurv_DIFFSLICES_3D(DATAINPUT)
% Function for definining the coordinates of a structure
% formed by concatenation of different  meshes  (DATAINPUT.NameFileMeshLOC)
% along a curve defined by the mesh DATAINPUT.NameMidLineMesh
% 3D version, copy of ParametricCurvature_3D.m
% JAHO, started 1-Nov-2020 (Sofia Airport, BLG)
% ---------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
% 1. Read midline mesh
DATAINPUT = DefaultField(DATAINPUT,'ReadMaterialColumn',1) ; 
[DATA1Dmsh]= ReadMeshFileStr(DATAINPUT.NameMidLineMesh,'READ_MATERIAL_COLUMN',DATAINPUT.ReadMaterialColumn)  ;
% 2. Determine initial point

if size(DATA1Dmsh.COOR,2) == 2
    DATA1Dmsh.COOR = [DATA1Dmsh.COOR , zeros(size(DATA1Dmsh.COOR,1),1)] ; 
end

endPOINTS1= setdiff(DATA1Dmsh.CN(:,1),DATA1Dmsh.CN(:,2) ) ;
endPOINTS2= setdiff(DATA1Dmsh.CN(:,2),DATA1Dmsh.CN(:,1) ) ;
endPOINTS = unique([endPOINTS1,endPOINTS2]) ; 
CandidateIniPoint =  endPOINTS(1) ; 
DATAINPUT = DefaultField(DATAINPUT,'InitialPoint',CandidateIniPoint) ;
if isempty(intersect(endPOINTS,DATAINPUT.InitialPoint)) 
    InitialPoint = CandidateIniPoint ; 
else
    InitialPoint = DATAINPUT.InitialPoint ;    
end

% 3. Sequence of points (coordinates)
nelem  = size(DATA1Dmsh.CN,1);
SEQ_POINTS = zeros(size(DATA1Dmsh.COOR,1),1) ;
SEQ_POINTS(1) = InitialPoint ;
SEQ_elem = zeros(size(DATA1Dmsh.CN,1),1) ;
for ielem = 1:nelem
    [indelem,innode ]= find(SEQ_POINTS(ielem) == DATA1Dmsh.CN)  ;
    if length(indelem ) > 1
        ISEQUAL =   ismember(indelem,SEQ_elem(ielem-1)) ;
        indselect = find(ISEQUAL == 0) ;
        indelem= indelem(indselect) ;
        innode= innode(indselect) ;
    end
    jnode= setdiff([1,2],innode) ;
    SEQ_POINTS(ielem+1) = DATA1Dmsh.CN(indelem,jnode) ;
    SEQ_elem(ielem) = indelem ;
end
COOR1D = DATA1Dmsh.COOR(SEQ_POINTS,:)  ;  % Coordinates 1D points
%MaterialType =; 
% if isempty()
% end

if  isempty(DATA1Dmsh.MaterialType)
    DATA1Dmsh.MaterialType = ones(size(DATA1Dmsh.CN,1),1) ; 
end 

TYPEslices =  DATA1Dmsh.MaterialType(SEQ_elem) ;
DATAIN.TYPEslices = TYPEslices ; 
% -------------------------------------------------------------------TYPEslices----TYPEslices
% 4. SPLINE INTERPOLATION
% https://es.mathworks.com/help/curvefit/cscvn.html
SPLcurve = cscvn(COOR1D') ;

% Maximum value of the curve parameter 
%diffCOOR = diff(COOR1D,1) ; 
%sKNOTS =  cumsum(sqrt(sum(diffCOOR.^2,2))) ;

sKNOTS = cumsum([0;((diff(COOR1D).^2)*ones(3,1)).^(1/4)]).';  % Breaking sequence, from Matlab 

 


DATAINPUT = DefaultField(DATAINPUT,'TORSION',[]) ; 
DATAINPUT.TORSION = DefaultField(DATAINPUT.TORSION,'TOTAL_REVOLUTIONS',0) ; 
DATAINPUT.TORSION = DefaultField(DATAINPUT.TORSION,'TYPE','LINEAR') ;
ANG_TORS = [] ;
switch  DATAINPUT.TORSION.TYPE
    case 'LINEAR'
        ANG_TORS  = sKNOTS./sKNOTS(end)*2*pi*DATAINPUT.TORSION.TOTAL_REVOLUTIONS ; 
        
        
    otherwise 
        error('Option not implemented')
     
end



figure(1)
hold on 
fnplt(SPLcurve); hold on,
plot3(COOR1D(:,1),COOR1D(:,2),COOR1D(:,3),'o') 
% What is the output  SPLcurve. Is given in parametric form ?
%Parametric variational, or natural, cubic spline curve (in ppform) passing
%through the given sequence points (:j), j = 1:end. The parameter value t(j)
% for the j-th point follows the Eugene Lee's [1] centripetal scheme, as accumulated square root of chord length:
% And what about the derivative.



DER_SPLcurve=  fnder(SPLcurve) ;  % Derivative
f = @(s) (fnval(SPLcurve,s)) ;
df =  @(s) (fnval(DER_SPLcurve,s)) ;

% CHECKING THAT THE KNOTS ARE CORRECTLY PLACED 
% --------------------------------------------
for ik = 1:length(sKNOTS)
    xKNOT = f(sKNOTS(ik)) ; 
    plot3(xKNOT(1),xKNOT(2),xKNOT(3),'g*')
    text(xKNOT(1),xKNOT(2),xKNOT(3),num2str(ik))

end


% -----------------
% 5. READ REFERENCE  MESH of the slide
% ------------------
DATA3D = cell(size(DATAINPUT.NameFileMeshLOC)) ;
COORref = cell(size(DATAINPUT.NameFileMeshLOC)) ;
    DATAINPUT =  DefaultField(DATAINPUT,'TOL_face_identification_NODES',0.1) ;  
DATAINgeom.TOL_face_identification_NODES = DATAINPUT.TOL_face_identification_NODES ; 
for islice = 1:length(DATAINPUT.NameFileMeshLOC)
    SLICE.NAME = DATAINPUT.NameFileMeshLOC{islice};
    DATA3D{islice} = GeometrySlice(SLICE,DATAINgeom) ; %  We read the slice mesh here
    
    % Coordinates with respect to the centroid of face 1
    COORref{islice} =  DATA3D{islice}.COOR ;
    for idim = 1:length(DATA3D{islice}.CENTRf1_real)
        COORref{islice}(:,idim) = DATA3D{islice}.COOR(:,idim) -DATA3D{islice}.CENTRf1_real(idim) ;
    end
    
end
ymax = max(DATA3D{1}.COOR(:,2)) ;  ymin = min(DATA3D{1}.COOR(:,2)) ;
h.y = 0.5*(ymax-ymin);
zmax = max(DATA3D{1}.COOR(:,3)) ;  zmin = min(DATA3D{1}.COOR(:,3)) ;
h.z = 0.5*(zmax-zmin);
xmax = max(DATA3D{1}.COOR(:,1)) ;  xmin = min(DATA3D{1}.COOR(:,2)) ;
L = xmax-xmin;



%6. Dummy reference mesh
nx = 3 ;  ny = 3 ; nz = 3 ; 
Xref = linspace(0,L,nx) ;      Yref = linspace(-h.y,h.y,ny) ;  Zref = linspace(-h.z,h.z,ny) ;
[X,Y,Z] = meshgrid(Xref,Yref,Zref) ;
COORrefDUMMY = [X(:),Y(:),Z(:)] ;


 
DATAINPUT = DefaultField(DATAINPUT,'PROP',[]) ; 

DATAINloc.TYPEslices = TYPEslices ; 
DATAINPUT = DefaultField(DATAINPUT,'VARIABLE_CROSS_SECTION',[]) ; 
DATAINloc.VARIABLE_CROSS_SECTION = DATAINPUT.VARIABLE_CROSS_SECTION ; 

[  COORtransf,CentrF2,dP_1_GLO,dANG,R_G1_GLO,R_12_GLO]= DivideSpline3D(COORref,f,df,L,h,COORrefDUMMY,...
    sKNOTS,DATAINPUT.PROP,ANG_TORS,DATAINloc) ;

%
% %
% DATAIN.angDOM =-dANG;  % Now this variable is not a scalar, but a vector
% DATAIN.dX = dX;   % Length in the x direction, local
% DATAIN.dY = dY ;
DATAIN.dP_1 =  dP_1_GLO ; 
DATAIN.COORtransf = COORtransf ;
DATAIN.angDOM = dANG ;

DATAIN.ANGiniGLO = R_G1_GLO ;
DATAIN.ANGrelGLO = R_12_GLO ; 



%
% % DATAIN.CIRCLE.RADIUS = LENGTH_X/2/sin(DATAIN.angDOM/2) ;  % Radius
% % DATAIN.CIRCLE.CENTER_CIR= [0,0] ;
% % disp('**************************************************')
% % disp(['if NUMBER OF CIRCULAR DIVISIONS  = ',num2str(NUMBER_DIVISIONS_CIRCLE),'...']) ;
% % disp(['then RADIUS CYLINDER = ',num2str(DATAIN.CIRCLE.RADIUS)]) ;
% % disp('**************************************************')
% DATAIN.ISTWIST_ANGLE = 0;
% %% HELICOID data. Specify the height of the structure after completing one revolution
% % For instance, if h is the height of the cross section, then, upon
% % completing a revolution, H should be H > h/2
% HEIGHT = [] ;
% DATAIN.ELEVATION_Z = [] ; % HEIGHT/NUMBER_DIVISIONS_CIRCLE ;