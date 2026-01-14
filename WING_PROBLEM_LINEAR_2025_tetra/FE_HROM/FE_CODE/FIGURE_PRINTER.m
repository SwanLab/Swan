%{
\begin{lstlisting}[basicstyle=\tiny]
%}
function FIGURE_PRINTER(Ax,d,mden,prisvol,TYPED_EL,NAMEPROJ)
clc;

% HYPOTHESES
%  ...
% square cells Ax = Ay

% DEFINITIONS
%
% Ax               [cube lenght] [mm]
% d                [fibre diameter] [mm]  
% mden             [mesh size, is used to establish the mesh density in the prisma face]
% prisvol          [establish the number of prisma volumes] 
% TYPED_EL         [Choose type of element, 1,2 need the inclusion diameter, 3,4,5

%FILE CREATION
fileID=fopen('rep_vol.bch','w');

if TYPED_EL==1
   PRINT_FIGURE1(Ax,d,mden,prisvol,fileID);
elseif TYPED_EL==2
   PRINT_FIGURE2(Ax,d,mden,prisvol,fileID); 
elseif TYPED_EL==3
   PRINT_FIGURE3(Ax,mden,prisvol,fileID);       
elseif TYPED_EL==4
   PRINT_FIGURE4(Ax,mden,prisvol,fileID);   
elseif TYPED_EL==5
   PRINT_FIGURE5(Ax,mden,prisvol,fileID);     
end

% CREATE PROJECT

fprintf(fileID,'Mescape Meshing CreateBoundary\n');
fprintf(fileID,'Yes\n');
fprintf(fileID,'Mescape Files WriteMesh\n');
 
SLH = '\' ;
if isunix
    SLH = '/' ; 
end

nameFILE = [pwd,SLH,NAMEPROJ,'.msh'] ; 

fprintf(fileID,'%s\n',nameFILE);


fclose(fileID);

display('END_OF_SCRIPT');

end

function PRINT_FIGURE1(Ax,d,mden,prisvol,fileID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AF = 0.25*pi*d^2;% [m^2]
AC = Ax*Ax;    % [m^2]
INC_PER = AF/AC;

Ay=Ax;

disp(Ax);
disp(d);
disp(INC_PER);

if(d>=Ax)
  display('WARNING inclusion bigger than representative volume dimmensions')  
 return;
else
      
 %CREATE GEOMETRY

 fprintf(fileID,'*****FIGURE1_CREATION\n');
 fprintf(fileID,'Mescape Geometry Create Object rectangle\n');
 fprintf(fileID,'%d %d\n',[-0.5*Ax,-0.5*Ay]);
 fprintf(fileID,'%d %d\n',[0.5*Ax,0.5*Ay]);
 fprintf(fileID,'Mescape\n');

 fprintf(fileID,'Mescape Geometry Create Object circle\n');
 fprintf(fileID,'%d %d\n',[0,0]);
 fprintf(fileID,'%2.1f %2.1f %2.1f\n',[0.0,0.0,1.0]);
 fprintf(fileID,'%d\n',0.5*d);
 fprintf(fileID,'Mescape\n');

 fprintf(fileID,'Mescape Geometry Create IntSolid2d substract\n');
 fprintf(fileID,'%d\n',1);
 fprintf(fileID,'escape\n');
 fprintf(fileID,'%d\n',2);
 fprintf(fileID,'Mescape\n');

 fprintf(fileID,'Mescape Geometry Create NurbsSurface\n');
 fprintf(fileID,'%d\n',5);
 fprintf(fileID,'escape\n');
 fprintf(fileID,'Mescape\n');

 fprintf(fileID,'Mescape Utilities Copy Surfaces DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,%5.4f\n',Ax);
 fprintf(fileID,'%d %d\n',[3,4]);
 fprintf(fileID,'escape\n');
 fprintf(fileID,'Mescape\n');

 %LAYERS DEFINITION

 fprintf(fileID,'''Layers New Fibre escape\n');
 fprintf(fileID,'''Layers Entities Fibre LowerEntities Volumes\n');
 fprintf(fileID,'%d\n',2);
 fprintf(fileID,'escape\n');
 fprintf(fileID,'Mescape\n');

 fprintf(fileID,'''Layers New Matrix escape\n');
 fprintf(fileID,'''Layers Entities Matrix LowerEntities Volumes\n');
 fprintf(fileID,'%d\n',1);
 fprintf(fileID,'escape\n');
 fprintf(fileID,'Mescape\n'); 

 fprintf(fileID,'''Layers Delete Layer0 escape\n'); 

 %MATERIALS DEFINITION

 fprintf(fileID,'Mescape Data Defaults ProblemType yes Examples/problem_type_solid1 escape\n');

 fprintf(fileID,'Mescape Data Materials AssignMaterial Steel Volumes\n');
 fprintf(fileID,'%d\n',2);
 fprintf(fileID,'Mescape\n');
 fprintf(fileID,'Mescape\n');

 fprintf(fileID,'Mescape Data Materials AssignMaterial Concrete Volumes\n');
 fprintf(fileID,'%d\n',1);
 fprintf(fileID,'Mescape\n');
 fprintf(fileID,'Mescape\n');
 
 % ROTATION

 fprintf(fileID,'Mescape Utilities Move Volumes MaintainLayers Rotation FNoJoin 0.0,1.0,0.0 FNoJoin 0.0,0.0,0.0 90\n');
 fprintf(fileID,'1 2\n');
 fprintf(fileID,'Mescape\n');

 %GEOMETRY MESH

 fprintf(fileID,'Mescape Meshing SemiStructured Volumes\n');
 fprintf(fileID,'%d\n',prisvol);
 fprintf(fileID,'%d\n',1);
 fprintf(fileID,'%d\n',2);
 fprintf(fileID,'escape\n');
 fprintf(fileID,'escape\n');

 fprintf(fileID,'Mescape Meshing ElemType Hexahedra\n');
 fprintf(fileID,'%d\n',1);
 fprintf(fileID,'%d\n',2);
 fprintf(fileID,'escape\n');

 fprintf(fileID,'Mescape Meshing Generate\n');
 fprintf(fileID,'%4.3f MeshingParametersFrom=Model\n',mden);
 fprintf(fileID,'Mescape Meshing MeshView\n'); 
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function PRINT_FIGURE2(Ax,d,mden,prisvol,fileID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AF = 0.25*pi*d^2;% [m^2]
%AC = 2*L*H/F;    % [m^2]
AC = Ax*Ax;    % [m^2]
INC_PER = 2*AF/AC;

Ay=Ax;

disp(Ax);
disp(d);
disp(INC_PER);

if( sqrt(2)*d >= Ax)
  display('WARNING inclusion bigger than representative volume dimmensions')  
 return;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CREATE GEOMETRY

fprintf(fileID,'*****FIGURE2_CREATION\n');
fprintf(fileID,'Mescape Geometry Create Object rectangle\n');
fprintf(fileID,'%d %d\n',[0,0]);
fprintf(fileID,'%d %d\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Object circle\n');
fprintf(fileID,'%d %d\n',[0,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'%2.1f %2.1f %2.1f\n',[0.0,0.0,1.0]);
fprintf(fileID,'%d\n',0.5*d);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Object circle\n');
fprintf(fileID,'%d %d\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'%2.1f %2.1f %2.1f\n',[0.0,0.0,1.0]);
fprintf(fileID,'%d\n',0.5*d);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create IntSolid2d substract\n');
fprintf(fileID,'%d\n',1);
fprintf(fileID,'escape\n');
fprintf(fileID,'%d %d\n',[2,3]);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax,0.5*Ay-0.5*d]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax-0.5*d,0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create NurbsSurface\n');
fprintf(fileID,'%d %d %d\n',[20,21,22]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces MaintainLayers Mirror FJoin 4 FJoin 7 TwoDim\n');
fprintf(fileID,'%d %d \n',[5,6]);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces MaintainLayers Mirror FJoin 5 FJoin 2 TwoDim\n');
fprintf(fileID,'%d %d %d %d\n',[5,6,7,8]);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create NurbsSurface\n');
fprintf(fileID,'%d %d %d %d\n',[11,24,31,37]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,%5.4f\n',Ax);
fprintf(fileID,'%d %d %d %d %d %d %d %d %d\n',[5,6,7,8,9,10,11,12,13]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

%LAYERS DEFINITION

fprintf(fileID,'''Layers New Fibre escape\n');
fprintf(fileID,'''Layers Entities Fibre LowerEntities Volumes\n');
fprintf(fileID,'%d\n',4);
fprintf(fileID,'%d\n',9);
fprintf(fileID,'%d\n',8);
fprintf(fileID,'%d\n',2);
fprintf(fileID,'%d\n',6);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'''Layers New Matrix escape\n');
fprintf(fileID,'''Layers Entities Matrix LowerEntities Volumes\n');
fprintf(fileID,'%d\n',1);
fprintf(fileID,'%d\n',3);
fprintf(fileID,'%d\n',5);
fprintf(fileID,'%d\n',7);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'''Layers Delete Layer0 escape\n');

 %MATERIALS DEFINITION

 fprintf(fileID,'Mescape Data Defaults ProblemType yes Examples/problem_type_solid1 escape\n');

 fprintf(fileID,'Mescape Data Materials AssignMaterial Steel Volumes\n');
 fprintf(fileID,'%d\n',4);
 fprintf(fileID,'%d\n',9);
 fprintf(fileID,'%d\n',8);
 fprintf(fileID,'%d\n',2);
 fprintf(fileID,'%d\n',6);
 fprintf(fileID,'Mescape\n');
 fprintf(fileID,'Mescape\n');

 fprintf(fileID,'Mescape Data Materials AssignMaterial Concrete Volumes\n');
 fprintf(fileID,'%d\n',1);
 fprintf(fileID,'%d\n',3);
 fprintf(fileID,'%d\n',5);
 fprintf(fileID,'%d\n',7);
 fprintf(fileID,'Mescape\n');
 fprintf(fileID,'Mescape\n');

 
 % ROTATION

 fprintf(fileID,'Mescape Utilities Move Volumes MaintainLayers Rotation FNoJoin 0.0,1.0,0.0 FNoJoin 0.0,0.0,0.0 90\n');
 fprintf(fileID,'1 2 3 4 5 6 7 8 9\n');
 fprintf(fileID,'Mescape\n');
 
 
%GEOMETRY MESH

fprintf(fileID,'Mescape Meshing SemiStructured Volumes\n');
fprintf(fileID,'%d\n',prisvol);
for i =1:9
fprintf(fileID,'%d\n',i);
end
fprintf(fileID,'escape\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Meshing ElemType Hexahedra\n');
for i =1:9
fprintf(fileID,'%d\n',i);
end
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Meshing Generate\n');
fprintf(fileID,'%4.3f MeshingParametersFrom=Model\n',mden);
fprintf(fileID,'Mescape Meshing MeshView\n');

end

function PRINT_FIGURE3(Ax,mden,prisvol,fileID)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
d=Ax;
Ay=Ax;
INC_PER = 0.25*pi;
%AC = d*d;
%AF = 0.25*pi*Ax^2;% [m^2]

disp(Ax);
disp(d);
disp(INC_PER);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE GEOMETRY

fprintf(fileID,'*****FIGURE2_CREATION\n');
fprintf(fileID,'Mescape Geometry Create Object rectangle\n');
fprintf(fileID,'%d %d\n',[0,0]);
fprintf(fileID,'%d %d\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Object circle\n');
fprintf(fileID,'%d %d\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'%2.1f %2.1f %2.1f\n',[0.0,0.0,1.0]);
fprintf(fileID,'%d\n',0.5*Ay);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create IntSolid2d substract\n');
fprintf(fileID,'%d\n',1);
fprintf(fileID,'escape\n');
fprintf(fileID,'%d %d\n',[2,3]);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d,%d,0to0.075\n',[0,0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create NurbsSurface\n');
fprintf(fileID,'%d %d %d\n',[9,10,11]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces MaintainLayers Mirror FJoin 4 FJoin 1 TwoDim\n');
fprintf(fileID,'%d %d \n',[3,4]);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces MaintainLayers Mirror FJoin 1 FJoin 2 TwoDim\n');
fprintf(fileID,'%d %d %d %d\n',[5,6,3,4]);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,%5.4f\n',Ax);
fprintf(fileID,'%d %d %d %d %d %d %d %d\n',[5,6,7,8,9,10,3,4]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

%LAYERS DEFINITION

fprintf(fileID,'''Layers New Fibre escape\n');
fprintf(fileID,'''Layers Entities Fibre LowerEntities Volumes\n');
fprintf(fileID,'%d\n',4);
fprintf(fileID,'%d\n',8);
fprintf(fileID,'%d\n',2);
fprintf(fileID,'%d\n',6);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'''Layers New Matrix escape\n');
fprintf(fileID,'''Layers Entities Matrix LowerEntities Volumes\n');
fprintf(fileID,'%d\n',1);
fprintf(fileID,'%d\n',3);
fprintf(fileID,'%d\n',5);
fprintf(fileID,'%d\n',7);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'''Layers Delete Layer0 escape\n');

%MATERIALS DEFINITION

fprintf(fileID,'Mescape Data Defaults ProblemType yes Examples/problem_type_solid1 escape\n');

fprintf(fileID,'Mescape Data Materials AssignMaterial Steel Volumes\n');
fprintf(fileID,'%d\n',4);
fprintf(fileID,'%d\n',8);
fprintf(fileID,'%d\n',2);
fprintf(fileID,'%d\n',6);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Data Materials AssignMaterial Concrete Volumes\n');
fprintf(fileID,'%d\n',3);
fprintf(fileID,'%d\n',7);
fprintf(fileID,'%d\n',5);
fprintf(fileID,'%d\n',1);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');


 % ROTATION

 fprintf(fileID,'Mescape Utilities Move Volumes MaintainLayers Rotation FNoJoin 0.0,1.0,0.0 FNoJoin 0.0,0.0,0.0 90\n');
 fprintf(fileID,'1 2 3 4 5 6 7 8\n');
 fprintf(fileID,'Mescape\n');


%GEOMETRY MESH

fprintf(fileID,'Mescape Meshing SemiStructured Volumes\n');
fprintf(fileID,'%d\n',prisvol);
for i =1:8
fprintf(fileID,'%d\n',i);
end
fprintf(fileID,'escape\n');
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape Meshing ElemType Hexahedra\n');
for i =1:8
fprintf(fileID,'%d\n',i);
end
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Meshing Generate\n');
fprintf(fileID,'%4.3f MeshingParametersFrom=Model\n',mden);
fprintf(fileID,'Mescape Meshing MeshView\n');


end

function PRINT_FIGURE4(Ax,mden,prisvol,fileID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 d=Ax;
 Ay=Ax;
 INC_PER = 0.25*pi;
 %AC = d*d;
 %AF = 0.25*pi*Ax^2;% [m^2]


disp(Ax);
disp(d);
disp(INC_PER);


%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CREATE GEOMETRY

fprintf(fileID,'*****FIGURE2_CREATION\n');
fprintf(fileID,'Mescape Geometry Create Object rectangle\n');
fprintf(fileID,'%d %d\n',[-0.5*Ax,0]);
fprintf(fileID,'%d %d\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Object circle\n');
fprintf(fileID,'%d %d\n',[0.5*Ax,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'%2.1f %2.1f %2.1f\n',[0.0,0.0,1.0]);
fprintf(fileID,'%d\n',0.5*Ay);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Object circle\n');
fprintf(fileID,'%d %d\n',[-0.5*Ax,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'%2.1f %2.1f %2.1f\n',[0.0,0.0,1.0]);
fprintf(fileID,'%d\n',0.5*Ay);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create IntSolid2d substract\n');
fprintf(fileID,'%d\n',1);
fprintf(fileID,'escape\n');
fprintf(fileID,'%d %d\n',[2,3]);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Point\n');
fprintf(fileID,'%d %d\n',[1.0,0.0]);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces MaintainLayers Mirror FJoin 8 FJoin 7 TwoDim\n');
fprintf(fileID,'%d\n',5);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax,0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d,%d,0to0.075\n',[0.5*Ax,-0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d,%d,0to0.075\n',[-0.5*Ax,0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d,%d,0to0.075\n',[-0.5*Ax,-0.5*Ay]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create NurbsSurface\n');
fprintf(fileID,'%d %d %d\n',[12,15,17]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Geometry Create NurbsSurface\n');
fprintf(fileID,'%d %d %d\n',[13,16,18]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,%5.4f\n',Ax);
fprintf(fileID,'%d %d %d %d\n',[5,6,7,8]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

%LAYERS DEFINITION

fprintf(fileID,'''Layers New Fibre escape\n');
fprintf(fileID,'''Layers Entities Fibre LowerEntities Volumes\n');
fprintf(fileID,'%d\n',4);
fprintf(fileID,'%d\n',3);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'''Layers New Matrix escape\n');
fprintf(fileID,'''Layers Entities Matrix LowerEntities Volumes\n');
fprintf(fileID,'%d\n',1);
fprintf(fileID,'%d\n',2);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');


%MATERIALS DEFINITION

fprintf(fileID,'Mescape Data Defaults ProblemType yes Examples/problem_type_solid1 escape\n');

fprintf(fileID,'Mescape Data Materials AssignMaterial Steel Volumes\n');
fprintf(fileID,'%d\n',4);
fprintf(fileID,'%d\n',3);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Data Materials AssignMaterial Concrete Volumes\n');
fprintf(fileID,'%d\n',2);
fprintf(fileID,'%d\n',1);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');


 % ROTATION

 fprintf(fileID,'Mescape Utilities Move Volumes MaintainLayers Rotation FNoJoin 0.0,1.0,0.0 FNoJoin 0.0,0.0,0.0 90\n');
 fprintf(fileID,'1 2 3 4\n');
 fprintf(fileID,'Mescape\n');


%GEOMETRY MESH

fprintf(fileID,'Mescape Meshing SemiStructured Volumes\n');
fprintf(fileID,'%d\n',prisvol);
for i =1:4
fprintf(fileID,'%d\n',i);
end
fprintf(fileID,'escape\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Meshing ElemType Hexahedra\n');
for i =1:4
fprintf(fileID,'%d\n',i);
end
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Meshing Generate\n');
fprintf(fileID,'%4.3f MeshingParametersFrom=Model\n',mden);
fprintf(fileID,'Mescape Meshing MeshView\n');

end

function PRINT_FIGURE5(Ad,mden,prisvol,fileID)

 
 INC_PER = 0.5*pi/sqrt(3);
 %AC = 2*Ad^2*sqrt(3);
 %AF = 0.25*pi*Ad^2;% [m^2]
 
disp(Ad);
disp(INC_PER);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CREATE GEOMETRY

fprintf(fileID,'*****FIGURE2_CREATION\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d %d 0to0.075\n',[0,0]);
fprintf(fileID,'%d %d 0to0.075\n',[0.5*Ad,0]);
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d %d 0to0.075\n',[0.5*Ad,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d %d 0to0.075\n',[0.25*Ad,0.5*Ad*0.5*1.73205080757]);
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d %d 0to0.075\n',[0.25*Ad,0.5*Ad*0.5*1.73205080757]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d %d 0to0.075\n',[0,Ad*0.5*1.73205080757]); %3.4641
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d %d 0to0.075\n',[0,Ad*0.5*1.73205080757]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d %d 0to0.075\n',[-0.25*Ad,0.5*Ad*0.5*1.73205080757]);
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d %d 0to0.075\n',[-0.25*Ad,0.5*Ad*0.5*1.73205080757]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d %d 0to0.075\n',[-0.5*Ad,0]);
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Line\n');
fprintf(fileID,'%d %d 0to0.075\n',[-0.5*Ad,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d %d 0to0.075\n',[0,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Arc\n');
fprintf(fileID,'%d %d 0to0.075\n',[0,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d %d 0to0.075\n',[0.5*Ad*0.5*1.4142-0.5*Ad,0.5*Ad*0.5*1.4142]);  %-0.5858
fprintf(fileID,'%d %d 0to0.075\n',[-0.25*Ad,0.5*Ad*0.5*1.73205080757]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Arc\n');
fprintf(fileID,'%d %d 0to0.075\n',[0,0]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d %d 0to0.075\n',[-0.5*Ad*0.5*1.4142+0.5*Ad,0.5*Ad*0.5*1.4142]);
fprintf(fileID,'%d %d 0to0.075\n',[0.25*Ad,0.5*Ad*0.5*1.73205080757]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create Arc\n');
fprintf(fileID,'%d %d 0to0.075\n',[-0.25*Ad,0.5*Ad*0.5*1.73205080757]);
fprintf(fileID,'old\n');
fprintf(fileID,'%d %d 0to0.075\n',[0,Ad*0.5*1.73205080757-0.5*Ad]);   
fprintf(fileID,'%d %d 0to0.075\n',[0.25*Ad,0.5*Ad*0.5*1.73205080757]);
fprintf(fileID,'old\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Geometry Create NurbsSurface\n');
fprintf(fileID,'%d %d %d\n',[5,6,7]);
fprintf(fileID,'escape\n');
fprintf(fileID,'%d %d %d\n',[1,2,8]);
fprintf(fileID,'escape\n');
fprintf(fileID,'%d %d %d\n',[3,4,9]);
fprintf(fileID,'escape\n');
fprintf(fileID,'%d %d %d\n',[7,8,9]);
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces MaintainLayers Mirror FJoin 1 FJoin 2 TwoDim\n');
fprintf(fileID,'%d %d %d %d\n',[1,2,3,4]);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Utilities Copy Surfaces DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,%5.4f\n',Ad);
fprintf(fileID,'%d %d %d %d %d %d %d %d\n',[1,2,3,4,5,6,7,8]);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

%LAYERS DEFINITION

fprintf(fileID,'''Layers New Fibre escape\n');
fprintf(fileID,'''Layers Entities Fibre LowerEntities Volumes\n');
fprintf(fileID,'%d\n',1);
fprintf(fileID,'%d\n',2);
fprintf(fileID,'%d\n',3);
fprintf(fileID,'%d\n',5);
fprintf(fileID,'%d\n',6);
fprintf(fileID,'%d\n',7);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'''Layers New Matrix escape\n');
fprintf(fileID,'''Layers Entities Matrix LowerEntities Volumes\n');
fprintf(fileID,'%d\n',4);
fprintf(fileID,'%d\n',8);
fprintf(fileID,'escape\n');
fprintf(fileID,'Mescape\n');

%MATERIALS DEFINITION

fprintf(fileID,'Mescape Data Defaults ProblemType yes Examples/problem_type_solid1 escape\n');

fprintf(fileID,'Mescape Data Materials AssignMaterial Steel Volumes\n');
fprintf(fileID,'%d\n',1);
fprintf(fileID,'%d\n',2);
fprintf(fileID,'%d\n',3);
fprintf(fileID,'%d\n',5);
fprintf(fileID,'%d\n',6);
fprintf(fileID,'%d\n',7);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');

fprintf(fileID,'Mescape Data Materials AssignMaterial Concrete Volumes\n');
fprintf(fileID,'%d\n',4);
fprintf(fileID,'%d\n',8);
fprintf(fileID,'Mescape\n');
fprintf(fileID,'Mescape\n');


 % ROTATION

 fprintf(fileID,'Mescape Utilities Move Volumes MaintainLayers Rotation FNoJoin 0.0,1.0,0.0 FNoJoin 0.0,0.0,0.0 90\n');
 fprintf(fileID,'1 2 3 4 5 6 7 8\n');
 fprintf(fileID,'Mescape\n');

% MESH DEFINITION

fprintf(fileID,'Mescape Meshing SemiStructured Volumes\n');
fprintf(fileID,'%d\n',prisvol);
for i =1:8
fprintf(fileID,'%d\n',i);
end
fprintf(fileID,'escape\n');
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Meshing ElemType Hexahedra\n');
for i =1:8
fprintf(fileID,'%d\n',i);
end
fprintf(fileID,'escape\n');

fprintf(fileID,'Mescape Meshing Generate\n');
fprintf(fileID,'%4.3f MeshingParametersFrom=Model\n',mden);
fprintf(fileID,'Mescape Meshing MeshView\n');

end
%{
\end{lstlisting}
%}