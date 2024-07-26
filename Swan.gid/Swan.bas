%==================================================================
%                        General Data File
% Title: *GenData(1)
% Units: *GenData(2)
% Dimensions: *GenData(3)
% Type of problem: *GenData(4)
% Type of Phisics: *GenData(5)
% Micro/Macro: *GenData(6)
%
%==================================================================

%% Data

Data_prb = {
'*GenData(1)';
'*GenData(2)';
'*GenData(3)';
'*GenData(4)';
'*GenData(5)';
'*GenData(6)';
};

%% Coordinates
% Node		X		Y		Z

gidcoord = [
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
];

%% Conectivities
*set elems(triangle)
*if(nelem(triangle)>0)
% Element	Node(1)		Node(2)		Node(3)		Material
*endif
*set elems(quadrilateral)
*if(nelem(quadrilateral)>0)
% Element	Node(1)		Node(2)		Node(3)		Node(4)		Material
*endif
*set elems(tetrahedra)
*if(nelem(tetrahedra)>0)
% Element	Node(1)		Node(2)		Node(3)		Node(4)		Material
*endif
*set elems(hexahedra)
*if(nelem(hexahedra)>0)
% Element	Node(1)		Node(2)		Node(3)		Node(4)		Node(5)		Node(6)		Node(7)		Node(8)		Material
*endif

gidlnods = [
*set elems(all)
*loop elems
*ElemsNum *ElemsConec *ElemsMat
*end elems
];

%% Variable Prescribed
% Node	    Dimension		Value

lnodes = [
*if(strcmp(GenData(3),"2D")==0)
*Set Cond Point_Constraint *nodes
*Add Cond Line_Constraint *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Ux) 
*NodesNum 2 *Cond(Uy) 
*end nodes
*Set Cond x_direction_point_constraint *nodes
*Add Cond x_direction_Line_Constraint *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Ux) 
*end nodes
*Set Cond y_direction_point_constraint *nodes
*Add Cond y_direction_Line_Constraint *nodes
*loop nodes *OnlyInCond
*NodesNum 2 *Cond(Uy) 
*end nodes
*tcl(proc evaluateConditionFunction { fnc x y z} { expr $fnc})*\
*tcl(proc evaluateElementConditionFunction { fnc i n} { expr $fnc})*\
*tcl(proc quitaDollar { fnc} { regsub -all {\$} $fnc {}})*\
*Set Cond Line_Constraint_function *nodes
*Add Cond Surface_Constraint_function *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *tcl(evaluateConditionFunction *cond(Ux) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*NodesNum 2 *tcl(evaluateConditionFunction *cond(Uy) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*end nodes
*elseif(strcmp(GenData(3),"3D")==0)
*Set Cond Point_Constraint *nodes
*Add Cond Line_Constraint *nodes
*Add Cond Surface_Constraint *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Ux) 
*NodesNum 2 *Cond(Uy)
*NodesNum 3 *Cond(Uz)
*end nodes
*Set Cond x_direction_point_constraint *nodes
*Add Cond x_direction_Line_Constraint *nodes
*Add Cond x_direction_Surface_Constraint *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Ux) 
*end nodes
*Set Cond y_direction_point_constraint *nodes
*Add Cond y_direction_Line_Constraint *nodes
*Add Cond y_direction_Surface_Constraint *nodes
*loop nodes *OnlyInCond
*NodesNum 2 *Cond(Uy) 
*end nodes
*Set Cond z_direction_point_constraint *nodes
*Add Cond z_direction_Line_Constraint *nodes
*Add Cond z_direction_Surface_Constraint *nodes
*loop nodes *OnlyInCond
*NodesNum 3 *Cond(Uz) 
*end nodes
*tcl(proc evaluateConditionFunction { fnc x y z} { expr $fnc})*\
*tcl(proc evaluateElementConditionFunction { fnc i n} { expr $fnc})*\
*tcl(proc quitaDollar { fnc} { regsub -all {\$} $fnc {}})*\
*Set Cond Line_Constraint_function *nodes
*Add Cond Surface_Constraint_function *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *tcl(evaluateConditionFunction *cond(Ux) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*NodesNum 2 *tcl(evaluateConditionFunction *cond(Uy) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*NodesNum 3 *tcl(evaluateConditionFunction *cond(Uz) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*end nodes
*endif
];

%% Force Prescribed
% Node		Dimension		Value

pointload_complete = [
*if(strcmp(GenData(3),"2D")==0)
*Set Cond Point_Force *nodes
*Add Cond Line_Load *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Fx) 
*NodesNum 2 *Cond(Fy) 
*end nodes
*Set Cond x_direction_point_force *nodes
*Add Cond x_direction_Line_Load *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Fx) 
*end nodes
*Set Cond y_direction_point_force *nodes
*Add Cond y_direction_Line_Load *nodes
*loop nodes *OnlyInCond
*NodesNum 2 *Cond(Fy) 
*end nodes
*tcl(proc evaluateConditionFunction { fnc x y z} { expr $fnc})*\
*tcl(proc evaluateElementConditionFunction { fnc i n} { expr $fnc})*\
*tcl(proc quitaDollar { fnc} { regsub -all {\$} $fnc {}})*\
*Set Cond Line_Load_function *nodes *CanRepeat
*Add Cond Surface_Load_function *nodes *CanRepeat
*loop nodes *OnlyInCond
*NodesNum 1 *tcl(evaluateConditionFunction *cond(Fx) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*NodesNum 2 *tcl(evaluateConditionFunction *cond(Fy) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*end nodes
*Set Elems(triangle)
*Add elems(quadrilateral)
*tcl(proc Computelength { P xL1 yL1 xL2 yL2 } { return [expr {$P*0.5*sqrt(($xL1-$xL2)*($xL1-$xL2)+($yL1-$yL2)*($yL1-$yL2))}] })*\
*Set Cond Line_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var resPx=tcl(Computelength *cond(PressureX,real) *xL1 *yL1 *xL2 *yL2)
*Set var resPy=tcl(Computelength *cond(PressureY,real) *xL1 *yL1 *xL2 *yL2)
*ElemsConec(*LocalNodes(1),int) 1 *resPx
*ElemsConec(*LocalNodes(1),int) 2 *resPy
*ElemsConec(*LocalNodes(2),int) 1 *resPx
*ElemsConec(*LocalNodes(2),int) 2 *resPy
*endif
*end elems
*Set Cond x_direction_Line_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var resPx=tcl(Computelength *cond(PressureX,real) *xL1 *yL1 *xL2 *yL2)
*ElemsConec(*LocalNodes(1),int) 1 *resPx
*ElemsConec(*LocalNodes(2),int) 1 *resPx
*endif
*end elems
*Set Cond y_direction_Line_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var resPy=tcl(Computelength *cond(PressureY,real) *xL1 *yL1 *xL2 *yL2)
*ElemsConec(*LocalNodes(1),int) 2 *resPy
*ElemsConec(*LocalNodes(2),int) 2 *resPy
*endif
*end elems
*elseif(strcmp(GenData(3),"3D")==0)
*Set Cond Point_Force *nodes
*Add Cond Line_Load *nodes
*Add Cond Surface_Load *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Fx) 
*NodesNum 2 *Cond(Fy)
*NodesNum 3 *Cond(Fz)
*end nodes
*Set Cond x_direction_point_force *nodes
*Add Cond x_direction_Line_Load *nodes
*Add Cond x_direction_Surface_Load *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Fx) 
*end nodes
*Set Cond y_direction_point_force *nodes
*Add Cond y_direction_Line_Load *nodes
*Add Cond y_direction_Surface_Load *nodes
*loop nodes *OnlyInCond
*NodesNum 2 *Cond(Fy) 
*end nodes
*Set Cond z_direction_point_force *nodes
*Add Cond z_direction_Line_Load *nodes
*Add Cond z_direction_Surface_Load *nodes
*loop nodes *OnlyInCond
*NodesNum 3 *Cond(Fz) 
*end nodes
*tcl(proc evaluateConditionFunction { fnc x y z} { expr $fnc})*\
*tcl(proc evaluateElementConditionFunction { fnc i n} { expr $fnc})*\
*tcl(proc quitaDollar { fnc} { regsub -all {\$} $fnc {}})*\
*Set Cond Line_Load_function *nodes *CanRepeat
*Add Cond Surface_Load_function *nodes *CanRepeat
*loop nodes *OnlyInCond
*NodesNum 1 *tcl(evaluateConditionFunction *cond(Fx) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*NodesNum 2 *tcl(evaluateConditionFunction *cond(Fy) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*NodesNum 3 *tcl(evaluateConditionFunction *cond(Fz) *NodesCoord(1) *NodesCoord(2) *NodesCoord(3))
*end nodes
*Set Elems(tetrahedra)
*Add elems(hexahedra)
*tcl(proc Computelength { P xL1 yL1 zL1 xL2 yL2 zL2} { return [expr {$P*0.5*sqrt(($xL1-$xL2)*($xL1-$xL2)+($yL1-$yL2)*($yL1-$yL2)+($zL1-$zL2)*($zL1-$zL2))}] })*\
*Set Cond Line_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var resPx=tcl(Computelength *cond(PressureX,real) *xL1 *yL1 *zL1 *xL2 *yL2 *zL2)
*Set var resPy=tcl(Computelength *cond(PressureY,real) *xL1 *yL1 *zL1 *xL2 *yL2 *zL2)
*Set var resPz=tcl(Computelength *cond(PressureZ,real) *xL1 *yL1 *zL1 *xL2 *yL2 *zL2)
*ElemsConec(*LocalNodes(1),int) 1 *resPx
*ElemsConec(*LocalNodes(1),int) 2 *resPy
*ElemsConec(*LocalNodes(1),int) 3 *resPz
*ElemsConec(*LocalNodes(2),int) 1 *resPx
*ElemsConec(*LocalNodes(2),int) 2 *resPy
*ElemsConec(*LocalNodes(2),int) 3 *resPz
*endif
*end elems
*Set Cond x_direction_Line_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var resPx=tcl(Computelength *cond(PressureX,real) *xL1 *yL1 *zL1 *xL2 *yL2 *zL2)
*ElemsConec(*LocalNodes(1),int) 1 *resPx
*ElemsConec(*LocalNodes(2),int) 1 *resPx
*endif
*end elems
*Set Cond y_direction_Line_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var resPy=tcl(Computelength *cond(PressureY,real) *xL1 *yL1 *zL1 *xL2 *yL2 *zL2)
*ElemsConec(*LocalNodes(1),int) 2 *resPy
*ElemsConec(*LocalNodes(2),int) 2 *resPy
*endif
*end elems
*Set Cond z_direction_Line_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var resPz=tcl(Computelength *cond(PressureZ,real) *xL1 *yL1 *zL1 *xL2 *yL2 *zL2)
*ElemsConec(*LocalNodes(1),int) 3 *resPz
*ElemsConec(*LocalNodes(2),int) 3 *resPz
*endif
*end elems
*Set elems(tetrahedra)
*tcl(proc ComputeComp1 { v1x v1y v1z v2x v2y v2z} { return [expr {$v1y*$v2z-$v1z*$v2y}] })*\
*tcl(proc ComputeComp2 { v1x v1y v1z v2x v2y v2z} { return [expr {$v1z*$v2x-$v1x*$v2z}] })*\
*tcl(proc ComputeComp3 { v1x v1y v1z v2x v2y v2z} { return [expr {$v1x*$v2y-$v1y*$v2x}] })*\
*tcl(proc ComputeArea { comp1 comp2 comp3} { return [expr {0.5*sqrt($comp1*$comp1+$comp2*$comp2+$comp3*$comp3)}] })*\
*tcl(proc ComputePressure { P Area} { return [expr {$P*$Area/3}] })*\
*Set Cond Surface_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var xL3=NodesCoord(*LocalNodes(3),1,real)
*Set var yL3=NodesCoord(*LocalNodes(3),2,real)
*Set var zL3=NodesCoord(*LocalNodes(3),3,real)
*Set var v1x=xL2-xL1
*Set var v1y=yL2-yL1
*Set var v1z=zL2-zL1
*Set var v2x=xL3-xL1
*Set var v2y=yL3-yL1
*Set var v2z=zL3-zL1
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var resPx=tcl(ComputePressure *cond(PressureX,real) *Area)
*Set var resPy=tcl(ComputePressure *cond(PressureY,real) *Area)
*Set var resPz=tcl(ComputePressure *cond(PressureZ,real) *Area)
*ElemsConec(*LocalNodes(1),int) 1 *resPx
*ElemsConec(*LocalNodes(1),int) 2 *resPy
*ElemsConec(*LocalNodes(1),int) 3 *resPz
*ElemsConec(*LocalNodes(2),int) 1 *resPx
*ElemsConec(*LocalNodes(2),int) 2 *resPy
*ElemsConec(*LocalNodes(2),int) 3 *resPz
*ElemsConec(*LocalNodes(3),int) 1 *resPx
*ElemsConec(*LocalNodes(3),int) 2 *resPy
*ElemsConec(*LocalNodes(3),int) 3 *resPz
*endif
*end elems
*Set Cond x_direction_Surface_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var xL3=NodesCoord(*LocalNodes(3),1,real)
*Set var yL3=NodesCoord(*LocalNodes(3),2,real)
*Set var zL3=NodesCoord(*LocalNodes(3),3,real)
*Set var v1x=xL2-xL1
*Set var v1y=yL2-yL1
*Set var v1z=zL2-zL1
*Set var v2x=xL3-xL1
*Set var v2y=yL3-yL1
*Set var v2z=zL3-zL1
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var resPx=tcl(ComputePressure *cond(PressureX,real) *Area)
*ElemsConec(*LocalNodes(1),int) 1 *resPx
*ElemsConec(*LocalNodes(2),int) 1 *resPx
*ElemsConec(*LocalNodes(3),int) 1 *resPx
*endif
*end elems
*Set Cond y_direction_Surface_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var xL3=NodesCoord(*LocalNodes(3),1,real)
*Set var yL3=NodesCoord(*LocalNodes(3),2,real)
*Set var zL3=NodesCoord(*LocalNodes(3),3,real)
*Set var v1x=xL2-xL1
*Set var v1y=yL2-yL1
*Set var v1z=zL2-zL1
*Set var v2x=xL3-xL1
*Set var v2y=yL3-yL1
*Set var v2z=zL3-zL1
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var resPy=tcl(ComputePressure *cond(PressureY,real) *Area)
*ElemsConec(*LocalNodes(1),int) 2 *resPy
*ElemsConec(*LocalNodes(2),int) 2 *resPy
*ElemsConec(*LocalNodes(3),int) 2 *resPy
*endif
*end elems
*Set Cond z_direction_Surface_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var xL3=NodesCoord(*LocalNodes(3),1,real)
*Set var yL3=NodesCoord(*LocalNodes(3),2,real)
*Set var zL3=NodesCoord(*LocalNodes(3),3,real)
*Set var v1x=xL2-xL1
*Set var v1y=yL2-yL1
*Set var v1z=zL2-zL1
*Set var v2x=xL3-xL1
*Set var v2y=yL3-yL1
*Set var v2z=zL3-zL1
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var resPz=tcl(ComputePressure *cond(PressureZ,real) *Area)
*ElemsConec(*LocalNodes(1),int) 3 *resPz
*ElemsConec(*LocalNodes(2),int) 3 *resPz
*ElemsConec(*LocalNodes(3),int) 3 *resPz
*endif
*end elems
*Set elems(hexahedra)
*tcl(proc ComputeComp1 { v1x v1y v1z v2x v2y v2z} { return [expr {$v1y*$v2z-$v1z*$v2y}] })*\
*tcl(proc ComputeComp2 { v1x v1y v1z v2x v2y v2z} { return [expr {$v1z*$v2x-$v1x*$v2z}] })*\
*tcl(proc ComputeComp3 { v1x v1y v1z v2x v2y v2z} { return [expr {$v1x*$v2y-$v1y*$v2x}] })*\
*tcl(proc ComputeArea { comp1 comp2 comp3} { return [expr {0.5*sqrt($comp1*$comp1+$comp2*$comp2+$comp3*$comp3)}] })*\
*tcl(proc ComputePressure { P Area1 Area2} { return [expr {$P*($Area1+$Area2)/4}] })*\
*Set Cond Surface_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var xL3=NodesCoord(*LocalNodes(3),1,real)
*Set var yL3=NodesCoord(*LocalNodes(3),2,real)
*Set var zL3=NodesCoord(*LocalNodes(3),3,real)
*Set var xL4=NodesCoord(*LocalNodes(4),1,real)
*Set var yL4=NodesCoord(*LocalNodes(4),2,real)
*Set var zL4=NodesCoord(*LocalNodes(4),3,real)
*Set var v1x=xL2-xL1
*Set var v1y=yL2-yL1
*Set var v1z=zL2-zL1
*Set var v2x=xL3-xL1
*Set var v2y=yL3-yL1
*Set var v2z=zL3-zL1
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area1=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var v1x=xL4-xL3
*Set var v1y=yL4-yL3
*Set var v1z=zL4-zL3
*Set var v2x=xL2-xL3
*Set var v2y=yL2-yL3
*Set var v2z=zL2-zL3
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area2=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var resPx=tcl(ComputePressure *cond(PressureX,real) *Area1 *Area2)
*Set var resPy=tcl(ComputePressure *cond(PressureY,real) *Area1 *Area2)
*Set var resPz=tcl(ComputePressure *cond(PressureZ,real) *Area1 *Area2)
*ElemsConec(*LocalNodes(1),int) 1 *resPx
*ElemsConec(*LocalNodes(1),int) 2 *resPy
*ElemsConec(*LocalNodes(1),int) 3 *resPz
*ElemsConec(*LocalNodes(2),int) 1 *resPx
*ElemsConec(*LocalNodes(2),int) 2 *resPy
*ElemsConec(*LocalNodes(2),int) 3 *resPz
*ElemsConec(*LocalNodes(3),int) 1 *resPx
*ElemsConec(*LocalNodes(3),int) 2 *resPy
*ElemsConec(*LocalNodes(3),int) 3 *resPz
*ElemsConec(*LocalNodes(4),int) 1 *resPx
*ElemsConec(*LocalNodes(4),int) 2 *resPy
*ElemsConec(*LocalNodes(4),int) 3 *resPz
*endif
*end elems
*Set Cond x_direction_Surface_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var xL3=NodesCoord(*LocalNodes(3),1,real)
*Set var yL3=NodesCoord(*LocalNodes(3),2,real)
*Set var zL3=NodesCoord(*LocalNodes(3),3,real)
*Set var xL4=NodesCoord(*LocalNodes(4),1,real)
*Set var yL4=NodesCoord(*LocalNodes(4),2,real)
*Set var zL4=NodesCoord(*LocalNodes(4),3,real)
*Set var v1x=xL2-xL1
*Set var v1y=yL2-yL1
*Set var v1z=zL2-zL1
*Set var v2x=xL3-xL1
*Set var v2y=yL3-yL1
*Set var v2z=zL3-zL1
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area1=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var v1x=xL4-xL3
*Set var v1y=yL4-yL3
*Set var v1z=zL4-zL3
*Set var v2x=xL2-xL3
*Set var v2y=yL2-yL3
*Set var v2z=zL2-zL3
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area2=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var resPx=tcl(ComputePressure *cond(PressureX,real) *Area1 *Area2)
*ElemsConec(*LocalNodes(1),int) 1 *resPx
*ElemsConec(*LocalNodes(2),int) 1 *resPx
*ElemsConec(*LocalNodes(3),int) 1 *resPx
*ElemsConec(*LocalNodes(4),int) 1 *resPx
*endif
*end elems
*Set Cond y_direction_Surface_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var xL3=NodesCoord(*LocalNodes(3),1,real)
*Set var yL3=NodesCoord(*LocalNodes(3),2,real)
*Set var zL3=NodesCoord(*LocalNodes(3),3,real)
*Set var xL4=NodesCoord(*LocalNodes(4),1,real)
*Set var yL4=NodesCoord(*LocalNodes(4),2,real)
*Set var zL4=NodesCoord(*LocalNodes(4),3,real)
*Set var v1x=xL2-xL1
*Set var v1y=yL2-yL1
*Set var v1z=zL2-zL1
*Set var v2x=xL3-xL1
*Set var v2y=yL3-yL1
*Set var v2z=zL3-zL1
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area1=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var v1x=xL4-xL3
*Set var v1y=yL4-yL3
*Set var v1z=zL4-zL3
*Set var v2x=xL2-xL3
*Set var v2y=yL2-yL3
*Set var v2z=zL2-zL3
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area2=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var resPy=tcl(ComputePressure *cond(PressureY,real) *Area1 *Area2)
*ElemsConec(*LocalNodes(1),int) 2 *resPy
*ElemsConec(*LocalNodes(2),int) 2 *resPy
*ElemsConec(*LocalNodes(3),int) 2 *resPy
*ElemsConec(*LocalNodes(4),int) 2 *resPy
*endif
*end elems
*Set Cond z_direction_Surface_Force *elems *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*Set var xL1=NodesCoord(*LocalNodes(1),1,real)
*Set var yL1=NodesCoord(*LocalNodes(1),2,real)
*Set var zL1=NodesCoord(*LocalNodes(1),3,real)
*Set var xL2=NodesCoord(*LocalNodes(2),1,real)
*Set var yL2=NodesCoord(*LocalNodes(2),2,real)
*Set var zL2=NodesCoord(*LocalNodes(2),3,real)
*Set var xL3=NodesCoord(*LocalNodes(3),1,real)
*Set var yL3=NodesCoord(*LocalNodes(3),2,real)
*Set var zL3=NodesCoord(*LocalNodes(3),3,real)
*Set var xL4=NodesCoord(*LocalNodes(4),1,real)
*Set var yL4=NodesCoord(*LocalNodes(4),2,real)
*Set var zL4=NodesCoord(*LocalNodes(4),3,real)
*Set var v1x=xL2-xL1
*Set var v1y=yL2-yL1
*Set var v1z=zL2-zL1
*Set var v2x=xL3-xL1
*Set var v2y=yL3-yL1
*Set var v2z=zL3-zL1
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area1=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var v1x=xL4-xL3
*Set var v1y=yL4-yL3
*Set var v1z=zL4-zL3
*Set var v2x=xL2-xL3
*Set var v2y=yL2-yL3
*Set var v2z=zL2-zL3
*Set var comp1=tcl(ComputeComp1 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp2=tcl(ComputeComp2 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var comp3=tcl(ComputeComp3 *v1x *v1y *v1z *v2x *v2y *v2z)
*Set var Area2=tcl(ComputeArea *comp1 *comp2 *comp3)
*Set var resPz=tcl(ComputePressure *cond(PressureZ,real) *Area1 *Area2)
*ElemsConec(*LocalNodes(1),int) 3 *resPz
*ElemsConec(*LocalNodes(2),int) 3 *resPz
*ElemsConec(*LocalNodes(3),int) 3 *resPz
*ElemsConec(*LocalNodes(4),int) 3 *resPz
*endif
*end elems
*endif
];

pointload_adjoint = [
*if(strcmp(GenData(3),"2D")==0)
*Set Cond Point_AdjointLoad *nodes
*Add Cond Line_AdjointLoad *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Fx) 
*NodesNum 2 *Cond(Fy) 
*end nodes
*elseif(strcmp(GenData(3),"3D")==0)
*Set Cond Point_AdjointLoad *nodes
*Add Cond Line_AdjointLoad *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *Cond(Fx) 
*NodesNum 2 *Cond(Fy)
*NodesNum 3 *Cond(Fz)
*end nodes
*endif
];

%% Volumetric Force
% Element	Dim		Force_Dim

Vol_force = [
*if(strcmp(GenData(3),"2D")==0)
*Set Elems(all)
*Set Cond Surface_Volumetric_force *elements
*loop elems *OnlyInCond
*ElemsNum 1 *Cond(Vol_forceX)
*ElemsNum 2 *Cond(Vol_forceY)
*end elems
*Set Cond Surface_x_direction_Volumetric_force *elements
*loop elems *OnlyInCond
*ElemsNum 1 *Cond(Vol_forceX)
*end elems
*Set Cond Surface_y_direction_Volumetric_force *elements
*loop elems *OnlyInCond
*ElemsNum 2 *Cond(Vol_forceY)
*end elems
*elseif(strcmp(GenData(3),"3D")==0)
*Set elems(all)
*Set Cond Volume_Volumetric_force *elements
*loop elems *OnlyInCond
*ElemsNum 1 *Cond(Vol_forceX)
*ElemsNum 2 *Cond(Vol_forceY)
*ElemsNum 3 *Cond(Vol_forceZ)
*end elems
*Set Cond Volume_x_direction_Volumetric_force *elements
*loop elems *OnlyInCond
*ElemsNum 1 *Cond(Vol_forceX)
*end elems
*Set Cond Volume_y_direction_Volumetric_force *elements
*loop elems *OnlyInCond
*ElemsNum 2 *Cond(Vol_forceY)
*end elems
*Set Cond Volume_z_direction_Volumetric_force *elements
*loop elems *OnlyInCond
*ElemsNum 3 *Cond(Vol_forceZ)
*end elems
*endif
];

%% Group Elements
% Element	Group_num

Group = [
*Set Cond Group_value *elements
*loop elems *OnlyInCond
*ElemsNum *Cond(Group)
*end elems
];

%% Initial Holes
% Elements that are considered holes initially
% Element

Initial_holes = [
*if(strcmp(GenData(3),"2D")==0)
*Set elems(all)
*Set Cond Initial_holes *elements
*loop elems *OnlyInCond
*ElemsNum
*end elems
*elseif(strcmp(GenData(3),"3D")==0)
*Set elems(all)
*Set Cond Volume_Initial_holes *elements
*loop elems *OnlyInCond
*ElemsNum
*end elems
*endif
];

%% Boundary Elements
% Elements that can not be removed
% Element

Boundary_elements = [
*if(strcmp(GenData(3),"2D")==0)
*Set elems(all)
*Set Cond Boundary_elements *elems
*loop elems *OnlyInCond
*ElemsNum
*end elems
*elseif(strcmp(GenData(3),"3D")==0)
*Set elems(all)
*Set Cond Volume_Boundary_elements *elems
*loop elems *OnlyInCond
*ElemsNum
*end elems
*endif
];

%% Micro gauss post
%
% Element

Micro_gauss_post = [
*Set Cond Micro_gauss *elements
*loop elems *OnlyInCond
*ElemsNum
*end elems
];


%% Micro Slave-Master
% Nodes that are Slaves
% Nodes	     Value (1-Slave,0-Master)

Micro_slave = [
*if(strcmp(Gendata(6),"MICRO")==0)
*if(strcmp(Gendata(3),"2D")==0)
*Set Cond Line_Micro_Slave *nodes
*loop nodes *OnlyInCond
*NodesNum *Cond(Slave)
*end nodes
*elseif(strcmp(Gendata(3),"3D")==0)
*Set Cond Surface_Micro_Slave *nodes
*loop nodes *OnlyInCond
*NodesNum *Cond(Slave)
*end nodes
*endif
*endif
];

%% Nodes solid
% Nodes that must remain 
% Nodes

if ~isempty(pointload_complete)
    nodesolid = unique(pointload_complete(:,1));
end

%% External border Elements
% Detect the elements that define the edge of the domain
% Element    	   Node(1) 	  Node(2)

External_border_elements = [
*if(strcmp(GenData(3),"2D")==0)
*Set elems(all)
*Set Cond Line_Border_elements *elems *CanRepeat
*loop elems *OnlyInCond
*ElemsNum *GlobalNodes
*end elems
*elseif(strcmp(GenData(3),"3D")==0)
*Set Elems(all)
*Set Cond Surface_Border_elements *elems *CanRepeat
*loop elems *OnlyInCond
*ElemsNum *GlobalNodes
*end elems
*endif
];

%% External border Nodes
% Detect the nodes that define the edge of the domain
% Node

External_border_nodes = [
*if(strcmp(GenData(3),"2D")==0)
*Set elems(all)
*Set Cond Line_Border_nodes *nodes
*loop nodes *OnlyInCond
*NodesNum
*end nodes
*elseif(strcmp(GenData(3),"3D")==0)
*Set elems(all)
*Set Cond Surface_Border_nodes *nodes
*loop nodes *OnlyInCond
*NodesNum
*end nodes
*endif
];

%% Materials
% Materials that have been used
% Material_Num	      Mat_density	Young_Modulus	Poisson

Materials = [
*loop materials
*format "%4i%13.5e%13.5e%13.5e%13.5e"
*MatNum *MatProp(Density,real) *MatProp(Young,real) *MatProp(Poisson,real)
*end materials
];
