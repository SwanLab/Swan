GiD_Process Mescape Postprocess
GiD_Process Mescape Files Read $input_post_res escape
GiD_Process Mescape Results Options ExtractBoundaries escape
GiD_Process Mescape Results IsoSurfaces Exact DesignVar1 1 0.0 escape
GiD_Process Mescape Results IsoSurfaces TurnIntoCuts escape
GiD_Process Mescape DoCut ConvertToSets ALL_CUTSETS escape
foreach set_name {"Set0 IsoSurf DesignVar1 = 0 WORKPIECE" "WORKPIECE boundary"} {
GidUtils::CreateMeshFromSet $set_name Layer0
}
GiD_Process Mescape Preprocess Yes
GiD_Process Mescape Geometry Edit ReConstruction OneLineForEachElement 1:end escape
GiD_Process Mescape Geometry Edit JoinLines 1:end escape escape
GiD_Process Mescape Geometry Create NurbsSurface 1:end escape
GiD_Process Mescape Meshing Generate Yes $mesh_element_size MeshingParametersFrom=Preferences
GiD_Process Mescape Files WriteMesh $gidpath$mesh_name.msh