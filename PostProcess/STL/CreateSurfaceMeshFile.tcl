set input_post_res "/home/ton/test_micro23.flavia.res"
set output_gid_project_name "/home/ton/test_micro_project.gid"
set mesh_element_size "0.0707107"
set mesh_name "hmmmm"

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
GiD_Process Mescape Files WriteMesh /home/ton/GiDx64/gid-16.1.2d/$mesh_name.msh
GiD_Process Mescape Files Save $output_gid_project_name escape