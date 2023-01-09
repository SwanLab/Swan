proc CreateSurfaceSTL {input_post_res output_gid_project_name} {
    set input_post_res "/home/ton/test_micro23.flavia.res"
    set output_gid_project_name "/home/ton/test_micro_project.gid"
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
     GiD_Process Mescape Geometry Create NurbsSurface 463 540 541 618 619 escape
     GiD_Process Mescape Files Save $output_gid_project_name escape

}