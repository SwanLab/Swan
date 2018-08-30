# Macros file for GiD v1.0

# Created by GiD version 13.0.1

#[_ ""]
set macrosdata(GeomSTL_PG,Active) 1
set macrosdata(GeomSTL_PG,ModificationDate) {2017-04-08 20:25:41}
set macrosdata(GeomSTL_PG,Accelerators) {}
set macrosdata(GeomSTL_PG,Group) {}
set macrosdata(GeomSTL_PG,Number) 1
set macrosdata(GeomSTL_PG,IsStandard) 0
set macrosdata(GeomSTL_PG,InToolbar) 0
set macrosdata(GeomSTL_PG,Icon) {}
set macrosdata(GeomSTL_PG,PrePost) prepost
set macrosdata(GeomSTL_PG,CreationDate) {2017-04-08 20:15:17}
set macrosdata(GeomSTL_PG,Description) {}

proc GeomSTL_PG {res_file ascii_file geom_file} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    GiD_Process Mescape Results IsoSurfaces DisplayStyle ContourFillColor 
    GiD_Process Mescape Results IsoSurfaces Exact {Gamma FUNCTION NODAL } 1 0.0
    GiD_Process Mescape 'Redraw
    GiD_Process Mescape Results isosurfaces TurnIntoCuts Mescape DoCut TurnIntoMeshes
    GiD_Process Mescape Results Options ExtractBoundaries
    GiD_Process Mescape Utilities Delete SurfaceSets {WORKPIECE 1} Yes
    GiD_Process Mescape Utilities Delete CutSets {IsoSurf {Gamma FUNCTION NODAL } = 0.5 WORKPIECE 1} Yes
    GiD_Process Mescape files saveall allmeshessets ${ascii_file}
    GiD_Process Mescape Preprocess Yes
    GiD_Process Mescape Files MeshRead "${ascii_file}.msh"
    GiD_Process Mescape Files WriteForBAS {C:/Program Files/GiD/GiD 13.0.1/templates/DXF.bas} ${geom_file}
    GiD_Process Mescape Files DxfRead -collapse:1 -ignorelayer:1 -tolerance:0.0001 -autotol:1 -- ${geom_file}
    GiD_Process Mescape Geometry Create IntMultLines InvertSelection escape Yes
    GiD_Process Mescape Geometry Delete AllTypes
}

#[_ ""]
set macrosdata(PNGNoSym_PG,Description) {}
set macrosdata(PNGNoSym_PG,PrePost) prepost
set macrosdata(PNGNoSym_PG,Group) {}
set macrosdata(PNGNoSym_PG,InToolbar) 0
set macrosdata(PNGNoSym_PG,CreationDate) {2017-04-08 16:15:09}
set macrosdata(PNGNoSym_PG,ModificationDate) {2017-04-08 16:21:16}
set macrosdata(PNGNoSym_PG,Active) 1
set macrosdata(PNGNoSym_PG,IsStandard) 0
set macrosdata(PNGNoSym_PG,Icon) {}
set macrosdata(PNGNoSym_PG,Accelerators) {}
set macrosdata(PNGNoSym_PG,Number) 2

proc PNGNoSym_PG {res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    GiD_Process Mescape Results ContOptions NumberOfColor 50
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Standard 
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor Standard 
    GiD_Process Mescape Results ContOptions ColorRamp Tangent 
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor White 
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor Black 
    GiD_Process Mescape Results ContOptions ColorRamp Tangent 
    GiD_Process Mescape View ReadView "${view_file}"
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape Results ContourFill {Gamma FUNCTION NODAL }
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
    GiD_Process Mescape 'Hardcopy PNG "${png_file}"
    GiD_Process Mescape quit No escape
}

#[_ ""]
set macrosdata(PNGNoSym,IsStandard) 0
set macrosdata(PNGNoSym,CreationDate) {2017-04-08 15:17:47}
set macrosdata(PNGNoSym,Active) 1
set macrosdata(PNGNoSym,PrePost) prepost
set macrosdata(PNGNoSym,Accelerators) {}
set macrosdata(PNGNoSym,Number) 3
set macrosdata(PNGNoSym,Group) {}
set macrosdata(PNGNoSym,Description) {}
set macrosdata(PNGNoSym,Icon) {}
set macrosdata(PNGNoSym,InToolbar) 0
set macrosdata(PNGNoSym,ModificationDate) {2017-04-08 16:21:41}

proc PNGNoSym {res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    GiD_Process Mescape Results ContOptions NumberOfColor 1 escape
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Black
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor White
    GiD_Process Mescape Results ContOptions ColorRamp Tangent
    GiD_Process Mescape Results ContOptions SetMaxOptions OutMaxColor White
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape View ReadView "${view_file}"
    GiD_Process Mescape results contoptions setmaxoptions setvalue 0
    GiD_Process Mescape results contoptions setminoptions setvalue -1e-35
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape Results ContourFill phi
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
    GiD_Process Mescape 'Hardcopy PNG "${png_file}"
    GiD_Process Mescape quit No escape
}

#[_ ""]
set macrosdata(PNGSingleSym_PG,Description) {}
set macrosdata(PNGSingleSym_PG,IsStandard) 0
set macrosdata(PNGSingleSym_PG,CreationDate) {2017-03-10 14:32:10}
set macrosdata(PNGSingleSym_PG,InToolbar) 0
set macrosdata(PNGSingleSym_PG,ModificationDate) {2017-04-08 17:06:01}
set macrosdata(PNGSingleSym_PG,PrePost) prepost
set macrosdata(PNGSingleSym_PG,Icon) {}
set macrosdata(PNGSingleSym_PG,Accelerators) {}
set macrosdata(PNGSingleSym_PG,Active) 1
set macrosdata(PNGSingleSym_PG,Number) 4
set macrosdata(PNGSingleSym_PG,Group) {}

proc PNGSingleSym_PG {res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    GiD_Process Mescape Utilities transformation create yes Mirror FJoin $corner2 FJoin $corner1 TwoDim
    GiD_Process Mescape Results ContOptions NumberOfColor 50
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Standard 
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor Standard 
    GiD_Process Mescape Results ContOptions ColorRamp Tangent 
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor White 
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor Black 
    GiD_Process Mescape Results ContOptions ColorRamp Tangent 
    GiD_Process Mescape View ReadView "${view_file}"
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape Results ContourFill {Gamma FUNCTION NODAL }
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
    GiD_Process Mescape 'Hardcopy PNG "${png_file}"
    GiD_Process Mescape quit No escape
}

#[_ ""]
set macrosdata(PNGSingleSym,ModificationDate) {2017-04-08 17:05:18}
set macrosdata(PNGSingleSym,Description) {}
set macrosdata(PNGSingleSym,IsStandard) 0
set macrosdata(PNGSingleSym,Accelerators) {}
set macrosdata(PNGSingleSym,Active) 1
set macrosdata(PNGSingleSym,Number) 5
set macrosdata(PNGSingleSym,InToolbar) 0
set macrosdata(PNGSingleSym,PrePost) prepost
set macrosdata(PNGSingleSym,Group) {}
set macrosdata(PNGSingleSym,CreationDate) {2017-01-25 21:00:12}
set macrosdata(PNGSingleSym,Icon) {}

proc PNGSingleSym {res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    GiD_Process Mescape Utilities transformation create yes Mirror FJoin $corner2 FJoin $corner1 TwoDim
    GiD_Process Mescape Results ContOptions NumberOfColor 1 escape
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Black
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor White
    GiD_Process Mescape Results ContOptions ColorRamp Tangent
    GiD_Process Mescape Results ContOptions SetMaxOptions OutMaxColor White
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape View ReadView "${view_file}"
    GiD_Process Mescape results contoptions setmaxoptions setvalue 0
    GiD_Process Mescape results contoptions setminoptions setvalue -1e-35
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape Results ContourFill phi
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
    GiD_Process Mescape 'Hardcopy PNG "${png_file}"
    GiD_Process Mescape quit No escape
}

#[_ ""]
set macrosdata(PNGDoubleSym_PG,Description) {}
set macrosdata(PNGDoubleSym_PG,Active) 1
set macrosdata(PNGDoubleSym_PG,IsStandard) 0
set macrosdata(PNGDoubleSym_PG,ModificationDate) {2017-04-08 17:43:20}
set macrosdata(PNGDoubleSym_PG,CreationDate) {2017-01-26 13:49:33}
set macrosdata(PNGDoubleSym_PG,Group) {}
set macrosdata(PNGDoubleSym_PG,Number) 6
set macrosdata(PNGDoubleSym_PG,Accelerators) {}
set macrosdata(PNGDoubleSym_PG,Icon) {}
set macrosdata(PNGDoubleSym_PG,PrePost) prepost
set macrosdata(PNGDoubleSym_PG,InToolbar) 0

proc PNGDoubleSym_PG {res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    if {$symcase} {
    GiD_Process Mescape Utilities transformation create yes Mirror FJoin $corner2 FJoin $corner1 TwoDim
    GiD_Process Mescape Utilities transformation create yes Mirror FJoin $corner4  FJoin $corner1 TwoDim
    } else  {
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner2 FJoin $corner1 
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner4  FJoin $corner1
    }
    GiD_Process Mescape Results ContOptions NumberOfColor 50
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Standard 
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor Standard 
    GiD_Process Mescape Results ContOptions ColorRamp Tangent 
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor White 
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor Black 
    GiD_Process Mescape Results ContOptions ColorRamp Tangent 
    GiD_Process Mescape View ReadView "${view_file}"
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape Results ContourFill {Gamma FUNCTION NODAL }
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
    GiD_Process Mescape 'Hardcopy PNG "${png_file}"
    GiD_Process Mescape quit No escape
}

#[_ ""]
set macrosdata(PNGDoubleSym,Group) {}
set macrosdata(PNGDoubleSym,PrePost) prepost
set macrosdata(PNGDoubleSym,Active) 1
set macrosdata(PNGDoubleSym,Icon) {}
set macrosdata(PNGDoubleSym,Description) {}
set macrosdata(PNGDoubleSym,InToolbar) 0
set macrosdata(PNGDoubleSym,Number) 7
set macrosdata(PNGDoubleSym,IsStandard) 0
set macrosdata(PNGDoubleSym,CreationDate) {2017-02-24 15:17:10}
set macrosdata(PNGDoubleSym,ModificationDate) {2017-04-08 17:43:13}
set macrosdata(PNGDoubleSym,Accelerators) {}

proc PNGDoubleSym {res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    if {$symcase} {
    GiD_Process Mescape Utilities transformation create yes Mirror FJoin $corner2 FJoin $corner1 TwoDim
    GiD_Process Mescape Utilities transformation create yes Mirror FJoin $corner4  FJoin $corner1 TwoDim
    } else  {
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner2 FJoin $corner1 
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner4  FJoin $corner1
    }
    GiD_Process Mescape Results ContOptions NumberOfColor 1 escape
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Black
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor White
    GiD_Process Mescape Results ContOptions ColorRamp Tangent
    GiD_Process Mescape Results ContOptions SetMaxOptions OutMaxColor White
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape View ReadView "${view_file}"
    GiD_Process Mescape results contoptions setmaxoptions setvalue 0
    GiD_Process Mescape results contoptions setminoptions setvalue -1e-35
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape Results ContourFill phi
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
    GiD_Process Mescape 'Hardcopy PNG "${png_file}"
    GiD_Process Mescape quit No escape
}

#[_ ""]
set macrosdata(PNGTripleSym_PG,IsStandard) 0
set macrosdata(PNGTripleSym_PG,Icon) {}
set macrosdata(PNGTripleSym_PG,Group) {}
set macrosdata(PNGTripleSym_PG,Active) 1
set macrosdata(PNGTripleSym_PG,PrePost) prepost
set macrosdata(PNGTripleSym_PG,InToolbar) 0
set macrosdata(PNGTripleSym_PG,CreationDate) {2017-03-10 14:20:00}
set macrosdata(PNGTripleSym_PG,Number) 8
set macrosdata(PNGTripleSym_PG,ModificationDate) {2017-04-08 17:42:03}
set macrosdata(PNGTripleSym_PG,Description) {}
set macrosdata(PNGTripleSym_PG,Accelerators) {}

proc PNGTripleSym_PG {res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner2 FJoin $corner6
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner3 FJoin $corner1
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner4 FJoin $corner2
    GiD_Process Mescape Results ContOptions NumberOfColor 50
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Standard 
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor Standard 
    GiD_Process Mescape Results ContOptions ColorRamp Tangent 
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor White 
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor Black 
    GiD_Process Mescape Results ContOptions ColorRamp Tangent 
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape View ReadView "${view_file}"
    GiD_Process Mescape Results ContourFill {Gamma FUNCTION NODAL }
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
    GiD_Process Mescape 'Hardcopy PNG "${png_file}"
    GiD_Process Mescape quit No escape
}

#[_ ""]
set macrosdata(PNGTripleSym,CreationDate) {2017-03-10 12:11:55}
set macrosdata(PNGTripleSym,Number) 9
set macrosdata(PNGTripleSym,IsStandard) 0
set macrosdata(PNGTripleSym,Description) {}
set macrosdata(PNGTripleSym,Accelerators) {}
set macrosdata(PNGTripleSym,Group) {}
set macrosdata(PNGTripleSym,InToolbar) 0
set macrosdata(PNGTripleSym,Icon) {}
set macrosdata(PNGTripleSym,PrePost) prepost
set macrosdata(PNGTripleSym,Active) 1
set macrosdata(PNGTripleSym,ModificationDate) {2017-04-08 17:41:44}

proc PNGTripleSym {res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner2 FJoin $corner6
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner3 FJoin $corner1
    GiD_Process Mescape Utilities transformation create yes Translation FJoin $corner4 FJoin $corner2
    GiD_Process Mescape Results ContOptions NumberOfColor 1 escape
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Black
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor White
    GiD_Process Mescape Results ContOptions ColorRamp Tangent
    GiD_Process Mescape Results ContOptions SetMaxOptions OutMaxColor White
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape View ReadView "${view_file}"
    GiD_Process Mescape results contoptions setmaxoptions setvalue 0
    GiD_Process Mescape results contoptions setminoptions setvalue -1e-35
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape Results ContourFill phi
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
    GiD_Process Mescape 'Hardcopy PNG "${png_file}"
    GiD_Process Mescape quit No escape
}

#[_ ""]
set macrosdata(GeomSTL,PrePost) prepost
set macrosdata(GeomSTL,CreationDate) {2017-01-20 13:02:17}
set macrosdata(GeomSTL,InToolbar) 0
set macrosdata(GeomSTL,Icon) {}
set macrosdata(GeomSTL,ModificationDate) {2017-04-08 17:51:17}
set macrosdata(GeomSTL,Active) 1
set macrosdata(GeomSTL,Accelerators) {}
set macrosdata(GeomSTL,IsStandard) 0
set macrosdata(GeomSTL,Description) {}
set macrosdata(GeomSTL,Group) {}
set macrosdata(GeomSTL,Number) 10

proc GeomSTL {res_file ascii_file geom_file} {
    GiD_Process Mescape Postprocess
    GiD_Process Mescape Files Read -rebuildIndex:0 ${res_file}
    GiD_Process Mescape Results IsoSurfaces DisplayStyle ContourFillColor 
    GiD_Process Mescape Results IsoSurfaces Exact phi 1 0.0
    GiD_Process Mescape 'Redraw
    GiD_Process Mescape Results isosurfaces TurnIntoCuts Mescape DoCut TurnIntoMeshes
    GiD_Process Mescape Results Options ExtractBoundaries
    GiD_Process Mescape Utilities Delete SurfaceSets {WORKPIECE 1} Yes
    GiD_Process Mescape Utilities Delete CutSets {IsoSurf phi = 0 WORKPIECE 1} Yes
    GiD_Process Mescape files saveall allmeshessets ${ascii_file}
    GiD_Process Mescape Preprocess Yes
    GiD_Process Mescape Files MeshRead "${ascii_file}.msh"
    GiD_Process Mescape Files WriteForBAS {C:/Program Files/GiD/GiD 13.0.1/templates/DXF.bas} ${geom_file}
    GiD_Process Mescape Files DxfRead -collapse:1 -ignorelayer:1 -tolerance:0.0001 -autotol:1 -- ${geom_file}
    GiD_Process Mescape Geometry Create IntMultLines InvertSelection escape Yes
    GiD_Process Mescape Geometry Delete AllTypes
}

#[_ ""]
set macrosdata(ExportSTL,Group) {}
set macrosdata(ExportSTL,Active) 1
set macrosdata(ExportSTL,CreationDate) {2017-01-20 13:33:29}
set macrosdata(ExportSTL,PrePost) prepost
set macrosdata(ExportSTL,Number) 11
set macrosdata(ExportSTL,Accelerators) {}
set macrosdata(ExportSTL,Description) {}
set macrosdata(ExportSTL,InToolbar) 1
set macrosdata(ExportSTL,ModificationDate) {2017-03-03 16:27:51}
set macrosdata(ExportSTL,IsStandard) 0
set macrosdata(ExportSTL,Icon) {volumesets_OneColor.png imported_images volumesets_OneColor.png themed_image}

proc ExportSTL {} {
    GiD_Process Mescape Utilities Copy Surfaces DoExtrude Surfaces MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,0.16 InvertSelection escape
    GiD_Process Mescape Meshing Generate Yes 0.11 MeshingParametersFrom=Preferences
    GiD_Process Mescape Files WriteForBAS {C:/Program Files/GiD/GiD 13.0.1/templates/STL.bas} {C:/Users/ferry/Downloads/result.stl}
    GiD_Process Mescape Quit No
}

#[_ ""]
set macrosdata(CreateSurfaceSTL,Active) 1
set macrosdata(CreateSurfaceSTL,Group) {}
set macrosdata(CreateSurfaceSTL,Accelerators) {}
set macrosdata(CreateSurfaceSTL,IsStandard) 0
set macrosdata(CreateSurfaceSTL,Number) 12
set macrosdata(CreateSurfaceSTL,ModificationDate) {2017-03-03 16:37:28}
set macrosdata(CreateSurfaceSTL,InToolbar) 1
set macrosdata(CreateSurfaceSTL,Description) {}
set macrosdata(CreateSurfaceSTL,PrePost) prepost
set macrosdata(CreateSurfaceSTL,CreationDate) {2017-02-24 15:48:38}
set macrosdata(CreateSurfaceSTL,Icon) {surface.png imported_images surface.png themed_image}

proc CreateSurfaceSTL {} {
    GiD_Process Mescape Geometry Create NurbsSurface InvertSelection Mescape
}

#[_ ""]
set macrosdata(ShowTopologyLevelSet,IsStandard) 0
set macrosdata(ShowTopologyLevelSet,Group) {}
set macrosdata(ShowTopologyLevelSet,Active) 1
set macrosdata(ShowTopologyLevelSet,Icon) {surfacesets_on_OneColor.png imported_images surfacesets_on_OneColor.png themed_image}
set macrosdata(ShowTopologyLevelSet,Description) {}
set macrosdata(ShowTopologyLevelSet,Number) 13
set macrosdata(ShowTopologyLevelSet,CreationDate) {2016-12-23 13:09:50}
set macrosdata(ShowTopologyLevelSet,InToolbar) 1
set macrosdata(ShowTopologyLevelSet,PrePost) prepost
set macrosdata(ShowTopologyLevelSet,ModificationDate) {2017-04-07 23:33:40}
set macrosdata(ShowTopologyLevelSet,Accelerators) {}

proc ShowTopologyLevelSet {} {
    GiD_Process Mescape Results ContOptions NumberOfColor 1 escape
    GiD_Process Mescape Results ContOptions SetMinOptions MinColor Black
    GiD_Process Mescape Results ContOptions SetMaxOptions MaxColor White
    GiD_Process Mescape Results ContOptions ColorRamp Tangent
    GiD_Process Mescape Results ContOptions SetMaxOptions OutMaxColor White
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape results contoptions setmaxoptions setvalue 0
    GiD_Process Mescape results contoptions setminoptions setvalue -1e-35
    GiD_Process Mescape Utilities Redisplay
    GiD_Process Mescape Results ContourFill phi
    GiD_Process Mescape 'Hardcopy Options ShowLegends No
    GiD_Process Mescape 'Hardcopy Options ShowAxes No
    GiD_Process Mescape 'Hardcopy Options PrintLogo No
}

#[_ "Open a simple editor to save small notes with the model"]
set macrosdata(Notes,Active) 1
set macrosdata(Notes,PrePost) prepost
set macrosdata(Notes,CreationDate) {2006-01-20 17:17:27}
set macrosdata(Notes,Number) 14
set macrosdata(Notes,Accelerators) {}
set macrosdata(Notes,InToolbar) 1
set macrosdata(Notes,Description) {Open a simple editor to save small notes with the model}
set macrosdata(Notes,Group) Tools
set macrosdata(Notes,Icon) {note.png imported_images note.png themed_image}
set macrosdata(Notes,IsStandard) 1
set macrosdata(Notes,ModificationDate) {2011-11-04 11:03:51}

proc Notes {} {
    OpenNotes
}

#[_ "Resize Graphical Window"]
set {macrosdata(Resize Graphical Window,Accelerators)} {}
set {macrosdata(Resize Graphical Window,Group)} View
set {macrosdata(Resize Graphical Window,InToolbar)} 1
set {macrosdata(Resize Graphical Window,Active)} 1
set {macrosdata(Resize Graphical Window,PrePost)} prepost
set {macrosdata(Resize Graphical Window,Description)} {Resize Graphical Window}
set {macrosdata(Resize Graphical Window,Icon)} {PantallaResize.png imported_images PantallaResize.png themed_image}
set {macrosdata(Resize Graphical Window,CreationDate)} {2010-02-12 13:01:02}
set {macrosdata(Resize Graphical Window,IsStandard)} 1
set {macrosdata(Resize Graphical Window,Number)} 15
set {macrosdata(Resize Graphical Window,ModificationDate)} {2011-11-04 11:03:51}

proc {Resize Graphical Window} {} {
    set info [ wm geometry .gid]
    set info2 [.central.s configure]
    set height_list [ lindex $info2 0]
    set width_list [ lindex $info2 1]
    if { [ regexp {([0-9]*)x([0-9]*)} $info trozo top_width top_height] } { 
        set gid_width [ lindex $width_list 4]
        set gid_height [ lindex $height_list 4]
        set width_offset [ expr $top_width - [ lindex $width_list 4]]
        set height_offset [ expr $top_height - [ lindex $height_list 4]]
        
        set value [ tk_dialogEntryRAM .gid.wResize "get new size"  "Enter new size of the graphical window ( format WIDTHxHEIGHT, for instance 800x600)"  gidquestionhead word ${gid_width}x${gid_height}]
        
        if { "$value" != "--CANCEL--"} {
            if { [ regexp {([0-9]*)x([0-9]*)} $value  trozo new_width new_height] } {
                wm geometry .gid  [ expr $width_offset + $new_width]x[ expr $height_offset + $new_height]
            }
        }
    }
}

#[_ "Dimension angle beween three points"]
set {macrosdata(Dimension angle,CreationDate)} {2006-02-14 17:36:29}
set {macrosdata(Dimension angle,Description)} {Dimension angle beween three points}
set {macrosdata(Dimension angle,Number)} 16
set {macrosdata(Dimension angle,PrePost)} prepost
set {macrosdata(Dimension angle,InToolbar)} 0
set {macrosdata(Dimension angle,Accelerators)} {}
set {macrosdata(Dimension angle,Group)} Dimension
set {macrosdata(Dimension angle,Icon)} {dimension_angle.png imported_images dimension_angle.png themed_image}
set {macrosdata(Dimension angle,IsStandard)} 1
set {macrosdata(Dimension angle,Active)} 1
set {macrosdata(Dimension angle,ModificationDate)} {2011-11-04 11:03:51}

proc {Dimension angle} {} {
    GiD_Process escape escape escape escape Utilities Dimension Create Angle
}

#[_ "Dimension pointing to a vertex"]
set {macrosdata(Dimension vertex,Active)} 1
set {macrosdata(Dimension vertex,InToolbar)} 0
set {macrosdata(Dimension vertex,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Dimension vertex,Icon)} {dimension_vertex.png imported_images dimension_vertex.png themed_image}
set {macrosdata(Dimension vertex,Number)} 17
set {macrosdata(Dimension vertex,Description)} {Dimension pointing to a vertex}
set {macrosdata(Dimension vertex,IsStandard)} 1
set {macrosdata(Dimension vertex,PrePost)} prepost
set {macrosdata(Dimension vertex,CreationDate)} {2006-02-14 17:51:53}
set {macrosdata(Dimension vertex,Group)} Dimension
set {macrosdata(Dimension vertex,Accelerators)} {}

proc {Dimension vertex} {} {
    GiD_Process MEscape Utilities Dimension Create Vertex
}

#[_ "Dimension pointing to an arc center"]
set {macrosdata(Dimension arc,Number)} 18
set {macrosdata(Dimension arc,IsStandard)} 1
set {macrosdata(Dimension arc,CreationDate)} {2006-02-14 21:00:15}
set {macrosdata(Dimension arc,InToolbar)} 0
set {macrosdata(Dimension arc,Accelerators)} {}
set {macrosdata(Dimension arc,Description)} {Dimension pointing to an arc center}
set {macrosdata(Dimension arc,PrePost)} prepost
set {macrosdata(Dimension arc,Active)} 1
set {macrosdata(Dimension arc,Group)} Dimension
set {macrosdata(Dimension arc,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Dimension arc,Icon)} {dimension_arc.png imported_images dimension_arc.png themed_image}

proc {Dimension arc} {} {
    GiD_Process MEscape Utilities Dimension Create Arc
}

#[_ "Modify the dimension text"]
set {macrosdata(Dimension edit,InToolbar)} 1
set {macrosdata(Dimension edit,Accelerators)} {}
set {macrosdata(Dimension edit,Description)} {Modify the dimension text}
set {macrosdata(Dimension edit,PrePost)} prepost
set {macrosdata(Dimension edit,Active)} 1
set {macrosdata(Dimension edit,Group)} Dimension
set {macrosdata(Dimension edit,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Dimension edit,Icon)} {dimension_edit.png imported_images dimension_edit.png themed_image}
set {macrosdata(Dimension edit,Number)} 19
set {macrosdata(Dimension edit,IsStandard)} 1
set {macrosdata(Dimension edit,CreationDate)} {2004-11-18 21:53:33}

proc {Dimension edit} {} {
    GiD_Process MEscape Utilities Dimension Edit
}

#[_ "Delete a dimension"]
set {macrosdata(Dimension delete,Number)} 20
set {macrosdata(Dimension delete,Accelerators)} {}
set {macrosdata(Dimension delete,InToolbar)} 0
set {macrosdata(Dimension delete,Description)} {Delete a dimension}
set {macrosdata(Dimension delete,Group)} Dimension
set {macrosdata(Dimension delete,Icon)} {dimension_delete.png imported_images dimension_delete.png themed_image}
set {macrosdata(Dimension delete,IsStandard)} 1
set {macrosdata(Dimension delete,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Dimension delete,Active)} 1
set {macrosdata(Dimension delete,PrePost)} prepost
set {macrosdata(Dimension delete,CreationDate)} {2004-11-18 21:53:33}

proc {Dimension delete} {} {
    GiD_Process MEscape Utilities Dimension Delete
}

#[_ "Separate element in layers by its material"]
set macrosdata(CreateLayersFromElementMaterials,Group) Layers
set macrosdata(CreateLayersFromElementMaterials,Number) 21
set macrosdata(CreateLayersFromElementMaterials,Description) {Separate element in layers by its material}
set macrosdata(CreateLayersFromElementMaterials,IsStandard) 1
set macrosdata(CreateLayersFromElementMaterials,ModificationDate) {2012-12-14 22:44:00}
set macrosdata(CreateLayersFromElementMaterials,CreationDate) {2012-12-14 22:39:32}
set macrosdata(CreateLayersFromElementMaterials,Active) 1
set macrosdata(CreateLayersFromElementMaterials,Icon) {sendto.png imported_images sendto.png themed_image}
set macrosdata(CreateLayersFromElementMaterials,Accelerators) {}
set macrosdata(CreateLayersFromElementMaterials,InToolbar) 1
set macrosdata(CreateLayersFromElementMaterials,PrePost) pre

proc CreateLayersFromElementMaterials {} {
    proc _CreateLayersFromElementMaterials { } {
        set used_materials [list]
        foreach mesh [GiD_Info mesh elements any -array] {
            lappend used_materials {*}[lsort -integer -unique [lindex $mesh 3]]        
        }
        set used_materials [lsort -integer -unique $used_materials]
        if { [llength $used_materials] } {
            if { [GiD_Info project ViewMode] != "MESHUSE" } {
                GiD_Process Mescape Meshing MeshView
            }
            GidUtils::DisableGraphics
            set current_layers [GiD_Info Layers]
            foreach material_id $used_materials {
                set layer_name Mat_$material_id
                if { [lsearch $current_layers $layer_name] == -1 } {
                    GiD_Process 'Layer New $layer_name escape
                }
                set element_ids [list]
                foreach mesh [GiD_Info mesh elements any -array] {
                    foreach element_id [lindex $mesh 1] element_material [lindex $mesh 3] {
                        if { $material_id == $element_material } {
                            lappend element_ids $element_id
                        }
                    }                
                }
                GiD_Process 'Layer Entities $layer_name Elements {*}$element_ids escape
            } 
            GidUtils::EnableGraphics
        }
        FillInfoLayers
    }
    
    _CreateLayersFromElementMaterials
}

#[_ "Separate each volume in a different layer"]
set macrosdata(CreateLayersFromGeometryVolumes,Icon) {sendto.png imported_images sendto.png themed_image}
set macrosdata(CreateLayersFromGeometryVolumes,InToolbar) 1
set macrosdata(CreateLayersFromGeometryVolumes,Active) 1
set macrosdata(CreateLayersFromGeometryVolumes,CreationDate) {2013-10-30 16:53:57}
set macrosdata(CreateLayersFromGeometryVolumes,PrePost) pre
set macrosdata(CreateLayersFromGeometryVolumes,Group) Layers
set macrosdata(CreateLayersFromGeometryVolumes,ModificationDate) {2013-10-30 16:55:53}
set macrosdata(CreateLayersFromGeometryVolumes,IsStandard) 1
set macrosdata(CreateLayersFromGeometryVolumes,Accelerators) {}
set macrosdata(CreateLayersFromGeometryVolumes,Number) 22
set macrosdata(CreateLayersFromGeometryVolumes,Description) {Separate each volume in a different layer}

proc CreateLayersFromGeometryVolumes {} {
    proc _CreateLayersFromGeometryVolumes { } {   
        if { [GiD_Info Geometry NumVolumes]>1 } {
            GidUtils::DisableGraphics
            set current_layers [GiD_Info Layers]
            foreach volume_id [GiD_Geometry list volume 1:end] {
                set layer_name Vol_$volume_id
                if { [lsearch $current_layers $layer_name] == -1 } {
                    GiD_Process 'Layer New $layer_name escape
                } else {
                    GiD_Process 'Layer ToUse $layer_name escape
                }
                GiD_Process 'Layer Entities $layer_name LowerEntities Volumes $volume_id escape
            } 
            GidUtils::EnableGraphics
            FillInfoLayers
        }    
    }
    
    _CreateLayersFromGeometryVolumes
}

#[_ "Import a selection of multiple IGES files and send then to layers named with file names."]
set macrosdata(IgesFilesToLayers,IsStandard) 1
set macrosdata(IgesFilesToLayers,Active) 1
set macrosdata(IgesFilesToLayers,CreationDate) {2013-01-23 12:17:04}
set macrosdata(IgesFilesToLayers,Group) {Geometry Creation}
set macrosdata(IgesFilesToLayers,Accelerators) {}
set macrosdata(IgesFilesToLayers,Number) 23
set macrosdata(IgesFilesToLayers,ModificationDate) {2013-01-23 12:36:33}
set macrosdata(IgesFilesToLayers,Description) {Import a selection of multiple IGES files and send then to layers named with file names.}
set macrosdata(IgesFilesToLayers,PrePost) pre
set macrosdata(IgesFilesToLayers,InToolbar) 1
set macrosdata(IgesFilesToLayers,Icon) {folder.png imported_images folder.png themed_image}

proc IgesFilesToLayers {} {
    proc IgesFilesToLayers { } {
        set multiple 1
        set filenames [Browser-ramR file read .gid [_ "Read IGES files"]  {} {{{IGES} {.igs }} {{All files} {.*}}} {} $multiple  ""]
        if { $filenames != "" } {
            GidUtils::WaitState .gid
            GidUtils::DisableGraphics
            set old_value_IGES_CreateAllInLayerToUse [GiD_Set IGES(CreateAllInLayerToUse)]
            GiD_Set IGES(CreateAllInLayerToUse) 1
            set auto_collapse [GiD_Set AutoCollapseAfterImport]
            if { $auto_collapse } {
                GiD_Set AutoCollapseAfterImport 0 ;#temporary disable it to do the collapse only once
                set old_value_CollapseIgnoringLayers [GiD_Set CollapseIgnoringLayers]
                GiD_Set CollapseIgnoringLayers 0 ;#to not join entities of different files
            }
            foreach filename $filenames {
                set layername [file tail $filename]
                if { [lsearch [GiD_Info layers] $layername] == -1 } {
                    GiD_Process 'Layers New $layername
                } else {
                    GiD_Process 'Layers ToUse $layername
                }
                GiD_Process Mescape Files IgesRead $filename
            }
            GiD_Process 'Layers Delete Layer0 ;#delete this layer if empty
            if { $auto_collapse } {
                GiD_Process Mescape Utilities Collapse Model Yes
                GiD_Set CollapseIgnoringLayers $old_value_CollapseIgnoringLayers
            }
            GiD_Set IGES(CreateAllInLayerToUse) $old_value_IGES_CreateAllInLayerToUse
            GidUtils::EnableGraphics
            GidUtils::EndWaitState .gid
            FillInfoLayers
            GiD_Process 'Zoom Frame
            GidUtils::SetWarnLine [_ "Read %s files" [llength $filenames]]
        }
        return 0
    }
    
    IgesFilesToLayers
}

#[_ "Import a selection of multiple Wavefront Object files and send then to the current working layer."]
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,Group) {Mesh Creation}
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,Active) 1
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,CreationDate) {2013-01-23 12:17:04}
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,PrePost) pre
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,Number) 24
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,Accelerators) {}
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,Description) {Import a selection of multiple Wavefront Object files and send then to the current working layer.}
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,InToolbar) 0
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,ModificationDate) {2013-01-23 12:36:33}
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,IsStandard) 1
set macrosdata(MultipleWavefrontObjectsToCurrentLayer,Icon) {folder.png imported_images folder.png themed_image}

proc MultipleWavefrontObjectsToCurrentLayer {} {
    proc MultipleWavefrontObjectsToCurrentLayer { } {
        set multiple 1
        set filenames [Browser-ramR file read .gid [_ "Read OBJ files"]  {} {{{Wavefront Object mesh files} {.obj }} {{All files} {.*}}} {} $multiple  ""]
        if { $filenames != "" } {
            set num_files [ llength $filenames]
            if { $num_files > 10} {
                # ask for confirmation
                # set resp [ tk_messageBox -parent .gid -title [_ "Confirm"]  #                -type yesno -default no -icon question  #                -message [_ "%d Wavefront Object files are to be imported. Continue ?" $num_files]]
                set resp [ tk_dialogRAM .gid.___confirm_obj_import [_ "Confirm"]  [_ "%d Wavefront Object files are to be imported. Continue ?" $num_files]  gidquestionhead 1 [_ "Yes"] [_ "No#C#I don't want to do that"]]
                if { $resp == 1} {
                    set resp "no"
                }
                if { [ string equal $resp "no"]} {
                    return
                }
            }
            GidUtils::WaitState .gid
            GidUtils::DisableGraphics
            set auto_collapse [GiD_Set AutoCollapseAfterImport]
            if { $auto_collapse } {
                GiD_Set AutoCollapseAfterImport 0 ;#temporary disable it to do the collapse only once
                set old_value_CollapseIgnoringLayers [GiD_Set CollapseIgnoringLayers]
                GiD_Set CollapseIgnoringLayers 0 ;#to not join entities of different files
            }
            # add collapsing
            set end_progress [ expr $num_files + 1]
            # GidUtils::EnableGraphics
            PostProgressBar 0 $end_progress [_ "Importing %d Wavefront Object files" $num_files]  [_ "Importing %d Wavefront Object files" $num_files]
            # update idletasks
            # GidUtils::DisableGraphics
            set idx_files 0
            foreach filename $filenames {
                # set layername [file tail $filename]
                # if { [lsearch [GiD_Info layers] $layername] == -1 } {
                    #     GiD_Process 'Layers New $layername
                    # } else {
                    #     GiD_Process 'Layers ToUse $layername
                    # }
                GiD_Process Mescape Files OBJRead $filename
                incr idx_files
                # GidUtils::EnableGraphics
                PostProgressBar $idx_files $end_progress [_ "Imported %d of %d Wavefront Objects" $idx_files $num_files]
                # update idletasks
                # GidUtils::DisableGraphics
            }
            # GiD_Process 'Layers Delete Layer0 ;#delete this layer if empty
            if { $auto_collapse } {
                # GidUtils::EnableGraphics
                PostProgressBar $num_files $end_progress [_ "Collapsing imported objects"]
                # update idletasks
                # GidUtils::DisableGraphics
                GiD_Process Mescape Utilities Collapse Model Yes
                GiD_Set AutoCollapseAfterImport $auto_collapse
                GiD_Set CollapseIgnoringLayers $old_value_CollapseIgnoringLayers
            }
            PostProgressBar $end_progress $end_progress
            GidUtils::EnableGraphics
            GidUtils::EndWaitState .gid
            FillInfoLayers
            GiD_Process 'Zoom Frame
            if { [GiD_Info project ViewMode] != "MESHUSE" } {
                GiD_Process escape escape escape escape Meshing MeshView escape escape escape escape escape
            }
            GidUtils::SetWarnLine [_ "Read %s Wavefront OBJect files" [llength $filenames]]
        }
        return 0
    }
    
    MultipleWavefrontObjectsToCurrentLayer
}

#[_ "Send entities to the back (hidden) part of a layer"]
set {macrosdata(Send to back,Icon)} {send_to_back.png imported_images send_to_back.png themed_image}
set {macrosdata(Send to back,Description)} {Send entities to the back (hidden) part of a layer}
set {macrosdata(Send to back,Active)} 1
set {macrosdata(Send to back,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Send to back,Group)} Layers
set {macrosdata(Send to back,InToolbar)} 0
set {macrosdata(Send to back,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Send to back,Number)} 25
set {macrosdata(Send to back,IsStandard)} 1
set {macrosdata(Send to back,PrePost)} pre
set {macrosdata(Send to back,Accelerators)} Ctrl-F1

proc {Send to back} {} {
    GiD_Process 'Layers SendToBack LowerEntities Surfaces
}

#[_ "Send unselected entities to the back (hidden) part of a layer"]
set {macrosdata(Send to back (Opposite),PrePost)} pre
set {macrosdata(Send to back (Opposite),Accelerators)} Ctrl-F2
set {macrosdata(Send to back (Opposite),Active)} 1
set {macrosdata(Send to back (Opposite),Description)} {Send unselected entities to the back (hidden) part of a layer}
set {macrosdata(Send to back (Opposite),IsStandard)} 1
set {macrosdata(Send to back (Opposite),ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Send to back (Opposite),Number)} 26
set {macrosdata(Send to back (Opposite),Icon)} {send_to_back_opposite.png imported_images send_to_back_opposite.png themed_image}
set {macrosdata(Send to back (Opposite),InToolbar)} 0
set {macrosdata(Send to back (Opposite),Group)} Layers
set {macrosdata(Send to back (Opposite),CreationDate)} {2004-11-18 21:53:33}

proc {Send to back (Opposite)} {} {
    GiD_Process 'Layers SendToBackOpposite LowerEntities Surfaces
}

#[_ "Send to the visible front part entities in back layers"]
set {macrosdata(Bring to front,IsStandard)} 1
set {macrosdata(Bring to front,Description)} {Send to the visible front part entities in back layers}
set {macrosdata(Bring to front,Active)} 1
set {macrosdata(Bring to front,PrePost)} pre
set {macrosdata(Bring to front,Group)} Layers
set {macrosdata(Bring to front,CreationDate)} {2004-11-30 18:57:31}
set {macrosdata(Bring to front,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Bring to front,Icon)} {bring_to_front.png imported_images bring_to_front.png themed_image}
set {macrosdata(Bring to front,Number)} 27
set {macrosdata(Bring to front,Accelerators)} {}
set {macrosdata(Bring to front,InToolbar)} 0

proc {Bring to front} {} {
    GiD_Process 'Layers BringToFrontAll Escape
}

#[_ "Send to back"]
set {macrosdata(Send to back E,Group)} Layers
set {macrosdata(Send to back E,Number)} 28
set {macrosdata(Send to back E,PrePost)} post
set {macrosdata(Send to back E,CreationDate)} {2004-11-30 20:20:29}
set {macrosdata(Send to back E,InToolbar)} 0
set {macrosdata(Send to back E,Icon)} {send_to_back.png imported_images send_to_back.png themed_image}
set {macrosdata(Send to back E,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Send to back E,Active)} 0
set {macrosdata(Send to back E,Accelerators)} Ctrl-F1
set {macrosdata(Send to back E,IsStandard)} 1
set {macrosdata(Send to back E,Description)} {Send to back}

proc {Send to back E} {} {
    set ::GidPriv(LayersBackOpposite) 0
    PDSendToBack ""  Elements
}

#[_ "Send to back (Opposite)"]
set {macrosdata(Send to back E (opposite),Description)} {Send to back (Opposite)}
set {macrosdata(Send to back E (opposite),Active)} 0
set {macrosdata(Send to back E (opposite),Accelerators)} Ctrl-F2
set {macrosdata(Send to back E (opposite),Number)} 29
set {macrosdata(Send to back E (opposite),ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Send to back E (opposite),Group)} Layers
set {macrosdata(Send to back E (opposite),Icon)} {send_to_back_opposite.png imported_images send_to_back_opposite.png themed_image}
set {macrosdata(Send to back E (opposite),PrePost)} post
set {macrosdata(Send to back E (opposite),InToolbar)} 0
set {macrosdata(Send to back E (opposite),CreationDate)} {2004-11-30 20:24:54}
set {macrosdata(Send to back E (opposite),IsStandard)} 1

proc {Send to back E (opposite)} {} {
    set ::GidPriv(LayersBackOpposite) 1
    PDSendToBack ""  Elements
}

#[_ "Bring to front"]
set {macrosdata(Bring to front E,CreationDate)} {2004-11-30 20:21:40}
set {macrosdata(Bring to front E,Accelerators)} {}
set {macrosdata(Bring to front E,Active)} 0
set {macrosdata(Bring to front E,Description)} {Bring to front}
set {macrosdata(Bring to front E,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Bring to front E,InToolbar)} 0
set {macrosdata(Bring to front E,Group)} Layers
set {macrosdata(Bring to front E,PrePost)} post
set {macrosdata(Bring to front E,Icon)} {bring_to_front.png imported_images bring_to_front.png themed_image}
set {macrosdata(Bring to front E,Number)} 30
set {macrosdata(Bring to front E,IsStandard)} 1

proc {Bring to front E} {} {
    PDSendToBack "" AllToFront
}

#[_ "Create a trimmed surface from a base surface and a new boundary"]
set {macrosdata(Trimmed surface,Accelerators)} {}
set {macrosdata(Trimmed surface,Active)} 1
set {macrosdata(Trimmed surface,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Trimmed surface,Number)} 31
set {macrosdata(Trimmed surface,PrePost)} pre
set {macrosdata(Trimmed surface,Description)} {Create a trimmed surface from a base surface and a new boundary}
set {macrosdata(Trimmed surface,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Trimmed surface,Group)} {Geometry Creation}
set {macrosdata(Trimmed surface,Icon)} {surface_trim.png imported_images surface_trim.png themed_image}
set {macrosdata(Trimmed surface,IsStandard)} 1
set {macrosdata(Trimmed surface,InToolbar)} 1

proc {Trimmed surface} {} {
    GiD_Process MEscape Geometry Create TrimmedNurb
}

#[_ "Create a hole in a surface from a line's loop "]
set {macrosdata(Hole surface,IsStandard)} 1
set {macrosdata(Hole surface,InToolbar)} 1
set {macrosdata(Hole surface,Description)} {Create a hole in a surface from a line's loop }
set {macrosdata(Hole surface,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Hole surface,Icon)} {surface_hole.png imported_images surface_hole.png themed_image}
set {macrosdata(Hole surface,Active)} 1
set {macrosdata(Hole surface,Accelerators)} {}
set {macrosdata(Hole surface,Group)} {Geometry Edit}
set {macrosdata(Hole surface,PrePost)} pre
set {macrosdata(Hole surface,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Hole surface,Number)} 32

proc {Hole surface} {} {
    GiD_Process MEscape Geometry Edit HoleNurb
}

#[_ "Create a new surface from a family of parallel lines"]
set {macrosdata(Surfaces from parallel lines,PrePost)} pre
set {macrosdata(Surfaces from parallel lines,Description)} {Create a new surface from a family of parallel lines}
set {macrosdata(Surfaces from parallel lines,Icon)} {surface_parallel_lines.png imported_images surface_parallel_lines.png themed_image}
set {macrosdata(Surfaces from parallel lines,Group)} {Geometry Creation}
set {macrosdata(Surfaces from parallel lines,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Surfaces from parallel lines,Active)} 1
set {macrosdata(Surfaces from parallel lines,CreationDate)} {2006-02-14 21:44:20}
set {macrosdata(Surfaces from parallel lines,Number)} 33
set {macrosdata(Surfaces from parallel lines,InToolbar)} 1
set {macrosdata(Surfaces from parallel lines,Accelerators)} {}
set {macrosdata(Surfaces from parallel lines,IsStandard)} 1

proc {Surfaces from parallel lines} {} {
    GiD_Process MEscape Geometry Create NurbsSurface ParallelLines
}

#[_ "Move a point (Warning: surfaces can't support big deformation)"]
set {macrosdata(move point,Number)} 34
set {macrosdata(move point,InToolbar)} 1
set {macrosdata(move point,Icon)} {point_move.png imported_images point_move.png themed_image}
set {macrosdata(move point,PrePost)} pre
set {macrosdata(move point,Accelerators)} {}
set {macrosdata(move point,IsStandard)} 1
set {macrosdata(move point,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(move point,Group)} {Geometry Edit}
set {macrosdata(move point,Description)} {Move a point (Warning: surfaces can't support big deformation)}
set {macrosdata(move point,Active)} 1
set {macrosdata(move point,CreationDate)} {2006-02-14 22:13:38}

proc {move point} {} {
    GiD_Process MEscape Geometry Edit MovePoint
}

#[_ "Divide a line in a number of divisions"]
set {macrosdata(Divide line num divisions,Description)} {Divide a line in a number of divisions}
set {macrosdata(Divide line num divisions,InToolbar)} 1
set {macrosdata(Divide line num divisions,Active)} 1
set {macrosdata(Divide line num divisions,CreationDate)} {2006-02-14 22:32:48}
set {macrosdata(Divide line num divisions,Number)} 35
set {macrosdata(Divide line num divisions,Accelerators)} {}
set {macrosdata(Divide line num divisions,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Divide line num divisions,IsStandard)} 1
set {macrosdata(Divide line num divisions,Group)} {Geometry Edit}
set {macrosdata(Divide line num divisions,Icon)} {line_divide_nparts.png imported_images line_divide_nparts.png themed_image}
set {macrosdata(Divide line num divisions,PrePost)} pre

proc {Divide line num divisions} {} {
    GiD_Process MEscape Geometry Edit DivideLine Multiple NumDivisions
}

#[_ "Divide a line near a point"]
set {macrosdata(Divide line near point,Number)} 36
set {macrosdata(Divide line near point,Accelerators)} {}
set {macrosdata(Divide line near point,Description)} {Divide a line near a point}
set {macrosdata(Divide line near point,InToolbar)} 1
set {macrosdata(Divide line near point,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Divide line near point,IsStandard)} 1
set {macrosdata(Divide line near point,Icon)} {line_divide_nearpoint.png imported_images line_divide_nearpoint.png themed_image}
set {macrosdata(Divide line near point,Group)} {Geometry Edit}
set {macrosdata(Divide line near point,Active)} 1
set {macrosdata(Divide line near point,CreationDate)} {2006-02-14 22:22:00}
set {macrosdata(Divide line near point,PrePost)} pre

proc {Divide line near point} {} {
    GiD_Process MEscape Geometry Edit DivideLine Multiple Point PointInLine
}

#[_ " geometry create IntersectLines"]
set {macrosdata(Intersect lines,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Intersect lines,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Intersect lines,Active)} 1
set {macrosdata(Intersect lines,Icon)} {intersect_lines.png imported_images intersect_lines.png themed_image}
set {macrosdata(Intersect lines,Accelerators)} {}
set {macrosdata(Intersect lines,InToolbar)} 1
set {macrosdata(Intersect lines,PrePost)} pre
set {macrosdata(Intersect lines,Group)} {Geometry Edit}
set {macrosdata(Intersect lines,Number)} 37
set {macrosdata(Intersect lines,Description)} { geometry create IntersectLines}
set {macrosdata(Intersect lines,IsStandard)} 1

proc {Intersect lines} {} {
    GiD_Process Mescape Geometry Create IntMultLines
}

#[_ " geometry create IntersectLines NoDivideLines "]
set {macrosdata(Intersect lines (no divide),Number)} 38
set {macrosdata(Intersect lines (no divide),PrePost)} pre
set {macrosdata(Intersect lines (no divide),Accelerators)} {}
set {macrosdata(Intersect lines (no divide),Description)} { geometry create IntersectLines NoDivideLines }
set {macrosdata(Intersect lines (no divide),Group)} {Geometry Edit}
set {macrosdata(Intersect lines (no divide),IsStandard)} 1
set {macrosdata(Intersect lines (no divide),ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Intersect lines (no divide),InToolbar)} 0
set {macrosdata(Intersect lines (no divide),Active)} 1
set {macrosdata(Intersect lines (no divide),Icon)} {intersect_lines_nodivide.png imported_images intersect_lines_nodivide.png themed_image}
set {macrosdata(Intersect lines (no divide),CreationDate)} {2004-11-18 21:53:33}

proc {Intersect lines (no divide)} {} {
    GiD_Process MEscape Geometry Create IntersectLines NoDivideLines
}

#[_ "Divide a surface in n parts in u or v direction"]
set {macrosdata(Divide surface,Number)} 39
set {macrosdata(Divide surface,PrePost)} pre
set {macrosdata(Divide surface,CreationDate)} {2006-02-15 23:30:30}
set {macrosdata(Divide surface,Icon)} {surface_divide_nparts.png imported_images surface_divide_nparts.png themed_image}
set {macrosdata(Divide surface,InToolbar)} 1
set {macrosdata(Divide surface,Description)} {Divide a surface in n parts in u or v direction}
set {macrosdata(Divide surface,Accelerators)} {}
set {macrosdata(Divide surface,Group)} {Geometry Edit}
set {macrosdata(Divide surface,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Divide surface,Active)} 1
set {macrosdata(Divide surface,IsStandard)} 1

proc {Divide surface} {} {
    GiD_Process MEscape Geometry Edit DivideSurf NumDivisions
}

#[_ "Split a surface from a path of lines"]
set {macrosdata(Split surface,Description)} {Split a surface from a path of lines}
set {macrosdata(Split surface,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Split surface,PrePost)} pre
set {macrosdata(Split surface,CreationDate)} {2006-02-15 23:19:26}
set {macrosdata(Split surface,IsStandard)} 1
set {macrosdata(Split surface,Active)} 1
set {macrosdata(Split surface,Accelerators)} {}
set {macrosdata(Split surface,Icon)} {surface_divide_split.png imported_images surface_divide_split.png themed_image}
set {macrosdata(Split surface,InToolbar)} 1
set {macrosdata(Split surface,Number)} 40
set {macrosdata(Split surface,Group)} {Geometry Edit}

proc {Split surface} {} {
    GiD_Process MEscape Geometry Edit SplitSurf
}

#[_ "Intersect a surface with lines"]
set {macrosdata(Intersect surface lines,Description)} {Intersect a surface with lines}
set {macrosdata(Intersect surface lines,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Intersect surface lines,PrePost)} pre
set {macrosdata(Intersect surface lines,CreationDate)} {2006-02-16 13:39:46}
set {macrosdata(Intersect surface lines,IsStandard)} 1
set {macrosdata(Intersect surface lines,Active)} 1
set {macrosdata(Intersect surface lines,Accelerators)} {}
set {macrosdata(Intersect surface lines,Icon)} {intersect_line_surface.png imported_images intersect_line_surface.png themed_image}
set {macrosdata(Intersect surface lines,InToolbar)} 1
set {macrosdata(Intersect surface lines,Number)} 41
set {macrosdata(Intersect surface lines,Group)} {Geometry Edit}

proc {Intersect surface lines} {} {
    GiD_Process MEscape Geometry Create IntSurfLine
}

#[_ "Intersect multiple surfaces"]
set {macrosdata(Intersect surfaces,IsStandard)} 1
set {macrosdata(Intersect surfaces,Number)} 42
set {macrosdata(Intersect surfaces,Icon)} {intersect_surfaces.png imported_images intersect_surfaces.png themed_image}
set {macrosdata(Intersect surfaces,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Intersect surfaces,PrePost)} pre
set {macrosdata(Intersect surfaces,Group)} {Geometry Edit}
set {macrosdata(Intersect surfaces,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Intersect surfaces,Description)} {Intersect multiple surfaces}
set {macrosdata(Intersect surfaces,Accelerators)} {}
set {macrosdata(Intersect surfaces,InToolbar)} 1
set {macrosdata(Intersect surfaces,Active)} 1

proc {Intersect surfaces} {} {
    GiD_Process MEscape Geometry Create IntMultSurfs
}

#[_ "Creates a construction horizontal line in a temporal layer, containing one point defined by the user.\n    Length of the line depends on current geometry size"]
set {macrosdata(Construction line horizontal,Active)} 1
set {macrosdata(Construction line horizontal,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Construction line horizontal,Number)} 43
set {macrosdata(Construction line horizontal,Accelerators)} {}
set {macrosdata(Construction line horizontal,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Construction line horizontal,IsStandard)} 1
set {macrosdata(Construction line horizontal,Group)} {Constr. line}
set {macrosdata(Construction line horizontal,Icon)} {border_horizontal.png imported_images border_horizontal.png themed_image}
set {macrosdata(Construction line horizontal,PrePost)} pre
set {macrosdata(Construction line horizontal,Description)} {Creates a construction horizontal line in a temporal layer, containing one point defined by the user.
    Length of the line depends on current geometry size}
set {macrosdata(Construction line horizontal,InToolbar)} 0

proc {Construction line horizontal} {} {
    set layers [GiD_Info Layers]
    set ipos [lsearch -exact $layers *construction_lines*]
    if { $ipos == -1 } { GiD_Process 'Layers New *construction_lines* escape }
    if { $ipos == -1 } { GiD_Process 'Layers Color *construction_lines* 255255000 escape }
    if { $ipos != -1 } { set layers [lreplace $layers $ipos $ipos] }
    set current_layer [GiD_Info Project LayerToUse]
    GiD_Process 'Layers ToUse *construction_lines* escape
    set bbox [eval GiD_Info layers -bbox $layers]
    set max_size [expr {60*[CEGiveSizeFromBBox $bbox]}]
    set p1 [GidUtils::GetCoordinates "Enter point for H line"]
    GiD_Process MEscape Utilities Id FFNoJoin $p1 escape
    GiD_Process MEscape Geometry Create line  @[expr {-.5*$max_size}],0 @$max_size,0 escape escape
    GiD_Process 'Layers ToUse $current_layer escape
}

#[_ "Creates a construction vertical line in a temporal layer, containing one point defined by the user.\n    Length of the line depends on current geometry size"]
set {macrosdata(Construction line vertical,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Construction line vertical,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Construction line vertical,Accelerators)} {}
set {macrosdata(Construction line vertical,Active)} 1
set {macrosdata(Construction line vertical,IsStandard)} 1
set {macrosdata(Construction line vertical,PrePost)} pre
set {macrosdata(Construction line vertical,Group)} {Constr. line}
set {macrosdata(Construction line vertical,Number)} 44
set {macrosdata(Construction line vertical,Description)} {Creates a construction vertical line in a temporal layer, containing one point defined by the user.
    Length of the line depends on current geometry size}
set {macrosdata(Construction line vertical,InToolbar)} 0
set {macrosdata(Construction line vertical,Icon)} {border_vertical.png imported_images border_vertical.png themed_image}

proc {Construction line vertical} {} {
    set layers [GiD_Info Layers]
    set ipos [lsearch -exact $layers *construction_lines*]
    if { $ipos == -1 } { GiD_Process 'Layers New *construction_lines* escape }
    if { $ipos == -1 } { GiD_Process 'Layers Color *construction_lines* 255255000 escape }
    if { $ipos != -1 } { set layers [lreplace $layers $ipos $ipos] }
    set current_layer [GiD_Info Project LayerToUse]
    GiD_Process 'Layers ToUse *construction_lines* escape
    set bbox [eval GiD_Info layers -bbox $layers]
    set max_size [expr {60*[CEGiveSizeFromBBox $bbox]}]
    set p1 [GidUtils::GetCoordinates "Enter point for V line"]
    GiD_Process MEscape Utilities Id FFNoJoin $p1 escape
    GiD_Process MEscape Geometry Create line  @0,[expr {-.5*$max_size}] @0,$max_size escape escape
    GiD_Process 'Layers ToUse $current_layer escape
}

#[_ "Split all triangles with edge size greater than an user size"]
set {macrosdata(Split Big Triangles,Icon)} {split_elements.png imported_images split_elements.png themed_image}
set {macrosdata(Split Big Triangles,Active)} 1
set {macrosdata(Split Big Triangles,CreationDate)} {2006-01-31 12:53:02}
set {macrosdata(Split Big Triangles,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Split Big Triangles,Description)} {Split all triangles with edge size greater than an user size}
set {macrosdata(Split Big Triangles,IsStandard)} 1
set {macrosdata(Split Big Triangles,Accelerators)} {}
set {macrosdata(Split Big Triangles,Number)} 45
set {macrosdata(Split Big Triangles,InToolbar)} 0
set {macrosdata(Split Big Triangles,Group)} {Mesh Edit}
set {macrosdata(Split Big Triangles,PrePost)} pre

proc {Split Big Triangles} {} {
    proc SplitTrianglesEdgeMaxSize { size } {
        if { $size <= 0 } return
        GiD_Process escape escape escape escape Meshing EditMesh SplitElems TriaToTria
        GiD_Process Filter:MAXLENGTH=[string trim $size]
        GiD_Process 1:
        GiD_Process MEscape Meshing EditMesh Smoothing
        GiD_Process 1:
        GiD_Process MEscape
    }
    proc SplitTrianglesEdgeMaxSizeW { } {
        if { ![info exists ::MaxEdgeSize] } {
            set ::MaxEdgeSize 1.0
        }
        set size [tk_dialogEntryRAM .gid.maxedgesize [_ "Split triangles"] [_ "Enter edge size"] gidquestionhead real+ $::MaxEdgeSize ""]
        if { $size == "--CANCEL--" } {
            return 1
        } else {
            set ::MaxEdgeSize $size
            SplitTrianglesEdgeMaxSize $size
        }
    }
    SplitTrianglesEdgeMaxSizeW
}

#[_ "Collapse all triangles with edge size smaller than an user size"]
set {macrosdata(Collapse Small Triangles,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Collapse Small Triangles,CreationDate)} {2006-01-31 15:41:36}
set {macrosdata(Collapse Small Triangles,InToolbar)} 0
set {macrosdata(Collapse Small Triangles,Group)} {Mesh Edit}
set {macrosdata(Collapse Small Triangles,IsStandard)} 1
set {macrosdata(Collapse Small Triangles,Active)} 1
set {macrosdata(Collapse Small Triangles,Description)} {Collapse all triangles with edge size smaller than an user size}
set {macrosdata(Collapse Small Triangles,Accelerators)} {}
set {macrosdata(Collapse Small Triangles,PrePost)} pre
set {macrosdata(Collapse Small Triangles,Icon)} {collapse_edges.png imported_images collapse_edges.png themed_image}
set {macrosdata(Collapse Small Triangles,Number)} 46

proc {Collapse Small Triangles} {} {
    proc CollapseTrianglesEdgeMinSize { size } {
        if { $size <= 0 } return
        GiD_Process escape escape escape escape Utilities Collapse Edges
        GiD_Process Tolerance [string trim $size]
        GiD_Process 1:
        GiD_Process MEscape Meshing EditMesh Smoothing
        GiD_Process 1:
        GiD_Process MEscape
    }
    proc CollapseTrianglesEdgeMinSizeW { } {
        if { ![info exists ::MinEdgeSize] } {
            set ::MinEdgeSize 1.0
        }
        set size [tk_dialogEntryRAM .gid.edgesize [_ "Edge size"] [_ "Enter edge size"] gidquestionhead real+ $::MinEdgeSize ""]
        if { $size == "--CANCEL--" } {
            return 1
        } else {
            set ::MinEdgeSize $size
            CollapseTrianglesEdgeMinSize $size
        }
    }
    CollapseTrianglesEdgeMinSizeW
}

#[_ "Align the quadratic nodes to the mid edges"]
set {macrosdata(Align Quadratic Nodes,Active)} 1
set {macrosdata(Align Quadratic Nodes,InToolbar)} 0
set {macrosdata(Align Quadratic Nodes,Group)} {Mesh Edit}
set {macrosdata(Align Quadratic Nodes,IsStandard)} 1
set {macrosdata(Align Quadratic Nodes,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Align Quadratic Nodes,PrePost)} pre
set {macrosdata(Align Quadratic Nodes,Number)} 47
set {macrosdata(Align Quadratic Nodes,Icon)} {align_quadratic_nodes.png imported_images align_quadratic_nodes.png themed_image}
set {macrosdata(Align Quadratic Nodes,CreationDate)} {2007-05-03 21:01:46}
set {macrosdata(Align Quadratic Nodes,Description)} {Align the quadratic nodes to the mid edges}
set {macrosdata(Align Quadratic Nodes,Accelerators)} {}

proc {Align Quadratic Nodes} {} {
    proc CenterQuadraticNodes { } {
        set isquadratic [GiD_Set Model(QuadraticType)]
        if { !$isquadratic } {
            return 1
        }
        ::GidUtils::WaitState .
        set t0 [clock seconds]
        foreach elemtype [lrange [GiD_Info mesh] 1 end] {
            set elems [GiD_Info mesh elements $elemtype]
            if { $elemtype == "Linear" } {
                foreach {num na nb ni mat} $elems {
                    set pa [lindex [GiD_Info Coordinates $na mesh] 0]
                    set pb [lindex [GiD_Info Coordinates $nb mesh] 0]
                    set pi [MathUtils::ScalarByVectorProd 0.5 [MathUtils::VectorSum $pa $pb]]
                    GiD_Mesh edit node $ni $pi
                }
            } elseif { $elemtype == "Triangle" } {
                foreach {num n0 n1 n2 n3 n4 n5 mat} $elems {
                    foreach ni [list $n3 $n4 $n5] na [list $n0 $n1 $n2] nb [list $n1 $n2 $n0] {
                        set pa [lindex [GiD_Info Coordinates $na mesh] 0]
                        set pb [lindex [GiD_Info Coordinates $nb mesh] 0]
                        set pi [MathUtils::ScalarByVectorProd 0.5 [MathUtils::VectorSum $pa $pb]]
                        GiD_Mesh edit node $ni $pi
                    }
                }
            } elseif { $elemtype == "Quadrilateral" } {
                if { $isquadratic == 1 } {
                    foreach {num n0 n1 n2 n3 n4 n5 n6 n7 mat} $elems {
                        foreach ni [list $n4 $n5 $n6 $n7] na [list $n0 $n1 $n2 $n3] nb [list $n1 $n2 $n3 $n0] {
                            set pa [lindex [GiD_Info Coordinates $na mesh] 0]
                            set pb [lindex [GiD_Info Coordinates $nb mesh] 0]
                            set pi [MathUtils::ScalarByVectorProd 0.5 [MathUtils::VectorSum $pa $pb]]
                            GiD_Mesh edit node $ni $pi
                        }
                    }
                } else {
                    foreach {num n0 n1 n2 n3 n4 n5 n6 n7 n8 mat} $elems {
                        foreach ni [list $n4 $n5 $n6 $n7] na [list $n0 $n1 $n2 $n3] nb [list $n1 $n2 $n3 $n0] {
                            set pa [lindex [GiD_Info Coordinates $na mesh] 0]
                            set pb [lindex [GiD_Info Coordinates $nb mesh] 0]
                            set pi [MathUtils::ScalarByVectorProd 0.5 [MathUtils::VectorSum $pa $pb]]
                            GiD_Mesh edit node $ni $pi
                        }
                        set pi [lindex [GiD_Info Coordinates $n0 mesh] 0]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n1 mesh] 0]]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n2 mesh] 0]]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n3 mesh] 0]]
                        set pi [MathUtils::ScalarByVectorProd 0.25 $pi]
                        GiD_Mesh edit node $n8 $pi
                    }
                }
            } elseif { $elemtype == "Tetrahedra" } {
                foreach {num n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 mat} $elems {
                    foreach ni [list $n5 $n6 $n7 $n8 $n9 $n10]  na [list $n1 $n2 $n1 $n1 $n2 $n3]  nb [list $n2 $n3 $n3 $n4 $n4 $n4] {
                        set pa [lindex [GiD_Info Coordinates $na mesh] 0]
                        set pb [lindex [GiD_Info Coordinates $nb mesh] 0]
                        set pi [MathUtils::ScalarByVectorProd 0.5 [MathUtils::VectorSum $pa $pb]]
                        GiD_Mesh edit node $ni $pi
                    }
                }
            } elseif { $elemtype == "Hexahedra" } {
                if { $isquadratic == 1 } {
                    foreach {num n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14 n15 n16 n17 n18 n19 n20 mat} $elems {
                        foreach ni [list $n9 $n10 $n11 $n12 $n13 $n14 $n15 $n16 $n17 $n18 $n19 $n20]  na [list $n1 $n2  $n3  $n1  $n1  $n2  $n3  $n4  $n5  $n6  $n7  $n5 ]  nb [list $n2 $n3  $n4  $n4  $n5  $n6  $n7  $n8  $n6  $n7  $n8  $n8 ] {
                            set pa [lindex [GiD_Info Coordinates $na mesh] 0]
                            set pb [lindex [GiD_Info Coordinates $nb mesh] 0]
                            set pi [MathUtils::ScalarByVectorProd 0.5 [MathUtils::VectorSum $pa $pb]]
                            GiD_Mesh edit node $ni $pi
                        }
                    }
                } else {
                    foreach {num n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14 n15 n16 n17 n18 n19 n20  n21 n22 n23 n24 n25 n26 n27 mat} $elems {
                        foreach ni [list $n9 $n10 $n11 $n12 $n13 $n14 $n15 $n16 $n17 $n18 $n19 $n20]  na [list $n1 $n2  $n3  $n1  $n1  $n2  $n3  $n4  $n5  $n6  $n7  $n5 ]  nb [list $n2 $n3  $n4  $n4  $n5  $n6  $n7  $n8  $n6  $n7  $n8  $n8 ] {
                            set pa [lindex [GiD_Info Coordinates $na mesh] 0]
                            set pb [lindex [GiD_Info Coordinates $nb mesh] 0]
                            set pi [MathUtils::ScalarByVectorProd 0.5 [MathUtils::VectorSum $pa $pb]]
                            GiD_Mesh edit node $ni $pi
                        }
                        foreach ni [list $n21 $n22 $n23 $n24 $n25 $n26]  na [list $n1  $n1  $n2  $n3  $n1  $n5 ]  nb [list $n2  $n2  $n3  $n7  $n4  $n6 ]  nc [list $n3  $n6  $n7  $n8  $n8  $n7 ]  nd [list $n4  $n5  $n6  $n4  $n5  $n8 ] {
                            set pi [lindex [GiD_Info Coordinates $na mesh] 0]
                            set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $nb mesh] 0]]
                            set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $nc mesh] 0]]
                            set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $nd mesh] 0]]
                            set pi [MathUtils::ScalarByVectorProd 0.25 $pi]
                            GiD_Mesh edit node $ni $pi
                        }
                        set pi [lindex [GiD_Info Coordinates $n1 mesh] 0]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n2 mesh] 0]]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n3 mesh] 0]]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n4 mesh] 0]]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n5 mesh] 0]]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n6 mesh] 0]]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n7 mesh] 0]]
                        set pi [MathUtils::VectorSum $pi [lindex [GiD_Info Coordinates $n8 mesh] 0]]
                        set pi [MathUtils::ScalarByVectorProd 0.125 $pi]
                    }
                }
            } elseif { $elemtype == "Tetrahedra" } {
                foreach {num n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14 n15 mat} $elems {
                    foreach ni [list $n7 $n8 $n9 $n10 $n11 $n12 $n13 $n14 $n15]  na [list $n1 $n2 $n1 $n1  $n2  $n3  $n4  $n5  $n4 ]  nb [list $n2 $n3 $n3 $n4  $n5  $n6  $n5  $n6  $n6 ] {
                        set pa [lindex [GiD_Info Coordinates $na mesh] 0]
                        set pb [lindex [GiD_Info Coordinates $nb mesh] 0]
                        set pi [MathUtils::ScalarByVectorProd 0.5 [MathUtils::VectorSum $pa $pb]]
                        GiD_Mesh edit node $ni $pi
                    }
                }
            } else {
                #other elements
            }
        }
        ::GidUtils::EndWaitState .
        ::GidUtils::SetWarnLine [_ "Quadratic elements set to linear shape. Time= %s seconds" [expr [clock seconds]-$t0]]
        GiD_Redraw
    }
    CenterQuadraticNodes
}

#[_ "Read preferences from a file (to avoid restart the system with other user)"]
set {macrosdata(Read Preferences,IsStandard)} 1
set {macrosdata(Read Preferences,Description)} {Read preferences from a file (to avoid restart the system with other user)}
set {macrosdata(Read Preferences,Active)} 1
set {macrosdata(Read Preferences,PrePost)} pre
set {macrosdata(Read Preferences,Group)} Tools
set {macrosdata(Read Preferences,CreationDate)} {2006-02-13 20:59:36}
set {macrosdata(Read Preferences,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Read Preferences,Icon)} {user_open.png imported_images user_open.png themed_image}
set {macrosdata(Read Preferences,Number)} 48
set {macrosdata(Read Preferences,Accelerators)} {}
set {macrosdata(Read Preferences,InToolbar)} 0

proc {Read Preferences} {} {
    proc ReadDefaultsFromFile { filename } {
        ReadDefaultValues $filename ;#tcl values
        readdefaults $filename 0 ;#kernel values
        if  { [winfo exists .gid.gidPreferences] } {
            FillPrefValuesFromGid
        }
    }
    proc ReadDefaultsFromFileW { } {
        set w .gid
        set types [list [list [_ "Preferences file"] ".ini"] [list [_ "All files"] ".*"]]
        set filename [Browser-ramR file read $w [_ "Read preferences file"] "" $types .ini 0]
        if { $filename == "" } { return }
        ReadDefaultsFromFile $filename
    }
    ReadDefaultsFromFileW
}

#[_ "Save preferences to a file (to avoid restart the system with other user)"]
set {macrosdata(Save Preferences,InToolbar)} 0
set {macrosdata(Save Preferences,Description)} {Save preferences to a file (to avoid restart the system with other user)}
set {macrosdata(Save Preferences,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Save Preferences,Number)} 49
set {macrosdata(Save Preferences,Group)} Tools
set {macrosdata(Save Preferences,IsStandard)} 1
set {macrosdata(Save Preferences,CreationDate)} {2006-02-13 21:49:11}
set {macrosdata(Save Preferences,PrePost)} pre
set {macrosdata(Save Preferences,Icon)} {user_save.png imported_images user_save.png themed_image}
set {macrosdata(Save Preferences,Active)} 1
set {macrosdata(Save Preferences,Accelerators)} {}

proc {Save Preferences} {} {
    proc SaveDefaultsToFile { filename } {
        writedefaults $filename;#kernel values
        WriteTCLDefaultsInFile $filename ;#tcl values
    }
    proc SaveDefaultsToFileW { } {
        set w .gid
        set types [list [list [_ "Preferences file"] ".ini"] [list [_ "All files"] ".*"]]
        set filename [Browser-ramR file write $w [_ "Save preferences file"] "" $types .ini 0]
        if { $filename == "" } { return }
        SaveDefaultsToFile $filename
    }
    SaveDefaultsToFileW
}

#[_ "Delete duplicated near points with a tolerance"]
set {macrosdata(Collapse points,CreationDate)} {2006-02-14 23:37:21}
set {macrosdata(Collapse points,Active)} 1
set {macrosdata(Collapse points,Icon)} {collapse_points.png imported_images collapse_points.png themed_image}
set {macrosdata(Collapse points,PrePost)} pre
set {macrosdata(Collapse points,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Collapse points,Accelerators)} {}
set {macrosdata(Collapse points,Number)} 50
set {macrosdata(Collapse points,Description)} {Delete duplicated near points with a tolerance}
set {macrosdata(Collapse points,InToolbar)} 0
set {macrosdata(Collapse points,Group)} {Geometry Edit}
set {macrosdata(Collapse points,IsStandard)} 1

proc {Collapse points} {} {
    GiD_Process MEscape Utilities Collapse Points
}

#[_ "Delete duplicated near lines with a tolerance"]
set {macrosdata(Collapse lines,Active)} 1
set {macrosdata(Collapse lines,Icon)} {collapse_lines.png imported_images collapse_lines.png themed_image}
set {macrosdata(Collapse lines,Number)} 51
set {macrosdata(Collapse lines,InToolbar)} 0
set {macrosdata(Collapse lines,CreationDate)} {2006-02-14 23:42:43}
set {macrosdata(Collapse lines,Group)} {Geometry Edit}
set {macrosdata(Collapse lines,Accelerators)} {}
set {macrosdata(Collapse lines,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Collapse lines,Description)} {Delete duplicated near lines with a tolerance}
set {macrosdata(Collapse lines,PrePost)} pre
set {macrosdata(Collapse lines,IsStandard)} 1

proc {Collapse lines} {} {
    GiD_Process MEscape Utilities Collapse Lines
}

#[_ "Delete duplicated entities of the model with a tolerance"]
set {macrosdata(Collapse  model,PrePost)} pre
set {macrosdata(Collapse  model,Active)} 1
set {macrosdata(Collapse  model,Description)} {Delete duplicated entities of the model with a tolerance}
set {macrosdata(Collapse  model,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Collapse  model,Number)} 52
set {macrosdata(Collapse  model,CreationDate)} {2006-02-14 23:39:16}
set {macrosdata(Collapse  model,InToolbar)} 1
set {macrosdata(Collapse  model,Group)} {Geometry Edit}
set {macrosdata(Collapse  model,IsStandard)} 1
set {macrosdata(Collapse  model,Icon)} {collapse_alltypes.png imported_images collapse_alltypes.png themed_image}
set {macrosdata(Collapse  model,Accelerators)} {}

proc {Collapse  model} {} {
    GiD_Process MEscape Utilities Collapse Model
}

#[_ "List some filtered results"]
set macrosdata(ListResults,ModificationDate) {2011-11-04 11:03:51}
set macrosdata(ListResults,PrePost) post
set macrosdata(ListResults,CreationDate) {2006-02-15 11:38:42}
set macrosdata(ListResults,Group) View
set macrosdata(ListResults,Active) 1
set macrosdata(ListResults,InToolbar) 1
set macrosdata(ListResults,Accelerators) {}
set macrosdata(ListResults,Description) {List some filtered results}
set macrosdata(ListResults,Number) 53
set macrosdata(ListResults,IsStandard) 1
set macrosdata(ListResults,Icon) {align_text_left.png imported_images align_text_left.png themed_image}

proc ListResults {} {
    namespace eval ::PostUtils {
    }
    #over: nodes or elements
    #results: list with this possible options: nodesnum, nodescoord, elemsnum,
    #         elemsconec, elemstype, elemslayer, <result name>
    #selection: list of id's (can be used a:b)
    #presentation: TEXT or LIST
    proc ::PostUtils::SelectResults { results over selection {presentation TEXT} } {
        set info [eval GiD_Info list_entities $over $selection]
        set info [split $info \n]
        
        set res ""
        if { $over == "nodes" } {
            if { [GiD_Info postprocess get cur_result] != "" } {
                set infoitems {numinfo usedbyinfo currentinfo moreinfo}
            } else {
                set infoitems {numinfo usedbyinfo moreinfo}
            }
            foreach $infoitems [lrange $info 1 end-1] {
                set res2 ""
                if { [lsearch $results "nodesnum"] != -1 } {
                    lappend res2 [lindex $numinfo 1]
                }
                if { [lsearch $results "nodescoord"] != -1 } {
                    lappend res2 [string trim [lindex $numinfo 2]]
                }
                foreach resultinfo [lindex $moreinfo 2] {
                    set name [lindex $resultinfo 0]
                    if { [lsearch $results $name] != -1 } {
                        if { [string trim [lindex $resultinfo 1]] != "" } {
                            set resultnames($name) ""
                            foreach item [string trim [lindex $resultinfo 1]] {
                                lappend resultnames($name) $item
                            }
                        }
                        lappend res2 [string trim [lindex [lindex $resultinfo 2] 0]]
                    }
                }
                lappend res $res2
            }
        } elseif { $over == "elements" } {
            if { [GiD_Info postprocess get cur_result] != "" } {
                set infoitems {numinfo currentinfo gausspoints moreinfo}
            } else {
                set infoitems {numinfo gausspoints moreinfo}
            }
            foreach $infoitems [lrange $info 1 end-1] {
                set res2 ""
                if { [lsearch $results "elemsnum"] != -1 } {
                    lappend res2 [lindex $numinfo 1]
                }
                if { [lsearch $results "elemsconec"] != -1 } {
                    lappend res2 [string trim [lindex $numinfo 2]]
                }
                if { [lsearch $results "elemstype"] != -1 } {
                    lappend res2 [string trim [lindex $numinfo 3]]
                }
                if { [lsearch $results "elemslayer"] != -1 } {
                    lappend res2 [string trim [lindex $numinfo 4]]
                }
                foreach resultinfo [lindex $moreinfo 2] {
                    set name [lindex $resultinfo 0]
                    if { [lsearch $results $name] != -1 } {
                        if { [string trim [lindex $resultinfo 2]] != "" } {
                            set resultnames($name) ""
                            foreach item [string trim [lindex $resultinfo 2]] {
                                lappend resultnames($name) $item
                            }
                        }
                        lappend res2 [string trim [lindex $resultinfo 3]]
                    }
                }
                lappend res $res2
            }
        }
        if { $presentation == "TEXT" } {
            set newres ""
            foreach item $results {
                if { [info exists resultnames($item)] && $resultnames($item) != "" } {
                    foreach component $resultnames($item) {
                        set txt $item:$component
                        regsub -all { } $txt _ txt
                        append newres $txt " "
                    }
                } else {
                    if { $item == "nodesnum" || $item == "elemsnum" } {
                        set txt [_ "Number"]
                        regsub -all { } $txt _ txt
                    } elseif { $item == "nodescoord" } {
                        set txt "X Y Z"
                    } elseif { $item == "elemsconec" } {
                        set txt [_ "connectivities"]
                    } else {
                        set txt $item
                        regsub -all { } $txt _ txt
                    }
                    append newres $txt " "
                }
            }
            set newres [string trim $newres]\n
            foreach item $res {
                #append newres [eval concat [eval concat $item]]\n
                append newres [eval concat $item]\n
            }
            set res $newres
        }
        WarnWinText $res
        return $res
    }
    
    proc ::PostUtils::ApplySelectResults { w } {
        global PostUtils
        set results ""
        if { $PostUtils(over) == "nodes" } {
            if { $PostUtils(output,number) } {
                lappend results nodesnum
            }
            if { $PostUtils(output,coordinates) } {
                lappend results nodescoord
            }
        } else {
            if { $PostUtils(output,number) } {
                lappend results elemsnum
            }
            if { $PostUtils(output,connectivities) } {
                lappend results elemsconec
            }
        }
        foreach item [GiD_Info postprocess get cur_results_list contour_fill] {
            regsub -all _ $item { } item
            if { $PostUtils(over,$item) == $PostUtils(over) && $PostUtils(output,$item) } {
                lappend results $item
            }
        }
        if { [llength $results] > 0 } {
            set SmallWinSelecting [GiD_Set SmallWinSelecting]
            FinishButton $w $w.buts [_ "Press 'Finish' to end selection"] "" disableall $SmallWinSelecting
            set entities [GidUtils::PickEntities $PostUtils(over) multiple]
            EndFinishButton $w "" $SmallWinSelecting
            if { $entities != "" } {
                ::PostUtils::SelectResults $results $PostUtils(over) $entities]
            } else {
                GidUtils::SetWarnLine [_ "Must select some entity"]
            }
        } else {
            GidUtils::SetWarnLine [_ "Must select some result"]
        }
    }
    
    proc ::PostUtils::SelectResultsW { { w .gid.selectresults } } {
        global PostUtils
        InitWindow $w [_ "Select Results"] SelectResultsWindowGeom ::PostUtils::SelectResultsW
        ttk::frame $w.f -style ridge.TFrame -borderwidth 2
        ttk::frame $w.f.f0
        ttk::radiobutton $w.f.f0.rb0 -text [_ "Nodes"] -value nodes -variable PostUtils(over)
        ttk::radiobutton $w.f.f0.rb1 -text [_ "Elements"] -value elements -variable PostUtils(over)
        if { ![info exists PostUtils(over)] || ($PostUtils(over) != "nodes" && $PostUtils(over) != "elements") } {
            set PostUtils(over) nodes
        }
        if { ![info exists PostUtils(output,number)] } {
            set PostUtils(output,number) 1
        }
        if { ![info exists PostUtils(output,coordinates)] } { set PostUtils(output,coordinates) 0 }
        ttk::checkbutton $w.f.f0.cbcoordinates -text [_ "Coordinates"] -variable PostUtils(output,coordinates)
        if { $PostUtils(over) == "nodes" } {
            $w.f.f0.cbcoordinates configure -state normal
        } else {
            $w.f.f0.cbcoordinates configure -state disable
        }
        if { ![info exists PostUtils(output,connectivities)] } { set PostUtils(output,connectivities) 0 }
        ttk::checkbutton $w.f.f0.cbconnectivities -text [_ "Connectivities"] -variable PostUtils(output,connectivities)
        if { $PostUtils(over) == "nodes" } {
            $w.f.f0.cbconnectivities configure -state disable
        } else {
            $w.f.f0.cbconnectivities configure -state normal
        }
        if { ![info exists PostUtils(output,number)] } { set PostUtils(output,number) 0 }
        ttk::checkbutton $w.f.number -text [_ "Number"] -variable PostUtils(output,number)
        #array unset PostUtils output,*
        set results [GiD_Info postprocess get cur_results_list contour_fill]
        set cont 0
        foreach item $results {
            set header [GiD_Result get -info [list $item [GiD_Info postprocess get cur_analysis] [GiD_Info postprocess get cur_step]]]
            regsub -all _ $item { } item
            if { [lindex [lindex $header 0] 5] == "OnNodes" } {
                set PostUtils(over,$item) nodes
            } else {
                set PostUtils(over,$item) elements
            }
            if { ![info exists PostUtils(output,$item)] } { set PostUtils(output,$item) 0 }
            ttk::checkbutton $w.f.cb$cont -text [TranslateResultName $item] -variable PostUtils(output,$item)
            if { $PostUtils(over,$item) == $PostUtils(over) } {
                $w.f.cb$cont configure -state normal
            } else {
                $w.f.cb$cont configure -state disable
            }
            incr cont
        }
        grid $w.f.f0.rb0 $w.f.f0.rb1 -sticky w
        grid $w.f.f0.cbcoordinates $w.f.f0.cbconnectivities -sticky w
        grid $w.f.f0 -sticky ew
        grid $w.f.number -sticky w
        for {set i 0} {$i < $cont} {incr i} {
            grid $w.f.cb$i -sticky w
        }
        grid rowconfigure $w.f [expr {$cont+2}] -weight 1
        grid columnconfigure $w.f 0 -weight 1
        grid $w.f -sticky nsew    
        ttk::frame $w.buts -style BottomFrame.TFrame  
        ttk::button $w.buts.apply -text [_ "Select"] -command "::PostUtils::ApplySelectResults $w" -style BottomFrame.TButton
        ttk::button $w.buts.close -text [_ "Close"] -command "destroy $w" -style BottomFrame.TButton
        
        grid $w.buts.apply $w.buts.close -sticky ew -padx 5 -pady 6
        grid $w.buts -sticky sew
        if { $::tcl_version >= 8.5 } { grid anchor $w.buts center }
        grid columnconfigure $w 0 -weight 1
        grid rowconfigure $w 0 -weight 1   
        bind $w <Destroy> [list +::PostUtils::DestroySelectResults %W $w] ;# + to add to previous script
        trace add variable PostUtils(over) write "::PostUtils::ChangeSelectResultsOver $w ;#"
    }
    
    proc ::PostUtils::DestroySelectResults { W w } {
        if { $W != $w } return
        #reenter multiple times, one by toplevel child, only interest w
        global PostUtils
        #catch { array unset PostUtils result,* } ;#remove also traces
        trace remove variable PostUtils(over) write "::PostUtils::ChangeSelectResultsOver $w ;#"
    }
    
    proc ::PostUtils::ChangeSelectResultsOver { w } {
        global PostUtils
        if { $PostUtils(over) == "nodes" } {
            $w.f.f0.cbcoordinates configure -state normal
            $w.f.f0.cbconnectivities configure -state disable
        } else {
            $w.f.f0.cbcoordinates configure -state disable
            $w.f.f0.cbconnectivities configure -state normal
        }
        set results [GiD_Info postprocess get cur_results_list contour_fill]
        set cont 0
        foreach item $results {
            regsub -all _ $item { } item
            if { $PostUtils(over,$item) == $PostUtils(over) } {
                $w.f.cb$cont configure -state normal
            } else {
                $w.f.cb$cont configure -state disable
            }
            incr cont
        }
    }
    ::PostUtils::SelectResultsW
}

#[_ "set all labels on/off"]
set {macrosdata(Label on/off,Description)} {set all labels on/off}
set {macrosdata(Label on/off,Active)} 1
set {macrosdata(Label on/off,IsStandard)} 1
set {macrosdata(Label on/off,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Label on/off,CreationDate)} {2006-02-15 23:50:41}
set {macrosdata(Label on/off,Group)} View
set {macrosdata(Label on/off,Number)} 54
set {macrosdata(Label on/off,Accelerators)} {}
set {macrosdata(Label on/off,Icon)} {view_labels.png imported_images view_labels.png themed_image}
set {macrosdata(Label on/off,PrePost)} prepost
set {macrosdata(Label on/off,InToolbar)} 1

proc {Label on/off} {} {
    proc SetLabelsOnOff { } {
        global GidPriv
        if { ![info exists GidPriv(labelson)] } {
            set GidPriv(labelson) 0
        }
        if { $GidPriv(labelson) } {
            set  GidPriv(labelson) 0
            GiD_Process 'Label Off
        } else {
            set  GidPriv(labelson) 1
            GiD_Process 'Label All
        }
    }
    SetLabelsOnOff
}

#[_ "Draw normals of surfaces by color"]
set {macrosdata(Draw normals surfaces color,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Draw normals surfaces color,Accelerators)} {}
set {macrosdata(Draw normals surfaces color,Active)} 1
set {macrosdata(Draw normals surfaces color,IsStandard)} 1
set {macrosdata(Draw normals surfaces color,PrePost)} pre
set {macrosdata(Draw normals surfaces color,Group)} View
set {macrosdata(Draw normals surfaces color,Number)} 55
set {macrosdata(Draw normals surfaces color,Description)} {Draw normals of surfaces by color}
set {macrosdata(Draw normals surfaces color,InToolbar)} 1
set {macrosdata(Draw normals surfaces color,Icon)} {view_normals_surface.png imported_images view_normals_surface.png themed_image}
set {macrosdata(Draw normals surfaces color,CreationDate)} {2006-02-16 13:04:57}

proc {Draw normals surfaces color} {} {
    proc DrawNormalsColor { } {
        if { [GiD_Info project ViewMode] == "GEOMETRYUSE" } {
            GiD_Process Mescape Utilities DrawNormals Surfaces Color
        } else {
            GiD_Process Mescape Utilities DrawNormals Color
        }
        foreach layer [GiD_Info Layers] {
            foreach {on frozen} [lrange [lindex [GiD_Info Layers $layer] 0]  0 1] break
            if { $on && !$frozen } {
                GiD_Process 'Layer:$layer
            }
        }
    }
    DrawNormalsColor
}

#[_ "Show number of higerentities (parents) of lines"]
set {macrosdata(Higherentities lines,PrePost)} pre
set {macrosdata(Higherentities lines,CreationDate)} {2006-02-16 13:23:38}
set {macrosdata(Higherentities lines,Active)} 1
set {macrosdata(Higherentities lines,Group)} View
set {macrosdata(Higherentities lines,InToolbar)} 1
set {macrosdata(Higherentities lines,Accelerators)} {}
set {macrosdata(Higherentities lines,Icon)} {view_higherentity_lines.png imported_images view_higherentity_lines.png themed_image}
set {macrosdata(Higherentities lines,Number)} 56
set {macrosdata(Higherentities lines,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Higherentities lines,IsStandard)} 1
set {macrosdata(Higherentities lines,Description)} {Show number of higerentities (parents) of lines}

proc {Higherentities lines} {} {
    if { [GiD_Info Project ViewMode] == "GEOMETRYUSE" } {
        GiD_Process Mescape Utilities DrawHigher Lines
    } else {
        GiD_Process Mescape Utilities DrawHigher Edges
    }
}

#[_ "Draw assigned mesh sizes"]
set {macrosdata(Draw mesh sizes,InToolbar)} 1
set {macrosdata(Draw mesh sizes,Icon)} {view_assigned_meshsizes.png imported_images view_assigned_meshsizes.png themed_image}
set {macrosdata(Draw mesh sizes,PrePost)} pre
set {macrosdata(Draw mesh sizes,Accelerators)} {}
set {macrosdata(Draw mesh sizes,IsStandard)} 1
set {macrosdata(Draw mesh sizes,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Draw mesh sizes,Group)} View
set {macrosdata(Draw mesh sizes,Description)} {Draw assigned mesh sizes}
set {macrosdata(Draw mesh sizes,Active)} 1
set {macrosdata(Draw mesh sizes,CreationDate)} {2006-02-16 16:52:44}
set {macrosdata(Draw mesh sizes,Number)} 57

proc {Draw mesh sizes} {} {
    GiD_Process Mescape Meshing DrawSizes All
}

#[_ "View XY plane"]
set {macrosdata(XY View,PrePost)} prepost
set {macrosdata(XY View,Active)} 1
set {macrosdata(XY View,Group)} View
set {macrosdata(XY View,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(XY View,Icon)} {view_rotate_xy.png imported_images view_rotate_xy.png themed_image}
set {macrosdata(XY View,Number)} 58
set {macrosdata(XY View,IsStandard)} 1
set {macrosdata(XY View,CreationDate)} {2006-02-16 20:33:42}
set {macrosdata(XY View,InToolbar)} 1
set {macrosdata(XY View,Accelerators)} {}
set {macrosdata(XY View,Description)} {View XY plane}

proc {XY View} {} {
    GiD_Process 'Rotate Angle 270 90
}

#[_ "View XZ plane"]
set {macrosdata(XZ View,PrePost)} prepost
set {macrosdata(XZ View,IsStandard)} 1
set {macrosdata(XZ View,Icon)} {view_rotate_xz.png imported_images view_rotate_xz.png themed_image}
set {macrosdata(XZ View,CreationDate)} {2006-02-16 20:36:44}
set {macrosdata(XZ View,Group)} View
set {macrosdata(XZ View,Active)} 1
set {macrosdata(XZ View,Accelerators)} {}
set {macrosdata(XZ View,Description)} {View XZ plane}
set {macrosdata(XZ View,Number)} 59
set {macrosdata(XZ View,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(XZ View,InToolbar)} 0

proc {XZ View} {} {
    GiD_Process 'Rotate Angle 270 0
}

#[_ "View YZ plane"]
set {macrosdata(YZ View,Icon)} {view_rotate_yz.png imported_images view_rotate_yz.png themed_image}
set {macrosdata(YZ View,Description)} {View YZ plane}
set {macrosdata(YZ View,Active)} 1
set {macrosdata(YZ View,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(YZ View,Group)} View
set {macrosdata(YZ View,InToolbar)} 0
set {macrosdata(YZ View,CreationDate)} {2006-02-16 20:37:05}
set {macrosdata(YZ View,Number)} 60
set {macrosdata(YZ View,IsStandard)} 1
set {macrosdata(YZ View,PrePost)} prepost
set {macrosdata(YZ View,Accelerators)} {}

proc {YZ View} {} {
    GiD_Process 'Rotate Angle 0 0
}

#[_ "View - XY plane"]
set {macrosdata(-XY View,CreationDate)} {2006-02-16 20:37:56}
set {macrosdata(-XY View,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(-XY View,Active)} 1
set {macrosdata(-XY View,Description)} {View - XY plane}
set {macrosdata(-XY View,IsStandard)} 1
set {macrosdata(-XY View,PrePost)} prepost
set {macrosdata(-XY View,Icon)} {view_rotate_-xy.png imported_images view_rotate_-xy.png themed_image}
set {macrosdata(-XY View,Group)} View
set {macrosdata(-XY View,Accelerators)} {}
set {macrosdata(-XY View,InToolbar)} 0
set {macrosdata(-XY View,Number)} 61

proc {-XY View} {} {
    GiD_Process 'Rotate Angle 90 -90
}

#[_ "View - YZ plane"]
set {macrosdata(-YZ View,CreationDate)} {2006-02-16 20:38:46}
set {macrosdata(-YZ View,InToolbar)} 0
set {macrosdata(-YZ View,Icon)} {view_rotate_-yz.png imported_images view_rotate_-yz.png themed_image}
set {macrosdata(-YZ View,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(-YZ View,Active)} 1
set {macrosdata(-YZ View,Accelerators)} {}
set {macrosdata(-YZ View,IsStandard)} 1
set {macrosdata(-YZ View,Description)} {View - YZ plane}
set {macrosdata(-YZ View,Group)} View
set {macrosdata(-YZ View,Number)} 62
set {macrosdata(-YZ View,PrePost)} prepost

proc {-YZ View} {} {
    GiD_Process 'Rotate Angle 180 0
}

#[_ "Change the render mode"]
set macrosdata(SwapRender,CreationDate) {2008-11-03 21:44:46}
set macrosdata(SwapRender,Active) 1
set macrosdata(SwapRender,PrePost) prepost
set macrosdata(SwapRender,Accelerators) F3
set macrosdata(SwapRender,Number) 63
set macrosdata(SwapRender,Group) View
set macrosdata(SwapRender,Description) {Change the render mode}
set macrosdata(SwapRender,Icon) {view_swap_render.png imported_images view_swap_render.png themed_image}
set macrosdata(SwapRender,InToolbar) 0
set macrosdata(SwapRender,ModificationDate) {2011-11-04 11:03:51}
set macrosdata(SwapRender,IsStandard) 1

proc SwapRender {} {
    set render_mode [ GiD_Info Project RenderMode]
    if { $render_mode != "postprocess"} {
        if { $render_mode !=  "normal" } {
            GiD_Process 'Render Normal
        } else {
            GiD_Process 'Render FlatLighting
        }
    } else {
        set postrender [ GiD_Info postprocess get cur_display_render]
        
        set lstRender [ list Normal Flat Smooth]
        set idxRender [ lsearch $lstRender $postrender]
        set nRender [ llength $lstRender]
        set nextRender [ expr ( $idxRender + $nRender + 1) % $nRender]
        GiD_Process 'Render [ lindex $lstRender $nextRender]    
    }
}

#[_ "converts scalar results to a vector"]
set {macrosdata(Scalar 2 vector,Icon)} {Scalar2Vector.png imported_images Scalar2Vector.png themed_image}
set {macrosdata(Scalar 2 vector,Number)} 64
set {macrosdata(Scalar 2 vector,IsStandard)} 1
set {macrosdata(Scalar 2 vector,Description)} {converts scalar results to a vector}
set {macrosdata(Scalar 2 vector,InToolbar)} 1
set {macrosdata(Scalar 2 vector,Group)} Tools
set {macrosdata(Scalar 2 vector,CreationDate)} {2007-02-28 17:41:13}
set {macrosdata(Scalar 2 vector,PrePost)} post
set {macrosdata(Scalar 2 vector,Accelerators)} {}
set {macrosdata(Scalar 2 vector,Active)} 1
set {macrosdata(Scalar 2 vector,ModificationDate)} {2011-11-04 11:03:51}

proc {Scalar 2 vector} {} {
    ::PostUtils::ResultScalarToVector
}

#[_ "View - XZ plane"]
set {macrosdata(-XZ View,Accelerators)} {}
set {macrosdata(-XZ View,Icon)} {view_rotate_-xz.png imported_images view_rotate_-xz.png themed_image}
set {macrosdata(-XZ View,PrePost)} prepost
set {macrosdata(-XZ View,Active)} 1
set {macrosdata(-XZ View,Group)} View
set {macrosdata(-XZ View,InToolbar)} 0
set {macrosdata(-XZ View,Number)} 65
set {macrosdata(-XZ View,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(-XZ View,CreationDate)} {2006-02-16 20:38:30}
set {macrosdata(-XZ View,Description)} {View - XZ plane}
set {macrosdata(-XZ View,IsStandard)} 1

proc {-XZ View} {} {
    GiD_Process 'Rotate Angle 90 0
}

#[_ "Pan dynamic\n    Useful to assign one keyboard accelerator"]
set {macrosdata(Pan dynamic,IsStandard)} 1
set {macrosdata(Pan dynamic,InToolbar)} 0
set {macrosdata(Pan dynamic,PrePost)} prepost
set {macrosdata(Pan dynamic,ModificationDate)} {2006-03-01 12:33:53}
set {macrosdata(Pan dynamic,Group)} View
set {macrosdata(Pan dynamic,CreationDate)} {2006-03-01 12:27:08}
set {macrosdata(Pan dynamic,Description)} {Pan dynamic
    Useful to assign one keyboard accelerator}
set {macrosdata(Pan dynamic,Active)} 1
set {macrosdata(Pan dynamic,Icon)} {pan.png imported_images pan.png themed_image}
set {macrosdata(Pan dynamic,Accelerators)} F10
set {macrosdata(Pan dynamic,Number)} 66

proc {Pan dynamic} {} {
    GiD_Process 'DynamicPan
}

#[_ "Previous view\n    Useful to assign one keyboard accelerator"]
set {macrosdata(Previous view,PrePost)} prepost
set {macrosdata(Previous view,InToolbar)} 0
set {macrosdata(Previous view,CreationDate)} {2006-03-01 12:30:27}
set {macrosdata(Previous view,IsStandard)} 1
set {macrosdata(Previous view,Active)} 1
set {macrosdata(Previous view,Accelerators)} F9
set {macrosdata(Previous view,Description)} {Previous view
    Useful to assign one keyboard accelerator}
set {macrosdata(Previous view,ModificationDate)} {2006-03-01 12:34:11}
set {macrosdata(Previous view,Number)} 67
set {macrosdata(Previous view,Icon)} {viewprev.png imported_images viewprev.png themed_image}
set {macrosdata(Previous view,Group)} View

proc {Previous view} {} {
    GiD_Process 'Zoom previous
}

#[_ "Zoom dynamic\n    Useful to assign one keyboard accelerator"]
set {macrosdata(Zoom dynamic,Active)} 1
set {macrosdata(Zoom dynamic,IsStandard)} 1
set {macrosdata(Zoom dynamic,InToolbar)} 0
set {macrosdata(Zoom dynamic,CreationDate)} {2006-03-01 12:35:26}
set {macrosdata(Zoom dynamic,Number)} 68
set {macrosdata(Zoom dynamic,Group)} View
set {macrosdata(Zoom dynamic,PrePost)} prepost
set {macrosdata(Zoom dynamic,Accelerators)} F8
set {macrosdata(Zoom dynamic,Icon)} {zoomin.png imported_images zoomin.png themed_image}
set {macrosdata(Zoom dynamic,Description)} {Zoom dynamic
    Useful to assign one keyboard accelerator}
set {macrosdata(Zoom dynamic,ModificationDate)} {2006-03-01 12:39:22}

proc {Zoom dynamic} {} {
    GiD_Process 'Zoom Dynamic
}

#[_ "Rotate trackball\n    Useful to assign one keyboard accelerator"]
set {macrosdata(Rotate trackball,Icon)} {rotate.png imported_images rotate.png themed_image}
set {macrosdata(Rotate trackball,Group)} View
set {macrosdata(Rotate trackball,Active)} 1
set {macrosdata(Rotate trackball,PrePost)} prepost
set {macrosdata(Rotate trackball,IsStandard)} 1
set {macrosdata(Rotate trackball,ModificationDate)} {2006-03-01 12:34:26}
set {macrosdata(Rotate trackball,CreationDate)} {2006-03-01 12:32:22}
set {macrosdata(Rotate trackball,InToolbar)} 0
set {macrosdata(Rotate trackball,Number)} 69
set {macrosdata(Rotate trackball,Accelerators)} F7
set {macrosdata(Rotate trackball,Description)} {Rotate trackball
    Useful to assign one keyboard accelerator}

proc {Rotate trackball} {} {
    GiD_Process 'Rotate Trackball
}

#[_ "Rotate center\n    Useful to assign one keyboard accelerator"]
set {macrosdata(Rotate center,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Rotate center,Icon)} {rotate.png imported_images rotate.png themed_image}
set {macrosdata(Rotate center,InToolbar)} 0
set {macrosdata(Rotate center,Active)} 1
set {macrosdata(Rotate center,Accelerators)} F6
set {macrosdata(Rotate center,Description)} {Rotate center
    Useful to assign one keyboard accelerator}
set {macrosdata(Rotate center,IsStandard)} 1
set {macrosdata(Rotate center,PrePost)} prepost
set {macrosdata(Rotate center,Number)} 70
set {macrosdata(Rotate center,ModificationDate)} {2006-03-01 12:35:02}
set {macrosdata(Rotate center,Group)} View

proc {Rotate center} {} {
    GiD_Process 'Rotate Center FJoin
}

#[_ "Info Graphical Window"]
set {macrosdata(Info Graphical Window,CreationDate)} {2010-02-12 12:57:44}
set {macrosdata(Info Graphical Window,InToolbar)} 1
set {macrosdata(Info Graphical Window,ModificationDate)} {2011-11-04 11:03:51}
set {macrosdata(Info Graphical Window,Accelerators)} {}
set {macrosdata(Info Graphical Window,Group)} View
set {macrosdata(Info Graphical Window,Description)} {Info Graphical Window}
set {macrosdata(Info Graphical Window,Active)} 1
set {macrosdata(Info Graphical Window,IsStandard)} 1
set {macrosdata(Info Graphical Window,Number)} 71
set {macrosdata(Info Graphical Window,Icon)} {PantallaInfo.png imported_images PantallaInfo.png themed_image}
set {macrosdata(Info Graphical Window,PrePost)} prepost

proc {Info Graphical Window} {} {
    set info [ wm geometry .gid]
    set info2 [.central.s configure]
    set height_list [ lindex $info2 0]
    set width_list [ lindex $info2 1]
    if { [ regexp {([0-9]*)x([0-9]*)} $info trozo top_width top_height] } { 
        set gid_width [ lindex $width_list 4]
        set gid_height [ lindex $height_list 4]
        set width_offset [ expr $top_width - [ lindex $width_list 4]]
        set height_offset [ expr $top_height - [ lindex $height_list 4]]
        WarnWinDirect "${gid_width}x${gid_height}"
    }
}

#[_ "Zoom all\n    Useful to assign one keyboard accelerator"]
set {macrosdata(Zoom frame,IsStandard)} 1
set {macrosdata(Zoom frame,Description)} {Zoom all
    Useful to assign one keyboard accelerator}
set {macrosdata(Zoom frame,InToolbar)} 0
set {macrosdata(Zoom frame,Active)} 1
set {macrosdata(Zoom frame,Number)} 72
set {macrosdata(Zoom frame,CreationDate)} {2006-03-01 12:46:29}
set {macrosdata(Zoom frame,PrePost)} prepost
set {macrosdata(Zoom frame,Accelerators)} {}
set {macrosdata(Zoom frame,ModificationDate)} {2006-03-01 12:47:33}
set {macrosdata(Zoom frame,Group)} View
set {macrosdata(Zoom frame,Icon)} {zframe.png imported_images zframe.png themed_image}

proc {Zoom frame} {} {
    GiD_Process 'Zoom Frame
}

#[_ "Dimension distance between two points"]
set {macrosdata(Dimension distance,Active)} 1
set {macrosdata(Dimension distance,InToolbar)} 1
set {macrosdata(Dimension distance,CreationDate)} {2004-11-18 21:53:33}
set {macrosdata(Dimension distance,Description)} {Dimension distance between two points}
set {macrosdata(Dimension distance,Group)} Dimension
set {macrosdata(Dimension distance,Number)} 73
set {macrosdata(Dimension distance,Accelerators)} {}
set {macrosdata(Dimension distance,PrePost)} prepost
set {macrosdata(Dimension distance,Icon)} {dimension_dist.png imported_images dimension_dist.png themed_image}
set {macrosdata(Dimension distance,IsStandard)} 1
set {macrosdata(Dimension distance,ModificationDate)} {2011-11-04 11:03:51}

proc {Dimension distance} {} {
    GiD_Process Mescape Utilities Dimension Create Distance
}

#[_ "Redraw\n    Useful to assign one keyboard accelerator"]
set macrosdata(Redraw,InToolbar) 0
set macrosdata(Redraw,Accelerators) Ctrl-less
set macrosdata(Redraw,PrePost) prepost
set macrosdata(Redraw,ModificationDate) {2006-03-01 12:42:09}
set macrosdata(Redraw,Icon) {redraw.png imported_images redraw.png themed_image}
set macrosdata(Redraw,Active) 1
set macrosdata(Redraw,Description) {Redraw
    Useful to assign one keyboard accelerator}
set macrosdata(Redraw,Group) View
set macrosdata(Redraw,Number) 74
set macrosdata(Redraw,IsStandard) 1
set macrosdata(Redraw,CreationDate) {2006-03-01 12:41:13}

proc Redraw {} {
    GiD_Redraw
}



set ::icon_chooser::ImportedDates(rotate.png) {2006-03-01 12:32:22}

image create photo ::icon_chooser::images::rotate.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAA29J
REFUSImllT+IG1cQxn86FjyFii1cbGHIFgGrcHEqDPvAArnKpTjwQYJzkJDrfCFFAm4uVQi4
uECKGBK4QIpLESJDCgccnKts44AcSKFACh242ICKDah4ARVzIJgU70m30p3s/BlY0O5o5ps3
3zfzGmbGzBKAE+eMHjTMjN9+dAbxZWYNXxQml2A/B8wMPyhsY6MwzMIzKAp73Ctsd6cwMwvx
Jx1nXAQygTtK73Ap8faWs7KKIHVH3dYA7jlnfzpnR984W3DcANIEum/D8VFwrgGQAi0BMkRq
EQcFMFT8rxXPfoe/nLM5uHPO8gyo4Lt+v3FuVUfO2U9ADnzQ7zfqvnnAL85ZFssQgO8BhaqC
ZwO4+X4IXJtFrieQQSg2F/gihTwla0GxDve+DKdL5ljNGq5XuAN40CqglKNA7pkzfPqRs9ar
AUljSY8HsD+CCpizuurpFoHl2bOSu1WWLH+423E2iC39+OliSxfaCvCk4ywHhruQCpQVrPfg
ci1w3tbjjrN1Qmu7XVg/gBtd0D14ffNUTfOS8vqXHeCSQK7kwIdvnJa0Rt0SoClwJYWvMxBB
BNIUtt+sEXdy3RkigbwUmEIQiMBUAUib9ZJEwq8pMAbGCtMKxoofQzUODZgHXHj4qAHw5NvA
croH+lzxPvxxWML9ydKhAfZ7sLMB2S7oAPwEyhLu/xzG7PIyDzPb2nSWNkMpBwq7wMMZFy/S
0VvXCiuu/Q8tta+2LUuElDCu/yTmpQA/dJylqmSJcEhYL6kqaSKBqAReO0ejKwG+ckGvWQKy
lMh/pogITEGn4CdKNRbKEbx7+3yQMwB/dJyJKhITz/WSAA+AmWx7gA++aqRUXhiUcOv2iiU2
sxMXB6WeGIJwU8KG/zxuWXwUKlQeylIZjoTDB9CPbTsz/wuJZ7OUKFyMoyJhk4fyBQSUMC4i
QpZCK9eFdOH+EUBjpcREUz0FUsLy7KXE7QcThUnYjaqREy/oRDi57uzCo34jjNqSCrY229bT
eJKphuRJrMpXMAk+jcmrseK9UFbCsIThSJldVi+U6d1P2tbKhbQZVkw40aKKVIUqcjAYKgcj
QQFN4JWnK662ZXvvnbZlWQCSBLwqTAUfWzN4rvQmgo+JK1XKRLi16u78N3bsnO0BVXKqnLr9
J4Dtq23zSai4/5KV8TcCL/io3Y9CzwAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(sendto.png) {2013-10-30 16:53:57}

image create photo ::icon_chooser::images::sendto.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAABLhJ
REFUSImVlH1I1VcYxz+/X1dn2nYNsld2u+YtFyg6y0qmzV6UgouzZRaNQlNzOXRjbNe5lCtt
WXNiDUuKyI3IkpEoyZrF0iiVujasxO6aL6vVpbFeZqnkW/fZHxev3q6u9YUDh3Oe83zP93ue
8ygiwgg0AENpaQKgAizRaOj5ZhbK2DDln6Qkifb05LvNm2FzQ4MMDn4tIgIigoiQEBYmI3MR
cZwfGtotd+928+aebobsdryPHlWciVevWSMz8vIAyNFqXRnHQgXo6uqSdwIDJTQ9fTQqLCxM
BlNTZTA1VTo6OiTRYBARQbNsmSfs17NucRlVe6Hx4UNHqrOX+kiyRLM2cwmU6lGmaBxGdNy4
oQC8Z7NJk2UOk8MsDukT3WoiaEYmwYGBciUyEg9VBSBcVak0mejt7SUkJERxmjvCEB8fLjbr
PZqWG53ZwlWVZrsdgEcDA+zo6RllaLn+kIyPFsPHegDq67vwyvyJJ40paLVe3Gq6Q/ZrGaMH
FkVEcGHqKm5Y3gKg7Yc2ht/wJcsai6KqPPr9JvnB/yE6IipKiIoibuZMcrKy3DX8X6ivFM0Y
WwGMq1fLRrudRl9fLvv5AXDt8GFlbIwKUFNTI7sjIqTK359NAQFc9vOj0mQa/0rt7e3SvGsX
pqAglw3TzZsOhoQEiV20yClUERE2bDBKcPAAgZfh/Tl6wPFogPPhttTXU9HRoWiio6Pl1Kko
tFovAJI3/chnynyap04F4GBbG9dmz6aio0NxMgDoF+ql1ZKMl5fDh5XLj/DukrmYi2Id+/rv
sNn+VpQzZ87I9xoNHt7eDPf3Y/lkO7d+TUZRHOYUFV3CYoGqqtpRhs7OTjFu20ZYQQEAA0+e
YDV/jiFhK1MiIwE4nZhIj82muLx0S0uLpJvNzM/Odq71P37Mnf37uXr+vINybFcYGbW1tRKT
nCy6oCB5ce+Va+lV8dLaq6mpkeLiYrdbhKanu/afCTCuArPZLKfLy9kTHMyqadOc6+fu3WOn
TgdApcnE+sJCt1p+Ec7PYDQaZbLNxrcLF5Lr40PuihUugeGqCjodlSaTs4wB1gUEiH3GDL7c
t4+lS5e6kbkoKCsrk715eWQaDKQaDGhUdwdHSn5EwUjpA3Q9fUrB9evc8fBgZ3ExMTExigpQ
WloqQUHzRKu9QNvtDD78JZanBbPIfv1PVl48x8X7951Jmu12mu121hcWuhD3DQ9T0dmJ9dkz
MtLSiImJcf06J06ckOLir4iL05GVtQwfH0+XBK2tf5FrOotfvwf5uhDm+PhQ3t7OwT+sJH4Q
QnxSMPn5dQwM+FJUVIa/v78rQV1dnWTm5DA3JQXtggW0V1fz6HwNmdtD2JH6NhqNq112u6Cq
o5YPDAxz7FgLZWWtxMaupaDggCvBppQUudrWhv+WLUwPDXVJNtTbi/XkSXpbGsn9IpKtCQYe
POijpKSJ6p9vozfG4Rm6nNYjR5ip0VB+4IC7grGwWq3yaU4Otr4+dElJaPV6l315/hxl0iS3
c92dnVw7dIigefM4e/z4xAQvoqGhQXaazdj9/Ji2cSPe06cD8HxwkN8qKui+coX8/HyS4uPd
/8R4veVlo7q6WkpKStz6znjjX6P4cngykriEAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(view_assigned_meshsizes.png) {2006-02-16 16:52:44}

image create photo ::icon_chooser::images::view_assigned_meshsizes.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAO1J
REFUSIm1VEsWgyAMTKQHtifhwvZNF8U0PyhWnY2PCWQyQUKkwO0LquW7AEAWmuIQb4QhF5dB
UoQAM7PnAkL+mKoWs0FOAIAOmoD2kmqwCmZ8gBQFgGInBwcy+BrHB7RL51on+xitBRmoFhw2
fRiZwjD7QjTfUjnQheuMwHTI8fMKe3anMi6JiLBuZs1auj0Tw8nGFswVMrONSxUys7vKoxdI
lemPm74dWj21qjv9D37e81mcE6glTIYR8imRjMo9uZ4k1wokI2tGZE6gMw97IpdeMtatP5Om
HAyq7zmZfwcH/hbB88W3P7Tb8QYMaWMJbXpHwwAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(collapse_points.png) {2006-02-14 23:37:21}

image create photo ::icon_chooser::images::collapse_points.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAcFJ
REFUSInNlDFuwkAQRd+i3CFN0iAhupzCDTW5AQUH4AxEgj4hRVJQW+IY0ASJSCC5SotESW3n
pzBrdm0vMakykmXP7MyfPzPrgbJICJApdFV8QKORBKglwcMDAMYFAIo3LYDVSgwGz0qSdQ3i
NaLAd558MskpO0Sg25Uk5BmLIKEoevMOKjWECDSx52RsayvkzsFn2/Eo9XpSmkr395Ix4QIK
hP1+qzTNJKHxWILa6VXTBehcLCxc7F/kejS3nSXbhRg/qAJg9dZJMSdUOw/SVAXQCcwA3Njw
6VQcDtDvw9cXDIcvSGDMGdATd7p2gDW0/CKTZC0JxbGUZdJstgtP2nbFcm63N1os5vVZSi0s
O1VbXo4PsQjW8y/EvRlNqDb18wbuLdPghTj5XTgPJqgkqUkWAneXgXfewhdjajjFsRiPpSwT
xojZbGfJ2MclY9zk5QTEsXh6Et/fOdjr647HR8PdHXQ6n0TRO7e3H8V2M8YD/bVlAMznUru9
URS9abGYF0yPR2kyybfvaJRv4EaAJXFXhFx9ucxXereb//iDwXPj2+SB1wR59v1+qyRZFxsl
EFMFv+DUFCTo0/Tn+n/yA1ejdq3eAZS/AAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(surface_parallel_lines.png) {2006-02-14 21:44:20}

image create photo ::icon_chooser::images::surface_parallel_lines.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAOZJ
REFUSIntVUsWAyEIC71YPbo9WboYP6hA7Uy76WvYjBqQwaDAAoIAqQYGR1ECWoNsscaY58Bi
KpSAIIQARacka1YCxQCQc76c0ZAdocyKLp1YPufE61zJ/1ZmOoNSCGfS+lFM9VPlPVande0w
6GMuUC/7a4k353ZwjtwrkQJQRhEeg1FHo+jNoOwEo0N8vEXWTn9ULNJa4RVXNnyjoBu2c8PF
zdNYSrUG8iMDAFJKQbNZwa0LcN5UbRxsxON355bx2iiy1mrO2bl9GvauD1cJ9QBTuttlUm/F
ZQzl+8BrHeCrwc/hCbHT9WKMywW/AAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(view_rotate_-xz.png) {2006-02-16 20:38:30}

image create photo ::icon_chooser::images::view_rotate_-xz.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAbZJ
REFUSImtVbuRwjAU3HdDLmISqwR3YHVACdABJUgd2B24BGekdgemA4iI3cFeABLyb2xztzMa
NHor7fv5ISTh8eM3IkKJLTJJc871aSAJay1JIkkSksTlcnmxRIQkBRH6kiLhrf67wxvOORZF
QQC43W4M4iRxPB6plKJ3RF62vvCyxtCjnka8siwjSSilmGUZAfAt+blQ1zXruiZJ5HlOpRS9
zWeK5LxPc5j3dQYjham6jBS24LsYYqRpSgAoioKPx4M+PmPMp9Z1XdNaG+oMgKfTKaQyTvEO
ALTWMMZMumCMYdd1IQn+gmit8d6TpJzPZ16vVwKvNnz/ymyLrw56CYJXY/kELKr8eNLXLi31
VAh6yrjqk1vK2KpeapqGxhgZnmmt0XUduq7r8ff7PdI0Fe/m4mrblgBoraXfl2XJIa8sy5Gt
R7DWhvmAaE7EowAA7/f76PEkSehnTbw2RZDnOauqIgBWVRUmtFKKbduOHh/NpKmCNU3Dpmng
nOudO+d4OBzwfD4n6+b5QSBupbVNuwZhSMaHawfmGiym6K/YPC22YrdEGKZra4S9CP47PcAX
fztb8QuXMfEVroo9KwAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(view_rotate_-yz.png) {2006-02-16 20:38:46}

image create photo ::icon_chooser::images::view_rotate_-yz.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAU1J
REFUSIm1lN2RhSAMhY939l07kBasAK3AUijBsRLswBIogRLUDrSC7IMXB1jwd/1mHB1ISHIM
ASJQbANIAKBpmoiJ2RFCOO/gqUlgI/GNwjHyPCcA0FqHM7FT3ClnBzr0VEqRlNIYEefcONEf
Q6XUtmirYip5xHGuB85RPndPfQQBq2p932+SfuVdUUpR0zRk/nNd12RLmabp9v0BAMYYyrJE
WZYAgHmeMU1TAgCcc1qWZb/3THitNXHOnehO3me4LGuCy/dk5VFKu86h6xva3+VUhCjDMNA4
jqiqyrG1x8I8z45PlmUoiuK0gBBCkN1OeZ6TfYMN31FgjwSXb7f4D4C1n82aPSrsoM4FuINp
W3vNVBebpL5OFFjbkFISY2z7H8MwUNd1Qdu2bRM/wK0OP2KvT/8tiM3zee7x+hD+OWHjV3VJ
Or+CV3R/lV91ja+wIf5MbwAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(point_move.png) {2006-02-14 22:13:38}

image create photo ::icon_chooser::images::point_move.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAALhJ
REFUSInVVEkSgDAIi/7/ST30Z3hwRJbSUqsHc+oMEJYGAItaKwHA7iwtN+p7AUAphb2EO6nQ
BM8SCLhbm4v6DJZdjV0N9yyccE6NZspKOU8x/gncV6TmpIYRiscRkNT7/WE8Yh+QSt8v6X1J
sOjSzD/V0LYYb9t2fLMJJKGN7dnGEDc7r/sRod2L+MjyGrgCZCvUfCpsQeuXViP7AP3KF+HJ
Xzy6mvwZcWtm5i8ezjWXYHkRv8cBxAV3qvbsi/8AAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(Scalar2Vector.png) {2007-02-28 17:41:13}

image create photo ::icon_chooser::images::Scalar2Vector.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAJZJ
REFUSInVVNEOgCAIPPi/vsn5TX4gPbgmlhUWzrrNB93tgAMB3CEi4q9aR8hnKmwJsL6EEPyz
pvZzUpGWE87XYO8rb3TAZisDAI3xIVvNd7RuxT0a2e+J5k63I3STPIueBOVY2VxlUt7PTOUR
US3psRIOGWrRGOODCpLoWXL4NlejNGQ9mQbcT2jmR/ApVTXhv80aixWlwjcPS0HumAAAAABJ
RU5ErkJggg==
}

set ::icon_chooser::ImportedDates(collapse_lines.png) {2006-02-14 23:42:43}

image create photo ::icon_chooser::images::collapse_lines.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAdVJ
REFUSInNlL1q40AUhb/xD2swBFV+jzguTCDdgsmbbJnOT+A2j5EmLhJCfrrAVm70BPYLBBbL
qWyBpLOFLGUkzyjOLix7QWJ05+jcc+ZeCQ5jK0D1xEcYAOkDIvuBLMuqeDty5K5eoRpRFPk3
AUiSRHUVLYA4biMhY6qbANrtVJHa2VMAhvVagFFpsHgLYDabVSR/Is8Ri8WiYLCZ60QHxIpj
CbaSKC+vAhv0/p6/2AiulRbA25sEaLVafd3nX4Wv3EG+VSySJFFt4D7VrKurH6XZ+pDg6IvG
Y+n+PhXVWff3QkKjUayHh7TeOI+mPWA4jPX0lDrBHZfJMPzGcBjTbqdA289eGCvWo1Gslxc5
JVVn3ZJydvbdadg7xvP5XPuPxwv25f592A09RtFRuJa1NsZAGIYEQUAQBJpOp96hd/emuQCA
6Xa7DAYDTk4GXF/3uL3N6qMqCYzJ8Q5Hzb8Zi0RRJF1cxIKt7u4Oxr0kq+eO82YVWq+l8/NY
sNPjY7XQMWobye1rs8k0HueOnp/TPyaukOM4is1GOj3NHb2+Zo2FOo6cs4nG2Pn8njtKgJ/c
3PxSrweTyYR+v28cvDl5gxrveV9elv8HLZdLL8exH9f/Fb8BmhCcmGFnBLMAAAAASUVORK5C
YII=
}

set ::icon_chooser::ImportedDates(split_elements.png) {2006-01-31 12:53:02}

image create photo ::icon_chooser::images::split_elements.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAM1J
REFUSInlVVsSgjAMTFoOjCepB8aJH1LaPEiDijPqftGy3Q3NAwAJIqK2KpkUQ/P5GRclU0D0
JXTREABAqsZ4nRgRBRsQEc+MbI0jfgEl03a//ND2nBh5XrTbA6heMOV9p+/A0ZhHJZvYqqpr
l1Y7e+qiqKhy26YV++WG8gBX7xLmf4skrmvLVTWEInygQf4PTkbeJn6wuQm60g6LOyauYDJ3
nUhpXkYtzwwnh2hC/ipADdARxEC1WvD5pBv33ldRIB8K/jiKoo2tH8Qd9bAMX9E88XIAAAAA
SUVORK5CYII=
}

set ::icon_chooser::ImportedDates(rotate.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::rotate.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAA29J
REFUSImllT+IG1cQxn86FjyFii1cbGHIFgGrcHEqDPvAArnKpTjwQYJzkJDrfCFFAm4uVQi4
uECKGBK4QIpLESJDCgccnKts44AcSKFACh242ICKDah4ARVzIJgU70m30p3s/BlY0O5o5ps3
3zfzGmbGzBKAE+eMHjTMjN9+dAbxZWYNXxQml2A/B8wMPyhsY6MwzMIzKAp73Ctsd6cwMwvx
Jx1nXAQygTtK73Ap8faWs7KKIHVH3dYA7jlnfzpnR984W3DcANIEum/D8VFwrgGQAi0BMkRq
EQcFMFT8rxXPfoe/nLM5uHPO8gyo4Lt+v3FuVUfO2U9ADnzQ7zfqvnnAL85ZFssQgO8BhaqC
ZwO4+X4IXJtFrieQQSg2F/gihTwla0GxDve+DKdL5ljNGq5XuAN40CqglKNA7pkzfPqRs9ar
AUljSY8HsD+CCpizuurpFoHl2bOSu1WWLH+423E2iC39+OliSxfaCvCk4ywHhruQCpQVrPfg
ci1w3tbjjrN1Qmu7XVg/gBtd0D14ffNUTfOS8vqXHeCSQK7kwIdvnJa0Rt0SoClwJYWvMxBB
BNIUtt+sEXdy3RkigbwUmEIQiMBUAUib9ZJEwq8pMAbGCtMKxoofQzUODZgHXHj4qAHw5NvA
croH+lzxPvxxWML9ydKhAfZ7sLMB2S7oAPwEyhLu/xzG7PIyDzPb2nSWNkMpBwq7wMMZFy/S
0VvXCiuu/Q8tta+2LUuElDCu/yTmpQA/dJylqmSJcEhYL6kqaSKBqAReO0ejKwG+ckGvWQKy
lMh/pogITEGn4CdKNRbKEbx7+3yQMwB/dJyJKhITz/WSAA+AmWx7gA++aqRUXhiUcOv2iiU2
sxMXB6WeGIJwU8KG/zxuWXwUKlQeylIZjoTDB9CPbTsz/wuJZ7OUKFyMoyJhk4fyBQSUMC4i
QpZCK9eFdOH+EUBjpcREUz0FUsLy7KXE7QcThUnYjaqREy/oRDi57uzCo34jjNqSCrY229bT
eJKphuRJrMpXMAk+jcmrseK9UFbCsIThSJldVi+U6d1P2tbKhbQZVkw40aKKVIUqcjAYKgcj
QQFN4JWnK662ZXvvnbZlWQCSBLwqTAUfWzN4rvQmgo+JK1XKRLi16u78N3bsnO0BVXKqnLr9
J4Dtq23zSai4/5KV8TcCL/io3Y9CzwAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(view_rotate_-xy.png) {2006-02-16 20:37:56}

image create photo ::icon_chooser::images::view_rotate_-xy.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAcFJ
REFUSImtVcutgzAQnI24kw5CCRQQKe4g6YASwjkX04HphJxyxR3gDpJCIs075NnPfPKAKCOt
BN7Fs7uzNkISHhtE6L3IsjCQhNaaJLHb7UgSaZpSSEJESFImPxcRikhg6TGOvqiqinVdEwCU
UgzkJHE8HpmmKX0i8vL1iec5hhn1OGI7HA6h0F8qdF1HY8xfUNu2bNuWJGGMYZqm9L7z+Rye
3+b0Du9zfYMRw5QuI4Y1+KyGGHmeEwDquubtdmNVVQReqQJAAgDWWlprAQDOORERFkWBsizl
crl4LSQwZFkGpRSUUqMUnHNyvV7/FoYq+3EuioL3+52x8kG4uVb+W/QcBEDo6xKWjQ/6OKW5
mQpFTzkXHbm5ji2aJWstlVIyXEuSBM/nE0Ofc44AkOe5LJJhu91CRFhVFZ1zFBE+Hg/s93ux
1kJEGG3KsiyR5/mLNJ4irTXxkjHYcMoAhAnz1nVdiG+apuebvJSG5jcwxrBpmsmNtNaMR9pb
0hNkQjB/EuN1kvg9tjydTgIASilkWTZqbxA5HqWlQ7sE4ZIcVvItgt6YrrkFlmL1bbEWyVzA
sF1rK+xV8O32AB/8dtbiB6mnw2kPaDLSAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(view_rotate_xz.png) {2006-02-16 20:36:44}

image create photo ::icon_chooser::images::view_rotate_xz.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAZFJ
REFUSImtVcGxgjAQfet4h7MX04F0IJ1gB9oB6YASPHvSDrAD6CA3z+ng/YMmkxD4gOObyUxI
Nnn7djeLkITDFgCKoiAASLgjiZnWOjUDST+OxyPdXEhCREhSEGDjGUQoIv6u+N6pEwCQ5zmj
DZJyu914uVzehu+1mHieY+hRIjAUeb1e2XUdAfhBMj3woaQ7WFUVjTFxpNZg2tcJJAxjeUkY
1uA7DSGUUnSulWXJRN8wpPf7nfv93ofRGMOqqvz3dozWWuvnRVHAWiuTDIfDgS55AFjXNeu6
ni/xxaLnIB9q594sy8YZfe3SXE150WObi57cXMRma+n5fFIpBaWUv6Tve+Z5DmttlGSHsiyn
Ez0cxhhmWcbz+cyu65hlWZR4N9q29YURvfY5grDCAETtINxzFZi0kzUKnJdN05AkmqYhALZt
O3p50pOGCev7no/HA6fTKcqB1pq73Q6v12s0b1prb+sJwlJaWrRL4JtkuLi0YS7BvyH6BVZ3
i7UY7ZAhhuFaqzBS8OvwAF/8dtbiDzIa/NPxk2f3AAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(surface_divide_nparts.png) {2006-02-15 23:30:30}

image create photo ::icon_chooser::images::surface_divide_nparts.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAQVJ
REFUSInNVUsWhDAIS+Z5/ytnFv3Ir6N242RTpaSBgggYEIAkgWwvGAYPQRKamdOHiDiJtFSo
WRgoUeYSkmVxxjqjJwCQKZEmzhyVcd+KKjmzG1twi/NoQrHq3V+Z1BOIBEOUPcGRKkKoS8zg
TbgaDxyF36hfeUqvg5yr3e0NM9nH0Hbi7nY4dlMd1MuWIqfpSBdo+so4lrMzYmYaS1QaKh8U
UHGhP+4uflIhvkpgufOX2O3v0BZbAmfN8wNm09qmWQkse+EW8jBev7SoT5NQjf5qZN0W8Pxy
AK6wyOSiyGE00y21yjOBjWzCTLrZpuufTKlhmuqOgMrHftQV56HA6z+ejC/DEX7/GQY/swAA
AABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(send_to_back.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::send_to_back.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAnJJ
REFUSIntlDFoE1EYx/9Xr3DgpdxwlRsidMjQwaGDpBkKXTo6JtLBQukgHXTL4BDQIaNIhiKl
QpGinSqU2iGDLuFwOjUUKQGv5ZAbTg36KKk8mtO/Q5qY3MWkgh0Ef9Pdu/u///d97/seEKWw
tkt1LEkAgJ7OkmTkF5IwKy5JYgQAUK31fiUJdURl5xkAEtM5QtOgrZSRmM5ROtsx9+EobY+e
RUU5jQTAwddjXn7+uvPHCAAkZhc4b9cg3QCJ2YVfcpLY+lCnuVTsk+pfo10Jkihu2jQ3bapj
SXavd4Jtc/vRSz57e8ASACkkjO0yxu8+5ejF8Xig7V2uXr9F51OD5oNdmjfvxxzUqLDmBZgv
rAKhhPS8c0j+fOhOjGSnWaMJ9xXo6SzdL43OsQ8U6Oksw5OQ+T2fM/t+jygm0NNZ+kcN3qi4
1NdtJnecHlEbBWgNGHQDCCWmFpdRBWClDAT3SkAoAQCNyhOlbw7+t5CFVy7n9l2ql1KxkHp6
CQCs07OvBrJvRWOtMTV/B4FhtF5kf9F//oBof0ZbMv+wNbTRq2FgXw8z0NNZFtZ2GZ6EfHFY
Z+qNS3PDHmh0JoPOxj9C+kcN5vd8Jncc6us2zYrLa4d1JjdsmkvFvpddGyVqkJjO0TIMTE5M
wgslAsOAMZNBICQ0KSEBaBMWUPMA3wOEAISArAeQThnN48+Koii/N4gyemGUq44HqWt4/M5D
TUhkMhaqK2WIrRKaH98rUU23QWzUokzMLWLxioWMoWE5ZUETYpikh9hoRgmEQFVIrDoeaoGA
1LTWnMuzGZ2pRNANQNW6wtJa90IjQPN7c2CJ/n1+Al6vBVnj2YnXAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(view_rotate_yz.png) {2006-02-16 20:37:05}

image create photo ::icon_chooser::images::view_rotate_yz.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAS1J
REFUSIm1lFuxhDAMQMPO/QcnkQAo6aAAHHRQgAUcIAEJkdBKwEHux1KmlNAt7HK+GPJ+NAAS
iMiiYI/WOqKW5/lBeviRCYIsVJIpy5L9OC8nqKoKiIi7rjsPnASnWDIAABHxMAzOIG7Ytu0m
VEqxMeZefpdzjRmf8ooJfw8RsbdZvM45Xh8i+kZgjGGlVLwh/pqH2yoWXRSF++Su68QHkDyD
y23N4NY7+TKlqLH0fCV5lKQIIsYYttZCXdc7PSLanC7LcrAL9aOsQ2bnFBHdTu6Y55kBgONX
8YTVOQMAT9N0cICInHiVz9FaHzJ3R2vNXiTsFQv/AODdAmstNE2TAbznM46j6LTv+82H7+zW
hn8itqc/C+Lz/T0PePwI/yXohFVdal1YwSN9f5R/JlSSZrg+sl4AAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates() {2017-01-20 13:02:17}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(line_divide_nearpoint.png) {2006-02-14 22:22:00}

image create photo ::icon_chooser::images::line_divide_nearpoint.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAMtJ
REFUSIm1VUkSwyAME5k+jWc7b1MPlGVcMCYQnZLxIoEwAH2Qg4AvV0QIQAfGXamSfPT88Sxh
ZWm1QlWJSCGnV24qKOvjUz2WTus/4aqf4RRzpXxU9acjb+ktAgCIMfa1tj5Y9GXRupPjBPkM
uuzwOZe3MZfinsCPj6nsuPvAjRVkw1rjWGeUo5wl6AYtwXZzTdJ8l+uhXhVzt4xNJYEQtNL+
NKVcv/zdN+Jc0WvNt4bWKl5vPLle9qGdpxHLKXm8Xafm9GPjEHgYX9ftnlq5oh7QAAAAAElF
TkSuQmCC
}

set ::icon_chooser::ImportedDates() {2017-03-10 14:32:10}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates() {2017-02-24 15:17:10}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(view_rotate_xy.png) {2006-02-16 20:33:42}

image create photo ::icon_chooser::images::view_rotate_xy.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAaZJ
REFUSImtVdGNgzAMfUb3n25ANqAbwAbtCN2AEcgGHaEjNBtQdYFkAzpBlQ3effTIQQMHVPek
SIAdPz/HMUISPTIA2O/3BAAZWiRxM8akbiAZV1mWPBwOdM5RSEJESFIwQBYZRCgiMdY47twO
ALDWMiF3zrEsS5KEvGxj4mWO94wSjl4kSVwuF97vd9Z1TZJQSpFkuuGHMiZHEkVRxOfZnOYw
n+sMEoapc0kYtuAzDUNordmnttvtCACPx4N91yUlvV6vzPM8ljHPc77cXu+TGkII8fl0OmFY
hGSDMQYhBJkr82yLrxa9BAEQqdewZL3Txykt9VQUPWVcdeWWKrbYS7fbjVpraK1jEO89syzD
8/lEVVWj4N57hhB+v0/d0OHquo5KKdZ1TecclVJsmoZ92wFg27bsx8HwNs+OgKlVFAUBsOs6
vtuUUgTA8/mc2DYpaNt2MlDTNKOZNFxffx2Y957WWnjv4xmQhDGG1loej0cBgKqqoLX++5CH
rbS2adcgDsl3Jf9FMP5bbpgCa7F5WmzF15LDe7m2Khwp+O/yAB/8drbiG0cZ9bG+7M/EAAAA
AElFTkSuQmCC
}

set ::icon_chooser::ImportedDates() {2017-03-10 14:20:00}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(intersect_lines.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::intersect_lines.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAANNJ
REFUSInVVEsWxCAIi31z1umZPC2zqLYKUdBVJ5vWDxBjEKDIInxhE8lkPVPSm0R9faiILPLU
0CfoFhfzr7GKJNbZ7vGjiiZ7AoBVbY7AqQX2cGzu7dCcu/GHBCxKeSUtclpJrwra4hnU6gBw
3JtYlSFqhZZG3NvtJj+AbejmjkDJlyPqgvbcS85h1tvBkEDPZvrIivAnpMZ4z4s2GTVd2yPl
3zGor2cTLN/KM/4YWolaabpxkWjEdtDMdpLeA9N/dCfLmPWz3+szBlHvb/fIf+AHIeOE+NyV
kIQAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(redraw.png) {2006-03-01 12:41:13}

image create photo ::icon_chooser::images::redraw.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAetJ
REFUSIm1VTtu20AQfUsLJNym0x1csUwK7yKCbOgIPkPUJJVOoCqV76AjGIZhYKEiLSvdQTeI
QALSuCAeObv8WJbhAQRxl7Pz3sy8WRoRAW0CAM45AQCj3wSLRLvh5/9fIiL4/vtOGrfZYSnH
xQ5VniEtyvbIcbGD996kRYkqz+oX3ntzWk0BAKfVFGlRhojaEj7MDkv58ee+8Zrw4bjYIdVH
RKShSaoiAmOtFe+90QlVeVaDc4PMtvON6bByzsnV002Dy3Q6bGPn1+tHQzidRQAdR9OB+K6p
ICNrB1qVZ6jyDM45Cep0Wk1x+/Ig5Hn78iAAMLn7Vm8U+7p+Qz9rrcR7g70bskm8MTsshRSr
PMO/v89BETpljQNwj//BAd0wlpfV4To4oDvrnJO0KJEWZdAbY60dbRrp8F0wLEArCRplQaRg
7l6vH02ceFqUYC5NDldPN73K1FS4HtS9c064pm3nm1Z8fcnSkvUe2/mm9vkSLfUpgPYeu3cB
9NQCrRiAttKxGM4G6AuubyQCEISmwTrTc0nwvsCNjTVtrHnn+gYlGmM9linQ9oPl6p3moQAf
VRGvDEA1+VL2Y33pfA8JoqXYx5Zjm6z3Y/jw3vd/2visZz9Z75uLML4TYqMvcMYcfCSDi2X6
Gfm+AU6xR5fuyTgjAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(dimension_vertex.png) {2006-02-14 17:51:53}

image create photo ::icon_chooser::images::dimension_vertex.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAYdJ
REFUSInlk8GK6jAYhU+KiBKle4vtyoXFR/EFXPlevoEvIbhSBJeCdmlLQURR0WBNmnNXd+HM
3Ll1xmHmcr9l8uecPzl/gGciAGA4HDIMQ7RaLaRpisVicVfEt06+ufhnj+VyyTiO4XkeXNdF
FEXFJIv7PBljDAEgiiJuNhtaaxnHMbXWVErxeDz+tbfCzbNo8bc9xycRADCbzWitRa1WQ5Ik
MMZgt9uh3++LwWDAbreL8XiMXq8nBABcLhdKKUUcx7TWolKpYLVaIQxDHA4HBEGA7XaLRqMh
vuFSo9Ho4TweOvD8+XEesf+Z3CWtteb1eoWUEo7jCGMMtdYgCcdxQBLValUopVgqlVAul0We
5zyfz6jX67jdbsjzHAAgpRSvDKbTKUmi2WwiTVMopeD7PpRScF0X6/Uavu9jMpmg0+lASon9
fg8pJfI8RxAEmM/n8DwPWZah3W6/P6rWWmqt79LMsoxJkvB0OhVO+iP8W8IEwI/8wcLiXyH8
W/zTvBwjvrP3n/ILe1LRjenZOakAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(dimension_angle.png) {2006-02-14 17:36:29}

image create photo ::icon_chooser::images::dimension_angle.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAIlJ
REFUSIntVFsOgCAM27z/objZ/JEEJpU9MGpiv2BspTQbREGUUkTHTgEYXIgk/6j8bslRhHQJ
KrSwweKXw6NZXA8MMZus1klXRX3yMYuCZhIxiV5vE408Of8iMm9qDYQ8M1OfgK/nXcS1GWnt
7wHHI3uByYrIJWaPTVbV9mqSPKJ4kM9w82OIHeINP9Vq5goTAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(border_horizontal.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::border_horizontal.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAKxJ
REFUSIndlVsKxCAMRW/EBYhk/yssQxYgdT6Kg80kLVIEZ85P6yO5TVIj0EEAICIAgAAXEbG3
nRw0UkrHQgyx9gtlL+RauLgW7gePi1/EPgi1Fy2jaTEPS883WBCzcpy5bq+N+rmF02r+lhZl
L8feNqGD7OHMH4cLxrA4McQqIo+e2idpgavy3cGZq87l9AM9XeArRU8d6hSdBnet2OojGn1F
/X4N/uFamMwbk6RldmBRGXgAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(surface_hole.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::surface_hole.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAUFJ
REFUSImtVUkSwyAMkzs5lP8/ltzUQ4yxwWSbqhdojCVvAMyoJMm+J7gw6V9qsiHY1hJOE4TU
hBnmYmKcFQQOANKkNiECAQSy9fXOtXdzd6JKsq8mKKAAADbnvEMqyG/fyn5wCDICWphHXWoI
e5sPSPNDgUzRnsQwCD23VPMbtbrEgmhMKaBplcHINYNlZs6S1r44Xq0DZxkfbxSM29+yX0Tl
imbT4tKbB60GBHs36GJdjrFt9cAnNU6RpXo0CCm+rvpo8PgAgFs9+BhhFJP58UgS4zsogyy3
ztmZE3c3+o4D2lVR1KYI4O8YvXTSTPr2Tb5bOyfCGkGPxFY1DFVzMqru4RQT09R7ggRqpPdZ
J6qByNLyH3C+T8bf0BgPRu3QGx7AG8ofEtwAC3wULweNx+uePuchKnlDECb1xfn/4gewIdsO
RTEJTgAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(intersect_surfaces.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::intersect_surfaces.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAUdJ
REFUSIm9VcsSwyAIXJz8/y9vDwmyIIltp1MOrYK81sUAAHGJAQBBmm8AkHFAViJk2hCkKAkA
I+KZFY/HcGvWEacILzgvZB3/sUnhvFTeJa9FpObUyFs7i4ITpbW9o6Z0+A02IfQ7wqWjN5Iz
9c3Miy1ryeWGmfEtqbg+H1AZZ+InWrBDyQQJb8gBMShKF+m0ITKanWoNuIdS6xu6OX8qHVdq
sKNGMDLbRvKr+Rt0DXKzriKoXEoi5JsH+BEHSlE+Fzv0spc6SA/YVcJ2uWGxKRNylup4Mqhe
Wfjp0xa+po6KvAbyNtdSbdlEnETpOv258s4Wut7P7Udri0Owc8Ka2vOEalD5xujAzbbUYYGl
mVzVLrUoi+DDcjcFucqaoPdp3xU1eML6dGvCdx7dp4lKtvKR2V72l2Mun1foDPR38AvpZ/Uf
8gIXxAOcBhhwGQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates() {2017-04-08 20:15:17}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(align_quadratic_nodes.png) {2007-05-03 21:01:46}

image create photo ::icon_chooser::images::align_quadratic_nodes.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAANNJ
REFUSInVVdsNwyAMvDzW6BwdotNkkkzTFTpHpc5RXT8SaGtsMCFS05NQJHMc5nAMkIJUgiau
t4fBnwd+fdXJBgiJPmZ0viuE4oY7ZASpYBx8HpjYVueSdkeaakbEOm27C7tCtbRXiMEAGTN0
NVc+YqPYlZi6Lp9pmkpQYxRxkN8LzNsv5O6bsOLGD32wsvghFuPCsBDnxwxJigJY+tcFJ7va
And6Oity7TlArl2nXJ8wSIB662kW17OmGCFM/+NVY0l15usCl/DxxMu1vZX753gB452clGqa
td8AAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(view_labels.png) {2006-02-15 23:50:41}

image create photo ::icon_chooser::images::view_labels.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAyFJ
REFUSImtVT1Im1EUPVFp8Y8gcQikWvKpSEBFBCFQqIN/pXGJohWRigqSKP5kkMTBQXDQTl3b
CgoKxZotqdRBTFCjgxBcXBoQqlQUQxJMNJ8hPR2siWli81k8473n3HvefffxgD9wuVwEAIRC
IeIuvN4fiQEAgE73jQCQMTIywoEBEQCSacPDw4nB3t63DIeZzDw7u+TMTIoEAKyvr8cSGXcT
RUVFyezNTZL8lVhqb08kAPj9ZEqrJtMKm5qaeHJycpMMBoPs7yd3du5x9Tfa29tZXFwsjQyA
NTU1XFlZSS04PDxkaWkpq6oM9HrJnp7UNmLzaGlpxObmd1RXX0KhkMk0mneYnJy6FTErKyu9
tcvL+Hi0Wm2yQK/X02AwSD5kAsLhMK1Wq3SxSqUiANTX10sXqdVqaeTFxUVmZmZybW2N29vb
94v8fj8BYHr6Jz982Pl39fHxcQLg5OTNhTmdkvYnO0ZqabkR2Gw27u7uJm6+xWLh6OgnyuVx
qd1eCgCIRCI4PT3F1NTNrWcAQGVlJXS65wgErmS3ArP5NdxusrW1FS6XCzabLb3BV68M1Ov1
9Pv9XFhYSD6Tz+elx0M2NJBWq7RH8397lA4Wi4WCINBsNlOlUnF+fv5xGwUCAYqiSI/HQ6VS
SZPJ9DgNRFFkKBSix+OhIAicmJhIWzgjHeEulpaWkJubi9raWnR1dSEajcJoNP7fiNxuN/Pz
82kyDXF/nxwaIkWR7OsjJa59apyfn3N5+Qs3NuzMy8ukVttLny9ecH+fNBrJq6t47Pr6ml6v
l5FIJKlx0ogKCwtlSmU7Pn9uRDD4BBUVURQUyGIvpLwcODoCRDGusdvtUCgUcDgcmJ2dpcPh
SP0r6XQ6AuDqqgUfPz6VNTe/REmJJsHA4GA/7HYZBgbexGJutxuCIEAul8NsNsPpdEqe2L0I
BMixsfcUhGesq6vj1tYW5+bmCIBtbW0Mh8MEAFmaOklYXiY3NoCcHCA7+yui0TVUV79AZ2fn
w2p1d3cTADUaEzs6yL098uiIvLh42Abd2/Xg4IDHx8dQq9UoKyt78Elv8RvKs64deSgXFAAA
AABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(surfacesets_on_OneColor.png) {2016-12-23 13:09:50}

image create photo ::icon_chooser::images::surfacesets_on_OneColor.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAPRJ
REFUSImlVUkSxCAIpK38/8s9h4mGTdGEGwX0EtEISRERkpQmKkYCgLjb/hUAI4OIPElvA0CS
CGAD8IYwwKaJZIAZUCShVYzCwFVFQ24mPETwEWTrQqbOx+WbNWU23HRxpi9o3WnsbE0nJbo3
vTV0+pWab678hH2qhsIOVLKCB13MBstzSBd51ezBzO36zNCBppe5CpK4Mtol46ddMgV36icq
dKQfSYN34J07tk3gyUSeMz8lmu52WKCXbsJtWzW/cRMc7LwBJ27MA1Op8SRa0GzeONh9uj3R
yk2bFU5JtMD0L/lGfUVkCL6CZ0QiL96W0/gBBO3fUMLNWpgAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(view_swap_render.png) {2008-11-03 21:44:46}

image create photo ::icon_chooser::images::view_swap_render.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAxZJ
REFUSImtlE9oXFUUxn/3zpSnQq21UkVaJGKXgq+CvkmqaW3pwhFC0RZEFN4gNQGFmI3gVgWh
dBOQ7FpEK0oXBTdiKGFKMhgqdDClBWkZlU5clYodU/Lu5M7nYpzpvMwf09Df7r13zvnOPee7
DzowAJLEP2ABjDEmt5wjzfy8lyTN/TEnAxC9Fcl+nbkbsTBc1fKBVbWeswAzM14v3mqQ/ThL
1IhkrKWLaDpS99teVCoVjY25dLSfbbYmSdFUR6kvXlmVJIVxqGGNKDwTNnsHCONQwWgAb0Py
bpLWCeNQUXlwV30/pk54ccy1ekOS3DtOH4wmOn/edxWQe6MdnCKaijSifYrUbMmszzx82Omb
/YbXbr4ER4BngYcNGZPBs9adkJrCzgAOAM8DOwyL9qe+8QBcvnlZYRxucGf3wMCKqbHGcV2V
7yvU3qxpfNSpWOweZ0ewk//2rklUk24fcTp2zHUnFQp1ufHee1BVXQdX9aGq9Hfv+GgqUqRc
O6k910OHnL5z8OjZLbCz+S5Xy2G2WiwWzxqLZtGkFnHpktf8fIP385aRiyPwAjAE1mYQQmr0
3nS0LxIHgZeB54DtTWuUzILpcR9hcWHRcAW4DvwJrIk7hZXuPXTSeL0BV4HfIHkvoXy6PNhH
0DTguevn7r+X7gcClM875fNOpdIAi22UOHb6ef+qfNGnnVKR6sfr+v2s18TE//i5X+FSSXIH
nbTS240t/C9e9Ym6rp7yCuNQczfm+ooJ0Ozssty4k//RD67cy/5TkSJFCuNQxRvFtlDLRgYw
M0dh5oqofebxn3i4Nfi0ewt7yWkYezJD5mQWS4bgVMC1v67RddfWMz1d19KS+HyPYWuxgR21
ZD/KtgsHHwbwOPAIsMVgjMFiEaKBByApJP0FOjlxoq5fhz5laekHbM3C08BTwC6aItuBB9aL
NEgKqxsTaBGeCRV8FcAQ8Mx/IruBJ4AdwINgbPPvcqewQvl02dyTQFsoDhUogD00xXYDTwKP
QTKZvvabEgCYLE7qwpcXCHYFJNWEba9uo3i0uOl6m+ZfiH0p1YhojqYAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(view_higherentity_lines.png) {2006-02-16 13:23:38}

image create photo ::icon_chooser::images::view_higherentity_lines.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAApFJ
REFUSInVlUtIG1EUhs9YjU/ELiQxioSQGKmIO8EupNIHpQp2qdRKFy5C3avQVaUQcNOgpLsS
WjRFXZW2VEF0pRsRpFijoGBAa5JNH9RI1czXxRWnySRVxJb2wHBn7vznP485/x2RdMM/As8F
9RAZx4TQREQAZHBQhIEB0HUzTESEqirjBVYrALwUgFQPXK7MFL81CgvNXty7D0ffVayRkeNK
dB02g5A8gOrq402vF6KTCul0Z6AqKztHVqdmHY9DaSknlisQDmcI39ICk09g9hYkf8DUVcPJ
4UBEJJdEAnE4RGIx5RV5JaK5RXIsIpU2EdcV0TbDmpl9dRXsdoMxP/8PFPt3jUAAamshFEqp
JScjuKkJ8XhEwmGR4i3BasveAMqrjU7NXIfIBMdjnDbJsRiUXoY3FvXRppvhKKHAHR0Qjf6i
Cb8furoM5hcCq8PqvrMzDdzaCsvLBnjUotYvS3C7AuLxtFTW10Fp1QADdHfDzq6p2BzN49E0
EU3cbpHdIbW7sCDS0CBaZYV5hlKihcNgs8HQ0H8/R2cwtrehvR1WVi62XPx+KCmBvj7Y31er
CHh7YW3tfMGIRqG+HhobYW7OGIWjBKwNw9tyWJ+ER4+hqAj6+yGZPD0YT4fhksDdSjjcSyWe
bob3jfD1Iynm9YLLBWNjWf4b0SjU1EBdHczPK6fPH2D2DryuU9dcG+gHBqmuK305HGYlnBD7
fFBcDD09mCwyDqEi2HkHn6YgKDBzDRIxpUWnM1WPKcR7e+pMt9thcdFMPGqBrZA56MMHUO+E
bbMcs/d8YgKsNrh5A55VwLcpM3FvLxQUZM/4TIF0Hfr6IS8P2tpgYwOWltR4+nwXPPeHhxAI
QDD4b54fPwFgGrssjM1BmwAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(dimension_edit.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::dimension_edit.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAKxJ
REFUSInVVUESxCAISzr+/8vZQ8cZqqBS10NzQwlBAaUkVFwweBi0bgAgAJJ0u0li50ZS7kao
MlTsGCRFUt1GPtR+VulI/yfYuxihWINEQCKAu7LpQ5+/pdIutIevjRcSPKdpSrXlqppVXVKw
9n5rePkvTUSE8owwL+L5SqcVskhndERgdXhfC+zg+wLucHqI6jB6KgCnTe3vl4XHvVqHUbaz
4B73eA2+P8k/9a5iPcMFVgYAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(dimension_delete.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::dimension_delete.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAjVJ
REFUSImtlD9oE3EUx78XTMFFxAhZMnQRJFDnuggWKqQgIhJwCdTJoRBBIk6h6CZ1yOJQBIcI
VbBDHDIUUjg0EJzUYBwKIUgoWEglhlpKbvg4XP7cXS6XlPTBD5L33ve993vf7/0krxFLIEmG
JPHmI9re9aSUSiVGcMMC0eWAqCQWVwl08OBZcAU3OrI0yKZ3xEp68FumadrOgumqe06SaLUx
Ll80gnscnwyRzoZjAb0tDv6XSiUaWhgBEk8OB+1bf2BnEpkcrjGcl3Al1+qo80/G9WsG5W9j
RyUsi7yySLjO0K/R0YKMv0en4M2F/H0IlWogmtwW7P1y8FAsw/c9XxCVqn/B+7oEtbp7va02
KUXo6RpJCvWD7/VH+vJDNA9sEXUttL6ptzqccKeNPLTaeNXtJdcNunADWm1fQGgkOZbA6Hwy
tL6pidIgnsT4+cGQJOPVU4POEXQtjLmwS/tIoqj0iCwkKCpNWJZbFmzkfZU5iN986PnePRf0
Bb3b8YB6+5/V2N4d2bAd2KlAdHmiDscWLpgQWYLG/oT3r2DajTwyHJvf2Lfzi+XJa3OelCI0
tcC85vDGJHFeISq6SkZR3/g0w9ldc1tw5S5OjkhlYSUNXetMeLOLPn9t7zie9CdxpuKVKsQS
UKvbr+3iaqCupy/cPID5274EksrCrbWpG7nI4fjEnjKTC346uxbcewJ3Hjv58CV64GTtxakJ
pGvZmFR2vJL4/HVmAgc3f/TybEUQZP8BVCPG/Cv1MtEAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(line_divide_nparts.png) {2006-02-14 22:32:48}

image create photo ::icon_chooser::images::line_divide_nparts.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAJlJ
REFUSInNVFsOwCAIQ7P7X5l9LJmL8mjVLeufUioUoogD9QIGVT1+f6mwKkRUjSSJDj4DboCI
3YHrNa1uvfBHz0DcWxguRSYAB4a72gJl4nWslGXylSHOpGtGmMV2wX1gh/VsA8qtcTj+xBqH
BpMUc/s2GQsgbmJRql+AQobENzxHBXhxdgb0mh7zBVndlOGw/Ct2oDtcxgmENUbWxddOJAAA
AABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(send_to_back.png) {2004-11-30 20:20:29}

image create photo ::icon_chooser::images::send_to_back.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAnJJ
REFUSIntlDFoE1EYx/9Xr3DgpdxwlRsidMjQwaGDpBkKXTo6JtLBQukgHXTL4BDQIaNIhiKl
QpGinSqU2iGDLuFwOjUUKQGv5ZAbTg36KKk8mtO/Q5qY3MWkgh0Ef9Pdu/u///d97/seEKWw
tkt1LEkAgJ7OkmTkF5IwKy5JYgQAUK31fiUJdURl5xkAEtM5QtOgrZSRmM5ROtsx9+EobY+e
RUU5jQTAwddjXn7+uvPHCAAkZhc4b9cg3QCJ2YVfcpLY+lCnuVTsk+pfo10Jkihu2jQ3bapj
SXavd4Jtc/vRSz57e8ASACkkjO0yxu8+5ejF8Xig7V2uXr9F51OD5oNdmjfvxxzUqLDmBZgv
rAKhhPS8c0j+fOhOjGSnWaMJ9xXo6SzdL43OsQ8U6Oksw5OQ+T2fM/t+jygm0NNZ+kcN3qi4
1NdtJnecHlEbBWgNGHQDCCWmFpdRBWClDAT3SkAoAQCNyhOlbw7+t5CFVy7n9l2ql1KxkHp6
CQCs07OvBrJvRWOtMTV/B4FhtF5kf9F//oBof0ZbMv+wNbTRq2FgXw8z0NNZFtZ2GZ6EfHFY
Z+qNS3PDHmh0JoPOxj9C+kcN5vd8Jncc6us2zYrLa4d1JjdsmkvFvpddGyVqkJjO0TIMTE5M
wgslAsOAMZNBICQ0KSEBaBMWUPMA3wOEAISArAeQThnN48+Koii/N4gyemGUq44HqWt4/M5D
TUhkMhaqK2WIrRKaH98rUU23QWzUokzMLWLxioWMoWE5ZUETYpikh9hoRgmEQFVIrDoeaoGA
1LTWnMuzGZ2pRNANQNW6wtJa90IjQPN7c2CJ/n1+Al6vBVnj2YnXAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(PantallaResize.png) {2010-02-12 13:01:02}

image create photo ::icon_chooser::images::PantallaResize.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAL9J
REFUSIntlcsNwyAQRGdRSonE9uIG0keuPnJNH2nAvRjJvZCDhUVIWH7Bp7wLFoyZndVKACkI
ANhYN12TGi+7P62TRLuq7MIcxEYwiw/ZWHeYh2UuG6B66qhErDqA4g2pY00hzkw9DDbW+Zbm
1o8fmtxyGmJj3Trrr9MW87hpUqG4BFVb/3iHP4V0jZ9wZ9thJcdTkJqnHrN11nTxH5Kw5dlY
tn1VwG9bEqOAfIJug5EJ3ir3Rj5Rr/HIzpzHC57HZsex/w5yAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(volumesets_OneColor.png) {2017-01-20 13:33:29}

image create photo ::icon_chooser::images::volumesets_OneColor.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAO9J
REFUSIm1VEcOwzAMI4v8/8vqSYZsU6MBylNiDVLDppnB8UHA+iFp7N0AYHMjaWbG5UXSAMAP
AeCJ3v5tZtxSXYTRezPE3O5Up1IggCtEqopFthTurNRfATGzCmqrKxkqCaf9Z4Z3kiJ1J28L
qPovA86gdJtO47nB0bYCYubxakRkDZCDqxbwYjh36rL74LKs5/kTJXVD2ximKAO62zdBOYZ4
8SZ7qfB0Dk6iVnJSnWyRak011YosbdEJdS/8vLpw5aJmRNmDqDCuQAmJpJk4+c53yb0tk60q
K+gSTEi2u+9BnapMSPm4/AuvhzzFF8bX2jsVsIexAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(collapse_edges.png) {2006-01-31 15:41:36}

image create photo ::icon_chooser::images::collapse_edges.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAMhJ
REFUSInlVUESwjAIXJI+WF9SH6yzHtrUQCC21Xqoe8osBAgQACxIEgBSI/FAcIbmx6wIacjr
QwCYW1/AmGmdd14yZi7xm1sHoU2Vhgq2KAeXCJTUBkoiIpZLVlluw+rwzgGnSbrwe3wSlcOr
cLXljpewcKZoLLoptKg5UYeoon5r9LLymw/yZ5gn0DHJrcbbCifEpoFsjO+Yo4vDdkJ3jPBy
j+ThCwaPjFDPw3eqPm3Wik3PZyvHyXvdRXvr0a7orZhW+onxBJFa4YDC/4qsAAAAAElFTkSu
QmCC
}

set ::icon_chooser::ImportedDates() {2017-01-26 13:49:33}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(pan.png) {2006-03-01 12:27:08}

image create photo ::icon_chooser::images::pan.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAc1J
REFUSImtlD1O60AUhb87YwmNFIkyDQXta9JkF2nYAhUroB0pKE12gTfwitfQsIY0XoGLNPSW
LBQm9xVmLDueBBs4lXV97v+5I7m1RJj4sXROM4AiBKWuhdxaGqqqiZTcZmIAdnUtANKNBrAK
QV+slaxrLELQAlgAAqq5zYQTDEJFZAmblt73PZbOaayqbWsVQmtsf8QSe8knV3WK+/ChFx2W
zinAvKropjdJ9ifxuap4m8169sE8ViHoFsBaWQC7uk47LJ3T56pqiedgYkN/Hx9ZnIwqhUlT
akU0hgxQeq8pEbTYG6MPhwO5zaT0Xm83G0mWdLIDgMt72NW1/Ht/H+zgrMPeGL05HgEoQuj9
G5TUISdHbKaQew5jyD2Hu6urL8kAoqoqAqljSWYQgfX6aQy3yTBWS2NwHz5UFeSz1txmclFK
Y9GVnEij0/X6iRKvkzv4o6qv1vJwOPTsW0ie36QEUfNdpPTfxajzjIg3Mq+q9ul5m80G8v9W
gr0xugpBb45Hto1plOy+HFFc4Iu1FCGMeua6ONtBquKpwZMJ5tfXunTux4EjJLe2fWPjORch
6E+CdmFicGgusfT+14JDZ8ml9wpwu9n8WnCA/xep5/o+khjtAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates() {2017-04-08 15:17:47}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(bring_to_front.png) {2004-11-30 18:57:31}

image create photo ::icon_chooser::images::bring_to_front.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAlhJ
REFUSInVVDFo20AUfbIl0CCDBhc0ePDgIaMHx/FW6FgyZHA2FzqUNGOmkiE4U0uHFErJ0MEO
oVNoM3gIJEOGYkyaggcNLnhQgiimCGNaQVU4iJvXwZFtWVIaKIH0TXfcf///e+/uA7PQimWS
HG3khDzZjEES2ov90QlJLJ4PJmEkoRXLzDyokCQCybabHWpvDqjttCIS3xgkkXv4hIFm/QOr
7zG9sjU+kAHg+eEXvhQCKJVw7+krApDGjOHlMMC4BfhKaMUyaycWN06scclplRKzRF1TYQN4
3e4CmhpKHCAoSYX28S7soQA0HfnVJf+28S2RhJyQQ+3cohr/CglAoL/UwjLnMga6PQcA8PPz
BwkAJEmKJpz9+MVd04au61gvZHHx+yJAkGdLZjUVtgw4sohsKUQAgMZAwIi5w9i4+c0665/O
uGc5QM+BCxXq/inmN+tUkkrQBN+YQrVG67tHbafFdNNi5l2Lkb912s1CtcaD80Eg+FqCT5oO
vsPP4u5jVlx/JGw3O2x/HdDqe7T6HjvfPGaehYWPMiDy5c1itTSHo56LNdOG2bVh5LMQrnuj
pkMz6TqoKgBPIGvoV5u/I3CD1MIygdHgcz0xSqKqOO25OLJd2AMBGAZMR0DkslALS0jdf0Rd
xiReCMCf8NMLYDIGlKTC/MZb7K09hmnbaHhAw3IAbySLYRhwDR15Q4VtOnCPP8KtrYdGRQiz
ZskJmYVqjcPLIa2+x0qzw0qzw8Wrr5he2Yo0OhZRL2K6kNX3mHvfjk0cVSBSojj4Y8mXIg6x
Ev2X+AOcTOsR2+iergAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(user_open.png) {2006-02-13 20:59:36}

image create photo ::icon_chooser::images::user_open.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAA/xJ
REFUSIm9lntQVGUYxn9nWRaXdSFRGbKiLcCwTDQbcMqpsHHU0WpkyBlycjJnwKZs0j9qZKYr
XWS6khdQ/2gaE4VyKcsJNc1acAqxApUZo7wMiBTCsrvssnv2nPP2B8EAG7ep6fnrm/O833me
5zvveeeD4ZDjzwqAAvDlnBny6yOLh5ecr5T+tamu9DsAml3b+x7WH3hGDr6EDN5hAph92cWs
lmSqHLEDpDK4KjU1TQ7u3kA4qEbYGorad++Ump2ZEkFIbaFIsGsoUV9UJJfW3SC1D04fQigA
X2fdIjGOuVRNjacgJ4N4u63P7rIfLyqLKqoUc7gDW+wk9lZfGMPVaDj14TyRph2iuTbL3o2m
SOd/w9y/UI9aaLq1m+Zjl0g4logYbaKYTMrwDQrAiSXpkpKuoQY0fJ0q3e1BVp1ys6rgaQrW
55OSkkqt61vC/k4i3tCP9varsqXwSfIfXwaAoelMjrGMnqvmqFOc5WXicDhGzPTfYyDDJxtN
kjwjjdi489xdMHK2AUjDOyJnSqT8uVgJh3xjfwfcrRAVQ979S1Fi7KMrdDc2SuUchzR/ukkO
pVqlsXjLiAoKQPC9uaKHwnQ0e3B36bjbQ9z+1gGSsrMjlEwAna1+Oq4E6VUVgn6dvJ88VBW/
TmbWAmlpbZHNLxbKIedHkpvzsJgAlOVlXLx8Pb//EsLtT2LDE+vI3r6VtBtjcTdVsP4hB44p
XlIdN418fM6vvpCspA48vgAooKs6ew//1mfpn3Clta3Ps0nBbEBS+mLqfz438sk5HA7p6bwo
O0ve/B/7CCIzqwGPnKx4is4Lh4m2WDHHmJlkT2b+yl3EJ6WP3WJjCXz+8kxZ8cAKoq5LBEXh
bN0xrnYeIRRSWPJCNxZr/IREIorldLEgKlhsYGgQ8IDvD874LMzJ3fHvE4jaI0ZFHqZpCYAJ
DB9drTbsj23DYh/jHxlNoPfaNTmXfx/zcydzYmcbjRlrkT/bmHnuADfPstHj1fB2qNjvWc2C
rdvGLTRQWLt6qcye1YIWFnTVQA3q9PrDqD6NgFcj5NcI+DTCfh3rXctZVL5/iEhZWZnY4+y4
XCcpK92uRAjUrVkl7oYj3HZvAmrQQFd1wr2C2qsTCoTp7dH5OJRI/u5dLMxaiOuHWlpa2/B2
XaPmxHGmxBpk3pFE4vSpzM9Iw2qJ5u09p8cx2QZhZU6O7N+3j/LPKmn4Zg+vPb+Gji4valgH
QEQQEQxNJ85mpdR5dtDAGweqnE6l+nB130XBZMbnDxJStQFedAMxBMWAOJuVaLN5YgIAr75S
xAcl79MUbSbKHEVCvA2vJ4AybR6ebg+bCt8Age9rTk6444Yg79GVkr929Zhz5y+Tp8ar2ci/
5wAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(send_to_back_opposite.png) {2004-11-30 20:24:54}

image create photo ::icon_chooser::images::send_to_back_opposite.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAhZJ
REFUSInllD9MU1EUh78SQ47b7fa62bHb04mSOODoxnOrgyEkLjjZSVkZwcEAkzZxsEyy0Q3H
1qnPiTK9Z+Jw2e7d3g3LcSgg6ZNSMSQav+29e/787snvXJhkvXOgYiKtADybn9dSBADRIFOA
ypUhAhcHcwBmIVGz0FJ2D5E4UZAriv8ZR67Qem842bylTwcj3LFFllpj7VC+4ofT08ptqLrE
xn5fTbevmKg0gLnLH2udQ/2UZrodAO+I9ntErz6qiLl6cs3lNR26Qs2bA41WNkuBdyZ/pLnl
ycsdIOC/5Te60t/A2FMzYuJEM1dotLJxfZKJW6qq2k6tNkeZmmlJEidqC6etQabS6WvUG2pz
ZEtJlfNgqVYJRcH95y8YAbWGwb7egqIAwH/p/tqRtlBdH2T6aJQppl6SdM0K/+Tc86U2jeW2
OlMdy9jbJoST216O/xeJE23vni2tlJ+GG2PiRNc7B6qqephbraeZRt2+mpVNFSn75rcUnxfO
nNP2yGq91x/bfpDp49xp1O2Pl2xKo5KnTJxotValETXICThjuPtwEecdhLOgezXkOCd8zxHv
Cd6Btfivn2HCp7OYVt+llmDgfWrJveNBs87R2x4ne1vg86k15qYdAtSXVlmNIxZFWGvUEB+A
QABkBnWlV3sS6xypD+yklqMTTxBhaEGCx8/QYKYRiYkALlQHuYt4CCGftcY/zA9GQhliq7/h
AAAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(bring_to_front.png) {2004-11-30 20:21:40}

image create photo ::icon_chooser::images::bring_to_front.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAlhJ
REFUSInVVDFo20AUfbIl0CCDBhc0ePDgIaMHx/FW6FgyZHA2FzqUNGOmkiE4U0uHFErJ0MEO
oVNoM3gIJEOGYkyaggcNLnhQgiimCGNaQVU4iJvXwZFtWVIaKIH0TXfcf///e+/uA7PQimWS
HG3khDzZjEES2ov90QlJLJ4PJmEkoRXLzDyokCQCybabHWpvDqjttCIS3xgkkXv4hIFm/QOr
7zG9sjU+kAHg+eEXvhQCKJVw7+krApDGjOHlMMC4BfhKaMUyaycWN06scclplRKzRF1TYQN4
3e4CmhpKHCAoSYX28S7soQA0HfnVJf+28S2RhJyQQ+3cohr/CglAoL/UwjLnMga6PQcA8PPz
BwkAJEmKJpz9+MVd04au61gvZHHx+yJAkGdLZjUVtgw4sohsKUQAgMZAwIi5w9i4+c0665/O
uGc5QM+BCxXq/inmN+tUkkrQBN+YQrVG67tHbafFdNNi5l2Lkb912s1CtcaD80Eg+FqCT5oO
vsPP4u5jVlx/JGw3O2x/HdDqe7T6HjvfPGaehYWPMiDy5c1itTSHo56LNdOG2bVh5LMQrnuj
pkMz6TqoKgBPIGvoV5u/I3CD1MIygdHgcz0xSqKqOO25OLJd2AMBGAZMR0DkslALS0jdf0Rd
xiReCMCf8NMLYDIGlKTC/MZb7K09hmnbaHhAw3IAbySLYRhwDR15Q4VtOnCPP8KtrYdGRQiz
ZskJmYVqjcPLIa2+x0qzw0qzw8Wrr5he2Yo0OhZRL2K6kNX3mHvfjk0cVSBSojj4Y8mXIg6x
Ev2X+AOcTOsR2+iergAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(send_to_back_opposite.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::send_to_back_opposite.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAhZJ
REFUSInllD9MU1EUh78SQ47b7fa62bHb04mSOODoxnOrgyEkLjjZSVkZwcEAkzZxsEyy0Q3H
1qnPiTK9Z+Jw2e7d3g3LcSgg6ZNSMSQav+29e/787snvXJhkvXOgYiKtADybn9dSBADRIFOA
ypUhAhcHcwBmIVGz0FJ2D5E4UZAriv8ZR67Qem842bylTwcj3LFFllpj7VC+4ofT08ptqLrE
xn5fTbevmKg0gLnLH2udQ/2UZrodAO+I9ntErz6qiLl6cs3lNR26Qs2bA41WNkuBdyZ/pLnl
ycsdIOC/5Te60t/A2FMzYuJEM1dotLJxfZKJW6qq2k6tNkeZmmlJEidqC6etQabS6WvUG2pz
ZEtJlfNgqVYJRcH95y8YAbWGwb7egqIAwH/p/tqRtlBdH2T6aJQppl6SdM0K/+Tc86U2jeW2
OlMdy9jbJoST216O/xeJE23vni2tlJ+GG2PiRNc7B6qqephbraeZRt2+mpVNFSn75rcUnxfO
nNP2yGq91x/bfpDp49xp1O2Pl2xKo5KnTJxotValETXICThjuPtwEecdhLOgezXkOCd8zxHv
Cd6Btfivn2HCp7OYVt+llmDgfWrJveNBs87R2x4ne1vg86k15qYdAtSXVlmNIxZFWGvUEB+A
QABkBnWlV3sS6xypD+yklqMTTxBhaEGCx8/QYKYRiYkALlQHuYt4CCGftcY/zA9GQhliq7/h
AAAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(viewprev.png) {2006-03-01 12:30:27}

image create photo ::icon_chooser::images::viewprev.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAqpJ
REFUSImllMtrE1EUxn93kj5sMBYj2IelFM2iIIqLiih0IeLKtoG68T+wuHTloll00W7apc0f
IG5KQkaEQqW4EDR9LKSlRUrEkj5jSctY6SNtk+MinWGSTNIEP7gMd+53zplzv++MEhFMaNjg
BphdXpf6Wg1VnrYQ35Tjk2yOppQSEVF2lgYQHu6lEFbid9PLsm/sAuBt9OVXzEullHI8KR1h
30RG+iQy0idWxPvJGWm/0cynr995+ugeiY3tXMRKIsnSmsHQQEAtrRmsJJLFNQqbLdqXupFS
0ABi0VHHQ/vHWxARRIRYdFREhImpWTEMQwbHoyIiDI5HZSe1JxNTsyIi+T18/BaXxYTBYfqM
tmsNbK2v0tLWwXrqkIY6N3faG3MSm+h56FfzC7q4gJfPAgruAhAM6ZI+gp4XXaqkEmWbNlGo
r5PeFVWwX7tWKpMdIqJMTvU9OIpzjhl9rOi9BtDa2e2Y7UHgdbFdTKXDw70yOB4VwzAsVU3l
d1J7lvJWD8GQLv1Puvh7dEaNS+NnYhOAW+2tnGayXL7kJjI9n690bCVls4EPgOnFpGUVKLil
t+Ev4qq/gqfeTTweB8Dv93NwfEbm+A+vnncXWyMY0gVgaCCgnPYVK+30vpKhswrY7WImrHRq
y/HzClSTtFyxon+MvYCJ/y1kwhrP8HAvIqLMVUlwqZHIg+lq+7OSZXJj0VEpF1fkoshIn/S/
+WB1YNrOhGk/O8pxHG06o4/J5PZNOjuauX+7g2w2x9E0xdzSKj9WtxkaCKhgSJeLOI4FgiFd
Ao+72D884eA4w2kmC0CNS8NT78LbUIv+eZ5KOO6i7OeYW/6Fx3uV3f00B+kMAJ46Fz5vHQer
GxVzyk5yMKRLy3UfrU1NAGwmk2z93s2744s4Vf+yq8U/hvITFDKnFYUAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(sendto.png) {2012-12-14 22:39:32}

image create photo ::icon_chooser::images::sendto.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAABLhJ
REFUSImVlH1I1VcYxz+/X1dn2nYNsld2u+YtFyg6y0qmzV6UgouzZRaNQlNzOXRjbNe5lCtt
WXNiDUuKyI3IkpEoyZrF0iiVujasxO6aL6vVpbFeZqnkW/fZHxev3q6u9YUDh3Oe83zP93ue
8ygiwgg0AENpaQKgAizRaOj5ZhbK2DDln6Qkifb05LvNm2FzQ4MMDn4tIgIigoiQEBYmI3MR
cZwfGtotd+928+aebobsdryPHlWciVevWSMz8vIAyNFqXRnHQgXo6uqSdwIDJTQ9fTQqLCxM
BlNTZTA1VTo6OiTRYBARQbNsmSfs17NucRlVe6Hx4UNHqrOX+kiyRLM2cwmU6lGmaBxGdNy4
oQC8Z7NJk2UOk8MsDukT3WoiaEYmwYGBciUyEg9VBSBcVak0mejt7SUkJERxmjvCEB8fLjbr
PZqWG53ZwlWVZrsdgEcDA+zo6RllaLn+kIyPFsPHegDq67vwyvyJJ40paLVe3Gq6Q/ZrGaMH
FkVEcGHqKm5Y3gKg7Yc2ht/wJcsai6KqPPr9JvnB/yE6IipKiIoibuZMcrKy3DX8X6ivFM0Y
WwGMq1fLRrudRl9fLvv5AXDt8GFlbIwKUFNTI7sjIqTK359NAQFc9vOj0mQa/0rt7e3SvGsX
pqAglw3TzZsOhoQEiV20yClUERE2bDBKcPAAgZfh/Tl6wPFogPPhttTXU9HRoWiio6Pl1Kko
tFovAJI3/chnynyap04F4GBbG9dmz6aio0NxMgDoF+ql1ZKMl5fDh5XLj/DukrmYi2Id+/rv
sNn+VpQzZ87I9xoNHt7eDPf3Y/lkO7d+TUZRHOYUFV3CYoGqqtpRhs7OTjFu20ZYQQEAA0+e
YDV/jiFhK1MiIwE4nZhIj82muLx0S0uLpJvNzM/Odq71P37Mnf37uXr+vINybFcYGbW1tRKT
nCy6oCB5ce+Va+lV8dLaq6mpkeLiYrdbhKanu/afCTCuArPZLKfLy9kTHMyqadOc6+fu3WOn
TgdApcnE+sJCt1p+Ec7PYDQaZbLNxrcLF5Lr40PuihUugeGqCjodlSaTs4wB1gUEiH3GDL7c
t4+lS5e6kbkoKCsrk715eWQaDKQaDGhUdwdHSn5EwUjpA3Q9fUrB9evc8fBgZ3ExMTExigpQ
WloqQUHzRKu9QNvtDD78JZanBbPIfv1PVl48x8X7951Jmu12mu121hcWuhD3DQ9T0dmJ9dkz
MtLSiImJcf06J06ckOLir4iL05GVtQwfH0+XBK2tf5FrOotfvwf5uhDm+PhQ3t7OwT+sJH4Q
QnxSMPn5dQwM+FJUVIa/v78rQV1dnWTm5DA3JQXtggW0V1fz6HwNmdtD2JH6NhqNq112u6Cq
o5YPDAxz7FgLZWWtxMaupaDggCvBppQUudrWhv+WLUwPDXVJNtTbi/XkSXpbGsn9IpKtCQYe
POijpKSJ6p9vozfG4Rm6nNYjR5ip0VB+4IC7grGwWq3yaU4Otr4+dElJaPV6l315/hxl0iS3
c92dnVw7dIigefM4e/z4xAQvoqGhQXaazdj9/Ji2cSPe06cD8HxwkN8qKui+coX8/HyS4uPd
/8R4veVlo7q6WkpKStz6znjjX6P4cngykriEAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(surface.png) {2017-02-24 15:48:38}

image create photo ::icon_chooser::images::surface.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAilJ
REFUSIm1VT1r21AUvenYJT+i/yIQOx/GUDQYQwhkMSiYIA/Z+gsCLSmebAjeMheyZavIIHC0
JGDIkvGRwYIWYkgqAtJ9fpwOD4lYepIlaM8kvfvuud/3EWmAsh8lwCtAROR5ChtERIsG48un
DS21bQYREULg8JBXCUcjmR6ER2He2vGxrOCC+qFARPSn+45//nG+oroinM0U4s+xmXs203QJ
7QqytClsW0JcC4RHIZwmw/MM2kQ6JnbYKOQe57NERNRqMRYNBn4jJ0QIOE0zIRHpNGcVxbUA
lRW41WLgFWCHkRSxCipdLL304f3P7dYcQggKdiN0OmvcuMtc4B7jtBnj5sZQi+lUQbl5AULg
8dJwzgfF5tWDwmCQkU8mChd7EeRXcxfLEwnfN7jWbjOev0ngV77allXgxXAoi3q6OGv9PgML
bQVvwP1OtL6ALw2G+qng+3n3jOhuBnDdoHIPJaitkCjBshiWxeb01oVtM+53IqjMSEMA8kTi
6UoXv3Dky4h9H+B9Bt7Kk6keFORA4vFSGzP2dOIYEcF1A7CjK1XHKYSAHEicNuPyiLqbAUbb
MV4ajOXZMu2jquAe4+lKrR9rIqLxWKLfZzx/110uz83jkYtGAMFuhMKxMWE4lEhmUJ6vf4KE
SJdpvQZot3VERRuCiEh5CrdbBa9OFYzHEhd7EZZny5REuQq8z7jr1NvuhZhMFEbbMfiAwQ5j
OlX/htiA+nn+H/gLeVfXYRFuK6gAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates() {2017-04-08 16:15:09}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(intersect_lines_nodivide.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::intersect_lines_nodivide.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAKdJ
REFUSInNVUkSgCAMK4z//3I8iE6FtBRFh5ygW7qhIi0AIrTMspaIdclUalnH6Ic9BNzrVu+D
sC8BFFpOjcMimtdZTGyI/1f7Dbp15Pu1tgfUHFyieMd6Q8ue8lc0aZJO+HsVad2iSI5Ol+TZ
udieOhLQhLxFSsWHZH/t/JTEhj47Nbq91W8pvZiFFR78PD14iMT/H8QDmTpKUPd0YPc1UZo9
m4WwAxaKUMom9t60AAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(user_save.png) {2006-02-13 21:49:11}

image create photo ::icon_chooser::images::user_save.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAABABJ
REFUSIm9VWtMW2UYfnoKhcKABR1RNLUQLmUiOmd0I5vOTTN0TGXDIG4a3ZQFXYzMW2aikuyH
LvG2S8aGxsyJCxvBXcStwrgM7EScEphDZw1pAxQYlxUO7Tk9Pe3jD0cVKsVN4/PrPd/7vt/z
PN/53nzAdLDhBQY+vsxKpCRLnFZy4XBgJayt7DR6OjppD9f/sXK2ejOPv4kpPQIAZNpbkNFj
wBFjVCCpmU6/v8bMsXExSNZUWN67hd/su5NBCVpeJ+XRqYmz27bRtvEGWlbMIwBkL1vxZ8HJ
u5LYUJBHr+ql+7LVKaqs/RcZHyGgovrYLKpC4fudC8iuPVRbtvLzEiFY+WWETQZKnQ5dyU5Y
622Ir08A/Q5qBCHowAQAaFppoiF5BJHHymG0n4Z+joLjaTF49LFCSh4PAaB4yyvB5v6K9c8U
8bOP9kJVFew3N6Ou/hRSrpkb2teTzxbxw6MnaXf0c8Nzm2f09d8i4KGiRKAhMRVRsRdwx6aZ
vQXAjnfJczt48MUoej3i7P8Bl3oBbQQK78mBJiImNIOzs5OHs4y0Vm3hVyl6dm5/J/SJyO/f
RtfbN9O24Ua2P3I9GxbFc6CxMdB0/4OrA7EAACO9Lgz1yZAUDWSXD4qkornoKZR/eoB9Q8Pc
+OrWqQyO+no2rl7CmswEnsjO4hMP5XHIeYkexUO3LJH0cf2m50PLNLf9SI9XodUxwLUvvcGS
XR8TALQzNejDdaXN9mGkzotF3vKl+O43O4Tx0dKQLM0dXezstrHmVOP/dI/wN+OguMd45lAx
Rrq/RrhOj7CIMETGGLAwrxxx15lmv2KzERx9K425y3KhnZsAaDT4qa0e/SO18Hg0WPmaEzp9
3BWRBBXzh+0EFUAXDfhVwD0GiIM4J+qQlb/n3zugMkH/oUII18YDEAC/iNHeaMQ8vhu6mFlm
JBSBNDzM80V3Y2H+HDTtc6Dz1qfBiw6kna/GTRnRmBhXMT6kICZ7HRbt2v2PiQKFlnU5zMzo
geolfIofiuyD5PJCEVW4x1V4XCrcogqvywf97auw/GClBgCMxiTuKNuLpdmLYe3pQ0HuA7DZ
bIF9hckgXBuL1qo+iEMynIMyxCEZ0ogKacwHxaVCnvBhtFdCsV2FxQs4xQkOiBP85ecu5Nx3
L6IidTAlG/BJZRVMGfOvbg5S0+fzTHsHZZ9KySPRLUuUZIk+n5eq6qW9f5AffGHmCUsbk4zJ
oZ+FmWA0GllRfQTJRgNEcQK17V2wdncjKyUJcfDC3PErwrQCcpcsRunLJVdOMIlaSyubrT1Y
tSAdByqrsHbNGpgMifD6/NhZVg5LnRltrd9qrpoAAB7OL6A2IgJZpnQ0Nbegqc4ctN/vLNP/
g8sbH4YAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(view_normals_surface.png) {2006-02-16 13:04:57}

image create photo ::icon_chooser::images::view_normals_surface.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAA8RJ
REFUSInNlX1I3HUcx9+36Tmn+ZjIvAcXPRiUFISsSZPlHsIJ29wsrVjEYIOEjU3aCgqOUtqj
2/yjMfKYLh+6Q2vZDv0jLiEirGCbQhiUdpHO7TzvnPmwU2+v/lBuu7k2hRF7w++fz+/9/nw+
vD+f7/cr3Ymid1YxPm5GkgQW5jG0fv16wAyYWOLf/48gJINhwBDBGhoKAlYkKUqScr7PZcWa
DD2qJyTdkbejo4v50fthbCyBqiprpAqsDA2ZABONjfWo/7d+nnwvi+z6VdjtKZSXV84qfH8P
8awzl+npp6itrZ9f/MSJZxbXUXX1B3cXtLW1hX/Y7WcJhdIAM1NTZpzOFC5fjqWr69dIsbk6
kyJeI/tCHr//sYJAwAxYCIXW0dT01SzZ3e5m3cQmNlPMForJbsoFTIAFsHDuXD0lJSV3b+vw
oby5UaZz/PiaRdr/v6ClpWHhbdlsuYCZS5dicbl+mi+cnp4OB9vbo+ZcMuNwpABW9u59K1KU
k5NDYNBPZ4dhzikLNlsCkAGYOX8+loGBwVuinW/vpGBmC5t5kwZHKpWVCYCFmRkTra0xYeIS
SXqxbDWD1T4tXxqnJUIZqRJIkkFu97g2bBi51Ury86kUUcJmitnKG7Q2JwPpgJnTp6OQpPz8
/Mj+P/7wI4rZwYXmZMCE378Wtzvp3pb2dffT2hoDWAgG1+LxXFnYDKqqCnG5fmBwcPBh3KUH
Ba93gl27XqG9PYWLF40cOVJBT88CPboXzpw5y4EDLzA6mhhezdkvE5criYoKcerUu/culJeX
x9GjR8Mk31U/NtsntLXF3ZZwdpW9XhN79sTjdKYC1rm4lUAgge3b4+ns7JpfzO12c+zYMX75
7mcKjhTy0p8FrC5fyef1KeHEw8MZ2GwJjIyYgHR8vkJCoUzARHPzcpqavgwnjng4pkamONte
q4P2g8pz5Wvpsiih2cchWvG6tP9HPXbRo08/S1VWllHSzTlljByONBkM+1RaWhqR0yBJ3ite
Xq0rVez7jyjWsEwzoRmFxmYUnRgtg4zq+8KjQ7E92rQ1XhIKBiWXa1LB4E2lpU1qdNSplhan
pJtyOBwRBaIkqburW/HRcbra7ZX1uUwZlkrGxOXqc/ylw8YeFbwepxs3VurataeVnv6tYmKW
yWgc1+RkjTZu3GGQtqi3t5eampr77YlU11jHin2P883X8UDmbUNN4/r1l5mYSGb37v+49hcD
n28Su72BkyejCARSgAy2bRNJSen4/cMP9pj7/eP09s4epLKyMjwez8Nzj/wLwitimp3YZ9gA
AAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates() {2017-03-10 12:11:55}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(zoomin.png) {2006-03-01 12:35:26}

image create photo ::icon_chooser::images::zoomin.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAA3ZJ
REFUSIm1lF9MW2UYxn+nLU1bKH82GWrHSMaQIcTIdInGBOKS3U0Ws13swhuCFxo1xlBNTDRZ
YgLRbN643pgs3s0b9odqNLJEAjcQIRRHmwqjJcA6bU9pKWiB/nu9YKfj0BbE6JN8N+f7nuf5
nve836uICBpMAMO/LslGOocBIPbDezRY4ig7jxkA3v6xWwBM43euSnno9e0tEUFE6L02KOpK
XHqvDYqIbFMAnj52nDe+vkBZ9dHHjLHbV0RjakvnuBN5qRtDkxKNrYpuw+lyS3nVIb6Y/4TJ
uYeiY9z1RglOLzLsCQJse6j3R+XzoVUArrzTpQAw0Ncl3pmpf36rUjDt/hB4oEoylSX3SMds
MnD9+/G8tY7gdLnFXtvAxIMPOdd5impbJb88XObZZ94v7vDy2nWW6j7iTNOXhANZVozwXI0F
79wsnD35OPX4navS2HqagG8C+4kztLa1KyVDaOmLVaTYOnCVDPsf2YfgdLklGluVaGxVbgxN
itPl1l1BdyWnyy09514ilcltqylgMxtpPFqbL4LO4dLZF1lQkyyoSaLmWwTVJEsrG+x00RGm
F9f4aUYlmPuK5VgYDt9idDZeOkM8FsVmNjLrvcRy7A/ujnRgtxp1hHyGm/3nZayyh7aWZkLx
TTJZsFuNKFsJPrj4Sj6DSTvccbEXh28C+5F6vP5lANaB7s56fRkH+rokMjcikbkRGejr2vdv
K96ZKVmf/5mQf5QLHw+W7iENB+mjf9VLB0XBCy2F3S2lIT80SmDPBE6XW54/8RSvtjeSzT2K
rREBRVEwGhSGPQGm538valbSwOlyS0tDHe3Nx9hM58hkhZwIt9XPsFZZKLdU0ilvYTIqWMoM
eGaX8C+GC0z2LFFTg4NQbItkKsu3U+9iq7Fx3PEkL1TVM5+I4k71Y6+o4dT6mzQ1OPAvhgs0
9jQIqklC8Q1W/kyzFvuLVCpNpKyMpXILkZUE8dQmyqaVe/EEjhprUY2SJfJ5PfLNyDJtLc2E
EylWk2m2MjnuL/SyZTZSUVFF8+FPqbaVUVdlxuufpbuzvmBeFhjc7D8vjpYOAEL+UcYqewCo
e+IQtbVHyGS3Z4/JaEBVI4SjMYCi4joDTbix9XR+M+CbyL9An9dTNOqeE14z2C2uCZ987fK+
Avshn8Dn9chv310G+E+ECwz+L/wN9w8tS3KDPcgAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(note.png) {2006-01-20 17:17:27}

image create photo ::icon_chooser::images::note.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAQlJ
REFUSInllV2OwyAMhD+iPRg3q3sz92TTB0hiCGmSKg+72pGQ+DFj4zEAAQkAM3nOZcbdxRhm
ApgAZvsU1uO+NJUpCQkk3H2H9woU/aw+zEDqgijOVU+0xk5KyV+vG6IZIXVj1cCCRWpspyPj
PnM/O2wAwzTHDaka7Yn8W9GnFdqC29hMzVIoxFCQ2ttQFp7P0mCT5tadetVgLrCc82I3e9jm
fsAOUbgTKkcPiVXpxTiG0m8Ywt11z7X9XxiLuUJmhtWB1db0bZ4Zc+2rZiakhXAmXfqPx9Kv
yl54CAr551Y+CdXSaZ7r28nPO/gucj45OBT57M3J5Su9djXvwNEJ4Hx1nOH6g3gDaqbl9V6a
o8wAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(surface_trim.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::surface_trim.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAUlJ
REFUSImtVTuyxCAMk3ZShPsflnR6BT+DvWGTeapCMLKwBQAeWRDUh5IUh5iZ7AflR/mmXSxJ
5BVkRk/rMnoFNgcAsGklWDkEkjwAkCRaZMg+6G5UMZptguZxAgAchr1PkRegNGKZIamJjJWN
vuRp24df0HiyBLnd3uxh1hmrWeLxQ6+2+JJnLSlQy2qFlSDpHGzVW9JZvxOPptMFKpkdjISf
lS1q2i2mpglaDeYaR5IlYHjY4uP+YPRaQRvCBTFKFb8u8Mc1hXEWq56tTcKAvQWfwyZicH4s
vAetgyKsx8YkGGR3JFOxtRS6XxXF1zaBWleiW2yyb1TpaudI2P7WU5rPw6q6Kz+HmKre7iBA
Dao3+iDIcyJDHuGhi+Z3gLwc+bqDB0etLLdN3jnqRYI9StLhyJcHrZQqfs4nE/zwmATsluDF
+v/FH1YL6FQPuoNWAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(collapse_alltypes.png) {2006-02-14 23:39:16}

image create photo ::icon_chooser::images::collapse_alltypes.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAr1J
REFUSInNlc1vDGEcxz/PZvqiWta2RAn/AHGoiJeIP8HFS128XsRNuFQaPeFQyYaD9KSRijg0
kUgk0psQRGijLkjY3RYpdrddOrS7OzNfh9ox051WSyK+l5n5vXx/v9/zfJ9nIAADIKGxMWEB
lI+UuL0lD5FhFYN6ruT8D9/T1eUolFfcNaXSyWllMlLftWzYuWhodhd+8Z9GpVKeAGIA4+kR
AL7mU9TUzATHAMxTC4CmhGhtNWH+wDP0vvjeU6mUBgYGFGTu6nKUzXqhCmZWnhRRy8xEGX8I
ZaWJs2PKfxxhvJTB7ZymuzvPpw8jTGRHePs2gkVC3kOfXy9fScPDCg7+l3u0WMxVrsoe8z0B
AfT3Oxoa8n7bsy+QCsH1687sXVZVQiWpUvH48XKIyIQzomcx5lecP0P2fVqTx/Kobyan92oB
u5Dm0EE7nAxQ2lzS5JMPkItB2pBbvo6W5lGSl+OsbHE5eSoRrlA7WGviHesxgzESreuQgcQq
h3jcY8PGZdEr5AVUKaHhF9KbN17UmYg0zLsHUc5/LLQ5ENzQyI56euZXyoKKFAoFJZNJZTKZ
SIK7d13t3VvSnTvun51Bz/NULBbluq4v3aBUCUza2VlWe3tJr197wcnnlnZwkqi7KYif9xQS
FAqirg4aGny66pcK3OeuymeKUAt1HUugDagNh42Oit7eKWpqwCnD6DuLixctmptjVXxWiPyq
q8lnWdwb0+AY7JUesXsWVn+Opp0J8gfg9JEya9d8pbFRfPtuaN9n8+hxPRMTichJ/QIak9xb
Lk2fW7Bv5igfncbYBm1yceLfcfbUkzNLadv+jROHv+B5kLy8gtbVDrZt/CWrWsoIm9zzLlO3
J6nf0Yh1yarEKJX2uHDOZtvWIhLcf7CEffsb2L07hhT+f8yF+aQ3W02L5liIpv+PKyGIH1NU
nxBHGqO/AAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(intersect_line_surface.png) {2006-02-16 13:39:46}

image create photo ::icon_chooser::images::intersect_line_surface.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAARZJ
REFUSInlVUsWAyEIC/N6/zPNnCxdoCMI6HTX17LoRzQECAoYE/0iSYSfk6mHBEkSfhv7nglZ
j4kISJrVTy3ghxj3oogcNrDNQBSLBAH1ScqpqMED63C0/125qCRDFsNPSlueuyH2ZE8kUJgS
myKMajUwZiy+1KZ0myQTM30gtBqmmEn7Y5XmeALXvCMj5xR3F1axzIFyyrIcEirNqxtGyFeJ
VGg9TVq8ayeNcXckJP/R3J1TrHe3E4LvT34zFUADcBysddOpGfVLtmVhDXvX736/wGtwmYFl
bedpOTT3cOqpZBI41E23CgC4risCzkwl7Yk+h+1hXNp5nu3p3M9/qpb70z4zhilZz/qDAElm
EXAj6V+yN1N3xkp1YpHyAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(surface_divide_split.png) {2006-02-15 23:19:26}

image create photo ::icon_chooser::images::surface_divide_split.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAS5J
REFUSInNVUlywzAMAzz5/5eRgyRuWpy60054cGJKALhJBoIRgCQBbC/2nk2QhOY/bIO76C6p
6QBk2T4RfGKGesUwGjejeFuR0NMjrx4cLdBFCM+iGkZnkD3NaUXwUoSa9Fw0oowFBcgGupIL
pOWLWuIGv1BslIqgKYT6/a4A/2oSWrxTzKPktBlhag5YanUBfdDDbIzmwUcpAVhZbK/9Osj6
QHrTWBBRIzVuyGjMEMtCBaCMRq7RZjQ6DANJuximNJ1qu/KVtk/kbDHNI8frhiee6cL9WWxx
jNGGtmDTtRFc4TKq91AVSBEb367X8XAp+1Yimzz9pM4kqyHei9wUMt+0T0SWJycg4vfIrpJz
SP6M0CNknklCGylmynXnq8Lyb6bdzN2PBR4fzL+zN3c8iRfkDEDtAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(folder.png) {2013-01-23 12:17:04}

image create photo ::icon_chooser::images::folder.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAABC5J
REFUSInVlWtsVFUQx39z93S3u+32bQu2VKg8tEVA5aEfjBGUkIhYQ0x8xBhMBL9gYtREiAYa
sV8MIAYTBJGH0RhRESKRACIJKEqD5f2ytjS0YNmWlqWF3b17z/ihD7Y0DcHEGP/JzT1nzsyZ
OTP/MwdSUVFWVFm3f5lu+3COGhHBn5mj9zzyUD8lpG7/Mj1Vc4zWtugZaTq+3lpraTp5Fqko
K5q5a9vqracO/ILjJKs8TwmEc1cZwLgJV4aNLkLEtxggLRhezKzpEzOAEanf7JkPZjIYJNKw
07qxLqk7uG+xqpJ0LdNeWl4ljb+v0WSiExHh5G+HcRyHRNJuMmnBAG78CuBw1+QJALz+1poa
U3/4hN4xdoyoqgKIiFSMH7NceiK5EQ2DRjVotBVlRZXfb1ryrXutHUDU8xS14HP6lAKhfHn8
6Xcqj9e3bDEAPuMjEomSTCT6lPKH5IJaAOJXInyzdv7mxiY9YgA6Wi/JsPJ7iXdFUFUZEIYo
oezh8mrVwmlGRBTj1+jliwAC6AADx8iKRe+/u2PPobZbPvS/Dzn206c2b0ieiIj6M/JYsuDt
6qdmTHJVhWQyCYDFaXtsztKVAHJk5yqbHowK0F1tLymIAw4KCGo1nriNcY/OdQDk0I6Vmhm+
nn8Rob2lnZbGc30ytUI8Hj+YVVz6jFHrKTjScq5ZRQVrrRQMLSA8bmS3h540ZxSMuf/hGfNy
jVpXzv95kdK7R0kvb6Hbe69TsEQvx7pemTPzsEEsWfmFKB4+f7DHAFQVEUFVcRyHTeu/qi0Z
XeYa91qXDp8wUeJdkT5+9xioiIiqqj89i+/2HJv/64qtalrONmpuabGo9QbNfWtzvYZzwqd7
4ywB0m5SLxdouonOfwT5qOrFyXcOzf46GAqXAP2aAeJotPWvE4lQ1oTZLy9N/iMHFWVFlV+u
WbI5M9uS6GzttxjMLsYES1EG3NmBGzlC85naq9UfbJy1eVftj71yA2C9pIr4cPzp/U4Qu9qq
XG3tLb2KiLixuPZyhxRmmkCmOELQp97QHz557fZAMNg59fnqqAFQdRGfH8eYPuKlbIJ6StuF
CLHOKGOnPtfXq25Efims3jBlo/GHpOFITbx8ROGTPQ48fGk+Lp3vQMkgr3gkQjerAVRihHNj
5BTmEGnYez0tKqhov7Faj8yCsWz4Ysvn9c2XartTZJNYC5cjVxj9wHjinRew1k1NMP5gqPsK
pgafOncAtYSyS/h55+5T2/cc+HjRG8+2me5Vy9mjJxg5aTLW7cTxGfCZmxb2RphAiPaWCNt3
7VtXXl52dEH1Z54jInS0nJfsgiJJC/gQnycmPV1Mejo9/9RxqqzfPC0YxLMip2tr/9h7tHnt
E9MnxQCkMDfzvmlTRs1T1YDagW/ILUGI7a6pW/fewhdq5r65avBm8r/C38w5z1GGtsRLAAAA
AElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(dimension_dist.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::dimension_dist.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAHRJ
REFUSIljYCADrNq/6j85+mgP/kMxMp+B4cCdA1jdy0QHF8EcAg0zXOGGKo7LwdgAyZ5gRFhJ
yBJGAvJDE+CMA+Q4ol3agEUuyTbQXgMjAzx0iEoawzN9DDWAHAvULpdHY3iQAKRSG0OMXIBi
Fh0rw6EKAO4YI1Ay/rkMAAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(folder.png) {2013-01-23 12:17:04}

image create photo ::icon_chooser::images::folder.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAABC5J
REFUSInVlWtsVFUQx39z93S3u+32bQu2VKg8tEVA5aEfjBGUkIhYQ0x8xBhMBL9gYtREiAYa
sV8MIAYTBJGH0RhRESKRACIJKEqD5f2ytjS0YNmWlqWF3b17z/ihD7Y0DcHEGP/JzT1nzsyZ
OTP/MwdSUVFWVFm3f5lu+3COGhHBn5mj9zzyUD8lpG7/Mj1Vc4zWtugZaTq+3lpraTp5Fqko
K5q5a9vqracO/ILjJKs8TwmEc1cZwLgJV4aNLkLEtxggLRhezKzpEzOAEanf7JkPZjIYJNKw
07qxLqk7uG+xqpJ0LdNeWl4ljb+v0WSiExHh5G+HcRyHRNJuMmnBAG78CuBw1+QJALz+1poa
U3/4hN4xdoyoqgKIiFSMH7NceiK5EQ2DRjVotBVlRZXfb1ryrXutHUDU8xS14HP6lAKhfHn8
6Xcqj9e3bDEAPuMjEomSTCT6lPKH5IJaAOJXInyzdv7mxiY9YgA6Wi/JsPJ7iXdFUFUZEIYo
oezh8mrVwmlGRBTj1+jliwAC6AADx8iKRe+/u2PPobZbPvS/Dzn206c2b0ieiIj6M/JYsuDt
6qdmTHJVhWQyCYDFaXtsztKVAHJk5yqbHowK0F1tLymIAw4KCGo1nriNcY/OdQDk0I6Vmhm+
nn8Rob2lnZbGc30ytUI8Hj+YVVz6jFHrKTjScq5ZRQVrrRQMLSA8bmS3h540ZxSMuf/hGfNy
jVpXzv95kdK7R0kvb6Hbe69TsEQvx7pemTPzsEEsWfmFKB4+f7DHAFQVEUFVcRyHTeu/qi0Z
XeYa91qXDp8wUeJdkT5+9xioiIiqqj89i+/2HJv/64qtalrONmpuabGo9QbNfWtzvYZzwqd7
4ywB0m5SLxdouonOfwT5qOrFyXcOzf46GAqXAP2aAeJotPWvE4lQ1oTZLy9N/iMHFWVFlV+u
WbI5M9uS6GzttxjMLsYES1EG3NmBGzlC85naq9UfbJy1eVftj71yA2C9pIr4cPzp/U4Qu9qq
XG3tLb2KiLixuPZyhxRmmkCmOELQp97QHz557fZAMNg59fnqqAFQdRGfH8eYPuKlbIJ6StuF
CLHOKGOnPtfXq25Efims3jBlo/GHpOFITbx8ROGTPQ48fGk+Lp3vQMkgr3gkQjerAVRihHNj
5BTmEGnYez0tKqhov7Faj8yCsWz4Ysvn9c2XartTZJNYC5cjVxj9wHjinRew1k1NMP5gqPsK
pgafOncAtYSyS/h55+5T2/cc+HjRG8+2me5Vy9mjJxg5aTLW7cTxGfCZmxb2RphAiPaWCNt3
7VtXXl52dEH1Z54jInS0nJfsgiJJC/gQnycmPV1Mejo9/9RxqqzfPC0YxLMip2tr/9h7tHnt
E9MnxQCkMDfzvmlTRs1T1YDagW/ILUGI7a6pW/fewhdq5r65avBm8r/C38w5z1GGtsRLAAAA
AElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(align_text_left.png) {2006-02-15 11:38:42}

image create photo ::icon_chooser::images::align_text_left.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAH1J
REFUSIntklEKgCAMQJ/hAXb/M0rsAFF9GWYmmAYmvS9l7rmNQYABUFUAJm5R1fSzk8AjIjlZ
KeV/lAfaYfwhrjhGRIAHNb2f0CHJJUjR8VhHSBgAE17sZLcamZud8St3OGuEXpqLN+0AYFmX
s7NWGHcQj+jv4MJ4HXyfHX9oNFcTr690AAAAAElFTkSuQmCC
}

set ::icon_chooser::ImportedDates(PantallaInfo.png) {2010-02-12 12:57:44}

image create photo ::icon_chooser::images::PantallaInfo.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAUJJ
REFUSInVVDGOgzAQHCO+gFzSnKLkDTzhanS/SEWX9tKl4hcoNU/gDdFZNJQoj9gUGGThXYNP
ueKmQcDO7K5n14qIMCMFgFx3VNYZFBEh1x0BmF5mJFVjCABVjaHEcuCFuUgAoGoMWUHKdUdV
Y0gREU5X49HCUhxSWxbu5ycAoKwztL1lzB89qbLOvB/q+P3DZherEosCpsb3BLe9JbjgCndr
TqTgYSzkktZK9/MTue7404htOpogevo2Qnq6Gvr82Bd8+zoo1odhLJZZX58U68MczGGzh7Xz
LEFyWSTkultIwR7KOlsCxNGQFoXD43JQfz9LsYgepVikANgr7h14XJxR3Zrv9UBxjrpo++np
7UJIFMCyJ+6/UDLRA04c4PdKig0m4DCMRXCtohJwbUu3wq+OaG+SLbODJu8R2ML/3+QXv7iZ
8PzWLHQAAAAASUVORK5CYII=
}

set ::icon_chooser::ImportedDates(border_vertical.png) {2004-11-18 21:53:33}

image create photo ::icon_chooser::images::border_vertical.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAJhJ
REFUSIntlGEKwyAMRp+lB5Di/U9YSg5QZn+lRJk6x1pK1weikM9EP4NgcAAiAsBAERFhHMYo
IqksSRCmEOdldo1kvSQ1FO99uUZ/4Hc4XeQntu4o1VuUOH7DBXnbBCVL4QyX/nLD3dBfrTZ/
orE5XV4g71FLmEIEaGnW17rHD3+1c9ui5bEdNY3N+bzBdwXUinzu1dyDDT6jfveUsjRKAAAA
AElFTkSuQmCC
}

set ::icon_chooser::ImportedDates() {2017-01-25 21:00:12}

image create photo ::icon_chooser::images:: -data {
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAAH6ji2bAAAABGdBTUEAAYagMeiWXwAAA6hJ
REFUOI2VlFtsVFUUhv+99zlnbmduLYWhaWm1jSkgVO2UUmkrRSlWEmK0EYiXWqnGJoZETTEq
RhOMUUzENCZTHiQixBJvD2BoiTZEpkSmnoFpdQxWCA3pnZZOZ5jO5cw524dKw0A7wf9x51/f
2uvPyiKcc8zL4/FgxRMNfpmxcgFqEldPfl8GACTN5jzUxUFMZst8EQAIq8vcZqCQrS1/ULhp
FQqKH+Ym8yUpGk6mQ25XYUNTOH/vAS44spxz5YSQWFmVVV2+AsIDlWUEAExPNb6tjw8HuSCq
mvdUJ5xfdk7fimGMAW1tbXc2cOTko2LT01hXuxuEOQFAEEUxzSNV1W22+r3LUm/sOyIND4Iy
VCZF6axwO4zmuIoTTW/WyMODoN6uQPTaWGcsmTgzn4btvbZfhOJV7oVGZYyB+Pv/MPz0w3cJ
wu6Az4sQgsyp3Wpe4M2Ites3wmTWMTnuwuXg1729vaDpZZTK1VueczHscM2Gt8tToyPihs1b
KdfTUUJJaSkIgevTw9xx9DTPe/LZYHZ1Xbu/r9+URhQZyy18cbeqWmyApoEmEvr1wUsfQTJI
acb4P8GfQ7XbZiRv11+W919txxr3/ayoZCW/MTOT1jrrWA8nss2cu/WZ02J5zQax5vHHiCiZ
FEX5bwbZ5nAePDG2WDSKogCkeFWVdc/+45kyVBQFZF/Hj+v53/3nQOmixubmZpCBgQH4fL5M
QEDXAXAQk2wEpZTKVqsWuj4BLcWRUgHKQAhBfX09hO7ubrS0tCwOI4Q6KjftEa32gukrA19o
4dC4oGnm3NLy/gihpyLxeJfa5wvQ8PRkr9+PxTcLALHa7dK6Rx6NBM8fXVJw7ydZ73z2p/HX
k9H4GrclEZuF+u4rz6uMBVjFxlo6OT6ElHo5I5Dl3VNIY1HutMiNQmlFbTyZQKSm3gJdB89a
BrHDeyG/7xz472fORvIKv4JkuLp4wgD41MQUKNXjFusVLttlS+A3cFGCGIuCmswwHvNcGPq4
lYxOTb6VCIcu8thsIiOQFa2sNO5q7aANu166dvybajI2NGz+fO+3E681kMiOqiwyMz2ydNvO
EzwaCSeVnh5iMC6QIROo5eXWg4bqLS/caPtg+3RTnRkAmNWePeKL1vLspTZx9UNurqXU0YCv
UQ+HpqDrc5eJc8DjaZ/jLM8vchzouGjffzhAl7jyMv18wWkYg6IoELQc13348JBPO99zJPT6
zpL/C7opTdOgadrdn5y71b8fv2qYnbMAoQAAAABJRU5ErkJggg==
}

set ::icon_chooser::ImportedDates(dimension_arc.png) {2006-02-14 21:00:15}

image create photo ::icon_chooser::images::dimension_arc.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAAZ9J
REFUSIm1VEuKwkAQfR0HYkgLLtwJRrfiwhN4Ag/htTyB53CjCLrXjaAiIkRQMMHYn5qFI5Nx
MjGd0YIQkrxfVXcaeGWx+MN+v6coilCpVIxEiJ5jpJS0WCxIKUVCCDqdThlYdwOjOD/qR4OD
wYDa7TaEENjtduh0OlgulygWi7BtG7PZLLd3liGksDM4pTqkCcSJRP+NCiA21slkQlprcM6x
2WyglMLhcECv12P9fp+63S6Gw+E3IQgCcl2Xrddr0lrDcRzM53M0m00cj0d4ngff9zM391gv
aC6PeuJ7KwXIEr4lEhKBJpUULdeUHvdOskhiXiEEXS4XuK4Ly7KYlJKu1+uNwG4Ux3FYGIZU
KBRg2zZTStH5fEapVEIURVBKgTGWbDAej0lrDc/zsN1uEYYharUagiBAuVzGarVCvV7HaDRC
q9UC5xy+74NzDiklGo0GptMpqtXq3/3HLmitSQhhtB9SiuLiWcgmBvlW+13Cd/GXGFDsbhTm
wwSMHL9hVoM04XhHv3BfhwVj7z704oEMjJ5jM8z0UYQZrcMn8u/aBO0pqfgAAAAASUVORK5C
YII=
}

set ::icon_chooser::ImportedDates(zframe.png) {2006-03-01 12:46:29}

image create photo ::icon_chooser::images::zframe.png -data {
iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAAGXcA1uAAAABGdBTUEAAYagMeiWXwAAA0VJ
REFUSIm1k09oVFcUxn/3zYAOdmJNkwY7KUNH8x6TJlKRBrpJKbQ1G0koLW6aLBS6aOnKRXFR
F5mFaynJRrupA0Wk0mClCQmJI4pYsZ1FnAlJNTF0qJP3MlOjdRwnL6eLyXuT+ZNMU+wHFy73
fud89zv3HCUiONAAfl14JIZhiBIRJoYHxA71oapoU8llAdCSU+fk+p3p4pWIICIMnr8qz1fX
ZPD8VRGR0sW52KKc+D4huq6XLhKTZ8UhOKtMcSM0Z3Pl1pwU7BJLA4hEY7JrdyPXZzPMpf8W
AK/DGJ+2KNjC/t0F2no6ixFfffgyO589xF9I81lPpyoyh/pl9Ocf/v2rNoNWeWA9fi6GYciC
lRPDMMQwDIlEY27WMoVINCb+5iB/ZHM4x7t2eOhs9fNJ115VpfDOyrc0aE8Iv/YSrzf6eKPZ
x8FgAzOzsyWSU/J8Oi6JybM1jW5c7qYe8cVVadsBkWhMCrZIwRa5cmuurKSu6Y3dqeu66Lou
8+ZT0XVdzJW81DQtIsw+fCK6rstkwpLJhCVTSUvuzP/ltreIlLoMIP5ghaOnL/NTfAkApaB1
j29zD9mMxYWTR/BoCqVgh1fD7/OUBbhlnRgekJsNx+kIG6Syz1i1we/zoPKPOPZBh3ICvA65
++MTBO7eZtFcJmu9AsBjoKtpGegoSYwP9Us+HZd8Oi7jQ/11f9trh/q4d/c2qeQ13v/8O0U9
bKeP/lMvbRfe+pQiqlpqHV9/+u6Wtrd0EInG5K39e3nv4D7stXXbTiCglMKjKaZ+u0f89z9r
im06PZFoTMLBFtpDrRw60M6CmWPezDG/tL7MHAtmjkMH2mkPtRIOttR0uWWJ2oIBUpk8Zy7e
YNHKURmtgDMXb5DK5GkLBkg+SFfl2FLgvvmUVDbHN18cpndwBHutXMKjKUZO9fLl0BiBihmr
K9DVtMwvyRk6wgYAI6d6ATh6+jIAF04ecbmBPTuZTs6sT02Fy8pPnhgekEC4G4BU8ho3G44D
0NLUSHPzq6zaa8WXeTRMc4m0lXEfdLjno6pPLhv+QLibfW++7V5unMCx0Us1261W0iqByuRO
YjvUVzdBPbgOxkYvief+jwAvJHGVwP+FfwCzjWJXaDcuRwAAAABJRU5ErkJggg==
}

array set preferences {toolbar_one_col 1}