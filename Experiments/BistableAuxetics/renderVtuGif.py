import glob
import os
import sys
import tempfile

from paraview.simple import (ColorBy, CreateRenderView, Delete, GetColorTransferFunction,
                             GetScalarBar, PointDatatoCellData, Render, ResetCamera,
                             SaveScreenshot, Show, Text, Threshold, WarpByVector,
                             XMLUnstructuredGridReader)
from PIL import Image

FRAME_SIZE = [800, 800]
CAMERA_MARGIN = 1.12
FRAME_DURATION_MS = 90
LAST_FRAME_DURATION_MS = 1800
SHEET_COLUMNS = 4
SHEET_ROWS = 2
SOLID_DENSITY = 0.5
MEDIUM_COLOR = [0.6, 0.6, 0.6]


def load_warped_displacement(file_name):
    reader = XMLUnstructuredGridReader(FileName=[file_name])
    warp = WarpByVector(Input=reader)
    warp.Vectors = ['POINTS', 'fValues']
    warp.ScaleFactor = 1.0
    return reader, warp


def delete_sources(*sources):
    for source in sources:
        Delete(source)


def create_view():
    view = CreateRenderView()
    view.ViewSize = FRAME_SIZE
    view.UseColorPaletteForBackground = 0
    view.Background = [1, 1, 1]
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1
    return view


def compute_camera_and_max_displacement(view, last_file):
    reader, warp = load_warped_displacement(last_file)
    warp.UpdatePipeline()
    u_max = warp.GetPointDataInformation().GetArray('fValues').GetRange(-1)[1]
    representation = Show(warp, view)
    ResetCamera(view)
    view.CameraParallelScale *= CAMERA_MARGIN
    camera = (list(view.CameraPosition), list(view.CameraFocalPoint), view.CameraParallelScale)
    delete_sources(representation, warp, reader)
    return camera, u_max


def create_strain_label(view):
    label = Text()
    representation = Show(label, view)
    representation.Color = [0, 0, 0]
    representation.FontSize = 18
    representation.WindowLocation = 'Upper Left Corner'
    return label


def has_density(source):
    return 'density' in source.PointData.keys()


def select_solid(source):
    cells = PointDatatoCellData(Input=source)
    cells.PassPointData = 1
    solid = Threshold(Input=cells)
    solid.Scalars = ['CELLS', 'density']
    solid.ThresholdMethod = 'Above Upper Threshold'
    solid.UpperThreshold = SOLID_DENSITY
    return [solid, cells]


def show_mesh(view, warp):
    representation = Show(warp, view)
    representation.Representation = 'Wireframe'
    representation.ColorArrayName = ['POINTS', '']
    representation.AmbientColor = MEDIUM_COLOR
    representation.DiffuseColor = MEDIUM_COLOR
    return representation


def show_deformed(view, warp, u_max):
    if not has_density(warp):
        return [show_displacement_magnitude(view, warp, u_max)]
    solid, cells = select_solid(warp)
    representation = show_displacement_magnitude(view, solid, u_max)
    representation.Representation = 'Surface With Edges'
    representation.EdgeColor = [0.15, 0.15, 0.15]
    return [show_mesh(view, warp), representation, solid, cells]


def show_displacement_magnitude(view, warp, u_max):
    representation = Show(warp, view)
    representation.Representation = 'Surface'
    ColorBy(representation, ('POINTS', 'fValues', 'Magnitude'))
    lut = GetColorTransferFunction('fValues')
    lut.RescaleTransferFunction(0.0, u_max)
    bar = GetScalarBar(lut, view)
    bar.Title = '|u| [mm]'
    bar.ComponentTitle = ''
    bar.TitleColor = [0, 0, 0]
    bar.LabelColor = [0, 0, 0]
    representation.SetScalarBarVisibility(view, True)
    return representation


def set_camera(view, camera):
    view.CameraPosition, view.CameraFocalPoint, view.CameraParallelScale = camera


def render_frame(view, file_name, strain, label, camera, u_max, png_name):
    reader, warp = load_warped_displacement(file_name)
    shown = show_deformed(view, warp, u_max)
    label.Text = 'strain = %.3f' % strain
    set_camera(view, camera)
    Render(view)
    SaveScreenshot(png_name, view, ImageResolution=FRAME_SIZE)
    delete_sources(*shown, warp, reader)
    return Image.open(png_name).convert('RGB').copy()


def render_frames(files, strains):
    view = create_view()
    camera, u_max = compute_camera_and_max_displacement(view, files[-1])
    label = create_strain_label(view)
    with tempfile.TemporaryDirectory() as tmp:
        return [render_frame(view, f, e, label, camera, u_max, os.path.join(tmp, 'f%03d.png' % i))
                for i, (f, e) in enumerate(zip(files, strains))]


def save_gif(frames, file_name):
    palette = [f.convert('P', palette=Image.ADAPTIVE, colors=128) for f in frames]
    durations = [FRAME_DURATION_MS]*len(palette)
    durations[-1] = LAST_FRAME_DURATION_MS
    palette[0].save(file_name, save_all=True, append_images=palette[1:], duration=durations, loop=0)


def save_frame_sheet(frames, file_name):
    width, height = FRAME_SIZE[0]//2, FRAME_SIZE[1]//2
    n = SHEET_COLUMNS*SHEET_ROWS
    sheet = Image.new('RGB', (SHEET_COLUMNS*width, SHEET_ROWS*height), 'white')
    for k in range(n):
        frame = frames[(k+1)*(len(frames)-1)//n].resize((width, height))
        sheet.paste(frame, ((k % SHEET_COLUMNS)*width, (k//SHEET_COLUMNS)*height))
    sheet.save(file_name)


def main(prefix, strain_step):
    files = sorted(glob.glob(prefix + '_u_step*.vtu'))
    strains = [strain_step*i for i in range(len(files))]
    frames = render_frames(files, strains)
    save_gif(frames, prefix + '_open.gif')
    save_frame_sheet(frames, prefix + '_open_frames.png')
    print('saved', prefix + '_open.gif', 'frames', len(frames))


if __name__ == '__main__':
    main(sys.argv[1], float(sys.argv[2]))
