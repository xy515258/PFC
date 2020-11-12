# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
filenames = []
for i in range(0,41):
	filename = "./data_"+str(i*5000).zfill(8)+".vtk"
	filenames.append(filename)
	
print(filenames)

data_000 = LegacyVTKReader(FileNames=filenames)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2150, 1174]

# show data in view
data_000Display = Show(data_000, renderView1)

# get color transfer function/color map for 'n1'
n1LUT = GetColorTransferFunction('n1')

# get opacity transfer function/opacity map for 'n1'
n1PWF = GetOpacityTransferFunction('n1')

# trace defaults for the display properties.
data_000Display.Representation = 'Slice'
data_000Display.ColorArrayName = ['POINTS', 'n1']
data_000Display.LookupTable = n1LUT
data_000Display.OSPRayScaleArray = 'n1'
data_000Display.OSPRayScaleFunction = 'PiecewiseFunction'
data_000Display.SelectOrientationVectors = 'None'
data_000Display.ScaleFactor = 191.9
data_000Display.SelectScaleArray = 'n1'
data_000Display.GlyphType = 'Arrow'
data_000Display.GlyphTableIndexArray = 'n1'
data_000Display.GaussianRadius = 9.595
data_000Display.SetScaleArray = ['POINTS', 'n1']
data_000Display.ScaleTransferFunction = 'PiecewiseFunction'
data_000Display.OpacityArray = ['POINTS', 'n1']
data_000Display.OpacityTransferFunction = 'PiecewiseFunction'
data_000Display.DataAxesGrid = 'GridAxesRepresentation'
data_000Display.SelectionCellLabelFontFile = ''
data_000Display.SelectionPointLabelFontFile = ''
data_000Display.PolarAxes = 'PolarAxesRepresentation'
data_000Display.ScalarOpacityUnitDistance = 17.574109285000795
data_000Display.ScalarOpacityFunction = n1PWF
data_000Display.IsosurfaceValues = [0.642685]
data_000Display.Slice = 959

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
data_000Display.DataAxesGrid.XTitleFontFile = ''
data_000Display.DataAxesGrid.YTitleFontFile = ''
data_000Display.DataAxesGrid.ZTitleFontFile = ''
data_000Display.DataAxesGrid.XLabelFontFile = ''
data_000Display.DataAxesGrid.YLabelFontFile = ''
data_000Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
data_000Display.PolarAxes.PolarAxisTitleFontFile = ''
data_000Display.PolarAxes.PolarAxisLabelFontFile = ''
data_000Display.PolarAxes.LastRadialAxisTextFontFile = ''
data_000Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [10000.0, 959.5, 959.5]
renderView1.CameraFocalPoint = [0.0, 959.5, 959.5]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
data_000Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [10000.0, 959.5, 959.5]
renderView1.CameraFocalPoint = [0.0, 959.5, 959.5]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 926.8068527402393

# save animation
SaveAnimation('./test.avi', renderView1, ImageResolution=[2148, 1172],
    FrameRate=3,
    FrameWindow=[0, 40], 
    # FFMPEG options
    Quality='1')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [10000.0, 959.5, 959.5]
renderView1.CameraFocalPoint = [0.0, 959.5, 959.5]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 926.8068527402393

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
