import sys
sys.path.append("/usr/lib64/paraview/python3.10/site-packages")

colormap          = sys.argv[1]
pv_fileout_prefix = sys.argv[2]
pv_specular       = sys.argv[3]
pv_light_intensity= sys.argv[4]
caption = ''
zoomfactor = 1.4 # camera zoom in by this factor

structuredvtsname = pv_fileout_prefix+'.vts'
pdffile           = pv_fileout_prefix+'.pdf'

# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Structured Grid Reader'
structuredvts = XMLStructuredGridReader(FileName=[structuredvtsname])
structuredvts.PointArrayStatus = ['scalar']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [800, 800] 

# show data in view
structuredvtsDisplay = Show(structuredvts, renderView1)

# trace defaults for the display properties.
structuredvtsDisplay.Representation = 'Surface'
structuredvtsDisplay.ColorArrayName = [None, '']
structuredvtsDisplay.OSPRayScaleArray = 'scalar'
structuredvtsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
structuredvtsDisplay.SelectOrientationVectors = 'None'
structuredvtsDisplay.ScaleFactor = 0.1
structuredvtsDisplay.SelectScaleArray = 'None'
structuredvtsDisplay.GlyphType = 'Arrow'
structuredvtsDisplay.GlyphTableIndexArray = 'None'
structuredvtsDisplay.GaussianRadius = 0.005
structuredvtsDisplay.SetScaleArray = ['POINTS', 'scalar']
structuredvtsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
structuredvtsDisplay.OpacityArray = ['POINTS', 'scalar']
structuredvtsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
structuredvtsDisplay.DataAxesGrid = 'GridAxesRepresentation'
structuredvtsDisplay.SelectionCellLabelFontFile = ''
structuredvtsDisplay.SelectionPointLabelFontFile = ''
structuredvtsDisplay.PolarAxes = 'PolarAxesRepresentation'
structuredvtsDisplay.ScalarOpacityUnitDistance = 0.043296791257339616

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
structuredvtsDisplay.DataAxesGrid.XTitleFontFile = ''
structuredvtsDisplay.DataAxesGrid.YTitleFontFile = ''
structuredvtsDisplay.DataAxesGrid.ZTitleFontFile = ''
structuredvtsDisplay.DataAxesGrid.XLabelFontFile = ''
structuredvtsDisplay.DataAxesGrid.YLabelFontFile = ''
structuredvtsDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
structuredvtsDisplay.PolarAxes.PolarAxisTitleFontFile = ''
structuredvtsDisplay.PolarAxes.PolarAxisLabelFontFile = ''
structuredvtsDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
structuredvtsDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(structuredvtsDisplay, ('POINTS', 'scalar'))

# rescale color and/or opacity maps used to include current data range
structuredvtsDisplay.RescaleTransferFunctionToDataRange(True, False)


# get color transfer function/color map for 'scalar'
scalarLUT = GetColorTransferFunction('scalar')

# get opacity transfer function/opacity map for 'scalar'
scalarPWF = GetOpacityTransferFunction('scalar')

ImportPresets(filename=colormap)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
scalarLUT.ApplyPreset('inertial', True)

# hide color bar/color legend
#structuredvtsDisplay.SetScalarBarVisibility(renderView1, False)

# uncomment to show color bar legend --------------------------------------

# uncomment to show color bar/color legend
#structuredvtsDisplay.SetScalarBarVisibility(renderView1, True)

# get color legend/bar for scalarLUT in view renderView1
#scalarLUTColorBar = GetScalarBar(scalarLUT, renderView1)

# Properties modified on scalarLUTColorBar
#scalarLUTColorBar.AutoOrient = 0
#scalarLUTColorBar.Orientation = 'Horizontal'
#scalarLUTColorBar.Title = ''
#scalarLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
#scalarLUTColorBar.LabelFontSize = 12

# change scalar bar placement
##scalarLUTColorBar.WindowLocation = 'UpperRightCorner'
#scalarLUTColorBar.WindowLocation = 'AnyLocation'
#scalarLUTColorBar.Position = [0.515, 0.695]
#scalarLUTColorBar.ScalarBarLength = 0.200000000000002
# uncomment to show legend end ----------------------------------

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# Properties modified on Specular
structuredvtsDisplay.Specular = float(pv_specular)

# create a new 'Text'
text1 = Text()

# find source
structuredvts = FindSource('structured.vts')

# Properties modified on text1
text1.Text = caption
#text1.Text = """Dissipation\n(-0.3,0.4)\nE=1e-6"""

# axis properties: ------------------------------------------
# Properties modified on renderView1.AxesGrid
#renderView1.AxesGrid.Visibility = 1

#renderView1.AxesGrid.XTitle = 's'
#renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
#renderView1.AxesGrid.XTitleFontSize = 20
#renderView1.AxesGrid.YTitle = 'z   '
#renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
#renderView1.AxesGrid.YTitleFontSize = 20
#renderView1.AxesGrid.ZTitle = ''
#renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]

#renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
#renderView1.AxesGrid.XLabelFontSize = 16
#renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
#renderView1.AxesGrid.YLabelFontSize = 16
#renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
#renderView1.AxesGrid.ZLabelFontSize = 16

# Properties modified on renderView1.AxesGrid
##renderView1.AxesGrid.XTitleFontFamily = 'Times'
##renderView1.AxesGrid.XLabelFontFamily = 'Times'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.AxesToLabel = 3

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.CullBackface = 1
renderView1.AxesGrid.CullFrontface = 0
# axis properties end ---------------------------------------

# Create a new 'Light'
light1 = AddLight(view=renderView1)
light1.Intensity = 0.5*float(pv_light_intensity)

# Create a new 'Light'
light2 = AddLight(view=renderView1)
light2.Intensity = 0.5*float(pv_light_intensity)
light2.FocalPoint = [0.5, 0.5, 0.0]
light2.Position   = [1.5, 1.5, 1.0]

light3 = AddLight(view=renderView1)
light3.Intensity = 0.5*float(pv_light_intensity)
light3.FocalPoint = [1.5, 1.5, 0.0]
light3.Position   = [0.5, 0.5, 1.0]

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4999975, 0.5001943495, 2.8800003795972564]
renderView1.CameraFocalPoint = [0.4999975, 0.5001943495, -0.12897189432561812]
renderView1.CameraParallelScale = 0.500


# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1037, 652]

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# get layout
layout1 = GetLayout()

# show data in view
text1Display = Show(text1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on text1Display
#text1Display.WindowLocation = 'UpperRightCorner'

# Properties modified on text1Display
text1Display.Color = [0.0, 0.0, 0.0]

camera=GetActiveCamera()
camera.Dolly(zoomfactor) # zoom in
Render()

# export view
ExportView(pdffile, view=renderView1)
