#!/usr/bin/pvbatch
import sys
colormap          = sys.argv[1]
pv_fileout_prefix = sys.argv[2]
#print('colormap=',colormap)
#print('pv_fileout_prefix=',pv_fileout_prefix)
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

# show color bar/color legend
structuredvtsDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'scalar'
scalarLUT = GetColorTransferFunction('scalar')

# get opacity transfer function/opacity map for 'scalar'
scalarPWF = GetOpacityTransferFunction('scalar')

ImportPresets(filename=colormap)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
scalarLUT.ApplyPreset('inertial', True)

# hide color bar/color legend
#structuredvtsDisplay.SetScalarBarVisibility(renderView1, False)

# get color legend/bar for scalarLUT in view renderView1
scalarLUTColorBar = GetScalarBar(scalarLUT, renderView1)

# change scalar bar placement
scalarLUTColorBar.WindowLocation = 'AnyLocation'
scalarLUTColorBar.Position = [0.63, 0.78]
scalarLUTColorBar.ScalarBarLength = 0.400000000000002

# Properties modified on scalarLUTColorBar
scalarLUTColorBar.AutoOrient = 0
scalarLUTColorBar.Orientation = 'Horizontal'

# Properties modified on scalarLUTColorBar
scalarLUTColorBar.Title = ''

# Properties modified on scalarLUTColorBar
scalarLUTColorBar.LabelColor = [0.0, 0.0, 0.0]

# Properties modified on scalarLUTColorBar
# scalarLUTColorBar.WindowLocation = 'UpperRightCorner'

# Properties modified on scalarLUTColorBar
scalarLUTColorBar.LabelFontSize = 24

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# Create a new 'Light'
light1 = AddLight(view=renderView1)

# Properties modified on light1
light1.Intensity = 0.4

# create a new 'Text'
text1 = Text()

# find source
structuredvts = FindSource('structured.vts')

# Properties modified on text1
text1.Text = caption
#text1.Text = """Dissipation\n(-0.3,0.4)\nE=1e-6"""

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1037, 652]

# get layout
layout1 = GetLayout()

# show data in view
text1Display = Show(text1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on text1Display
text1Display.WindowLocation = 'UpperRightCorner'

# Properties modified on text1Display
text1Display.Color = [0.0, 0.0, 0.0]

camera=GetActiveCamera()
camera.Dolly(zoomfactor) # zoom in
Render()

# export view
ExportView(pdffile, view=renderView1)

# Properties modified on text1Display
text1Display.Color = [0.0, 0.0, 0.0]


