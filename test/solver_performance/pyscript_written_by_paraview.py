# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XDMF Reader'
out_data_sim_0xmf = XDMFReader(registrationName='out_data_sim_0.xmf', FileNames=['/home/masarunagaso/workspace/TomoATT/test/solver_performance/OUTPUT_FILES/out_data_sim_0.xmf'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
out_data_sim_0xmfDisplay = Show(out_data_sim_0xmf, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
out_data_sim_0xmfDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
out_data_sim_0xmfDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'procid'
procidLUT = GetColorTransferFunction('procid')

# get opacity transfer function/opacity map for 'procid'
procidPWF = GetOpacityTransferFunction('procid')

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=out_data_sim_0xmf)

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'

# hide data in view
Hide(out_data_sim_0xmf, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'T_res_src_0_inv_0000'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(procidLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'T_res_src_0_inv_0000'
t_res_src_0_inv_0000LUT = GetColorTransferFunction('T_res_src_0_inv_0000')

# get opacity transfer function/opacity map for 'T_res_src_0_inv_0000'
t_res_src_0_inv_0000PWF = GetOpacityTransferFunction('T_res_src_0_inv_0000')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# Properties modified on slice1Display.DataAxesGrid
slice1Display.DataAxesGrid.GridAxesVisibility = 1

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2329, 1650)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [6086.06543558693, 0.0, 7825.811009542308]
renderView1.CameraFocalPoint = [6086.06543558693, 0.0, -3.410605131648481e-13]
renderView1.CameraViewUp = [1.0, 2.220446049250313e-16, 0.0]
renderView1.CameraParallelScale = 2025.468932642534

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).