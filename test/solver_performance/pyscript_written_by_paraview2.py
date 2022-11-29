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

# Properties modified on out_data_sim_0xmf
out_data_sim_0xmf.PointArrayStatus = ['T_res_src_0_inv_0000']

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

# get color transfer function/color map for 'T_res_src_0_inv_0000'
t_res_src_0_inv_0000LUT = GetColorTransferFunction('T_res_src_0_inv_0000')

# get opacity transfer function/opacity map for 'T_res_src_0_inv_0000'
t_res_src_0_inv_0000PWF = GetOpacityTransferFunction('T_res_src_0_inv_0000')

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=out_data_sim_0xmf)

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

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

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# Properties modified on slice1Display.DataAxesGrid
slice1Display.DataAxesGrid.GridAxesVisibility = 1

# get color legend/bar for t_res_src_0_inv_0000LUT in view renderView1
t_res_src_0_inv_0000LUTColorBar = GetScalarBar(t_res_src_0_inv_0000LUT, renderView1)

# change scalar bar placement
t_res_src_0_inv_0000LUTColorBar.Orientation = 'Horizontal'
t_res_src_0_inv_0000LUTColorBar.WindowLocation = 'Any Location'
t_res_src_0_inv_0000LUTColorBar.Position = [0.33262344353799916, 0.34363636363636363]
t_res_src_0_inv_0000LUTColorBar.ScalarBarLength = 0.3299999999999998

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