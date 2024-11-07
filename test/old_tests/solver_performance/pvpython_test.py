from paraview.simple import *


if __name__ == '__main__':
    # get args
    args = sys.argv
    if len(args) != 3:
        print("Usage: pvpython draw_slice.py <data_file> <dataset_name>")
        sys.exit(1)

    # get data
    fpath = args[1]
    name_dset = args[2]


    #fpath = "OUTPUT_FILES/out_data_sim_0.xmf"
    #name_dset = 'fun_inv_0000'

    # create a new 'XDMF Reader'
    out_data_sim_0xmf = XDMFReader(registrationName=fpath.split('/')[-1], FileNames=[fpath])

    # Limit the data array to be read
    # Properties modified on out_data_sim_0xmf
    out_data_sim_0xmf.PointArrayStatus = [name_dset]

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    ## show data in view
    #out_data_sim_0xmfDisplay = Show(out_data_sim_0xmf, renderView1, 'UnstructuredGridRepresentation')

    ## trace defaults for the display properties.
    #out_data_sim_0xmfDisplay.Representation = 'Surface'

    # reset view to fit data
    #renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # get color transfer function/color map for 'T_res_src_0_inv_0000'
    name_dsetLUT = GetColorTransferFunction(name_dset)

    # get opacity transfer function/opacity map for 'T_res_src_0_inv_0000'
    name_dsetPWF = GetOpacityTransferFunction(name_dset)


    # show color bar/color legend
    #out_data_sim_0xmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    #renderView1.Update()

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

    # trace defaults for the display properties.
    slice1Display.Representation = 'Surface'


    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)


    #
    # change the position of color bar
    #
    # Properties modified on slice1Display.DataAxesGrid
    slice1Display.DataAxesGrid.GridAxesVisibility = 1

    # get color legend/bar for t_res_src_0_inv_0000LUT in view renderView1
    name_dsetLUTColorBar = GetScalarBar(name_dsetLUT, renderView1)

    # change scalar bar placement
    name_dsetLUTColorBar.Orientation = 'Horizontal'
    name_dsetLUTColorBar.WindowLocation = 'Any Location'
    name_dsetLUTColorBar.Position = [0.33262344353799916, 0.24363636363636363]
    name_dsetLUTColorBar.ScalarBarLength = 0.3299999999999998



    # update the view to ensure updated data information
    #renderView1.Update()


    #
    # axis options
    #

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=slice1.SliceType)

    # Properties modified on slice1Display.DataAxesGrid
    slice1Display.DataAxesGrid.GridAxesVisibility = 1


    #
    # window layout and camera settings
    #
    layout1 = GetLayout()

    #--------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(1280, 980)

    #-----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [6086.06543558693, 0.0, 7825.811009542308]
    renderView1.CameraFocalPoint = [6086.06543558693, 0.0, -3.410605131648481e-13]
    renderView1.CameraViewUp = [1.0, 2.220446049250313e-16, 0.0]
    renderView1.CameraParallelScale = 2025.468932642534

    # save screenshot
    SaveScreenshot('test_screenshot.png', renderView1, ImageResolution=[2560,1280])

    Render()
    Interact()