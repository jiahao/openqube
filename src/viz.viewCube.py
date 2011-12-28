#!/usr/bin/env python

############
# Parameters
############

DEBUG = True

#Specify 16-bit or 8-bit color depth
ColorDepth = 16

#Dimensions of output window (in pixels)
OutputHeight = 1200
OutputWidth  = 1000

#############
#Draws volume

DrawVolume = True

#Specify color scheme
#RGB triples with each value ranging from 0.0 - 1.0
BackgroundColor = (1.0, 1.0, 1.0)

##UIUC official colors
#NegativeColor = (0.96, 0.50, 0.14)
#PositiveColor = (0.00, 0.24, 0.49)

#MIT official colors
NegativeColor = (0.4, 0.4, 0.4)
PositiveColor = (0.6, 0.2, 0.2)

########################
#Draws in nodal surfaces

DrawNodes = True

NodeColor = (0.0, 1.0, 1.0)
NodeAlpha = 0.3

################
#Draws wireframe

DrawFrame = True

FrameColor = (0.9, 0.9, 0.9)
FrameAlpha = 1.0

############
#Draws atoms

DrawAtoms = True

############
#Draws bonds

DrawBonds = True

Interactive = True

########################################################

##################
# Code starts here
##################

#Set up VTK
import vtk

def RenderCubeInVTK(filename = 'test.cube', mindatum = 0.0, maxdatum = 0.0):
    global renWin

    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    #######################
    # Read in Gaussian cube
    #######################

    CubeData = vtk.vtkGaussianCubeReader()
    CubeData.SetFileName(filename)

    CubeData.Update()

    #Get intrinsic scale from data
    scale = sum([x**2 for x in CubeData.GetTransform().GetScale()])

    CubeData.SetHBScale(scale) #scaling factor to compute bonds with hydrogens
    CubeData.SetBScale(scale)  #scaling factor for other bonds

    CubeData.Update()

    ###################
    #Calculate scalings

    #VTK only knows how to render integer data in the interval [0,255] or [0,65535]
    #Here, we calculate scaling factors to map the cube data to the interval.

    if mindatum == maxdatum == 0.0:
        if DEBUG:
            print "Autodetecting range"
            mindatum, maxdatum = CubeData.GetGridOutput().GetPointData().GetScalars().GetRange()

    # Find the remapped value that corresponds to zero
    zeropoint = int(2**ColorDepth*(-mindatum)/(maxdatum-mindatum))
    absmaxdatum = max(-mindatum, maxdatum)

    maxnegativeintensity = min(1.0, 1.0 - (absmaxdatum - abs(mindatum))/absmaxdatum)
    minnegativeintensity = 0.0
    if zeropoint < 0:
        minpositiveintensity = - zeropoint/(2**ColorDepth*absmaxdatum)
    else:
        minpositiveintensity = 0.0
        maxpositiveintensity = min(1.0, 1.0 - (absmaxdatum - abs(maxdatum))/absmaxdatum)
    if DEBUG:
        print "Range plotted = [%f,%f]" % (mindatum, maxdatum)
        print "Negative colors = [0,%d)" % max(0,zeropoint)
        print "Negative intensities = [%f,%f]" % (maxnegativeintensity,minnegativeintensity)
        print "Positive colors = (%d,%d)" % (max(0,zeropoint), 2**ColorDepth)
        print "Positive intensities = [%f,%f]" % (minpositiveintensity,maxpositiveintensity)
        print "On this scale, zero = %d" % zeropoint

    ################################
    # Calculate opacity transfer map

    #The code here differentiates between two cases:
    #1. the scalar data are all positive, so it's just a simple linear ramp
    #2. the scalar data are signed, so do two linear ramps

    opacityTransferFunction = vtk.vtkPiecewiseFunction()

    if zeropoint < 0:
        opacityTransferFunction.AddPoint(        0, minpositiveintensity)
    else:
        opacityTransferFunction.AddPoint(        0, maxnegativeintensity)
        opacityTransferFunction.AddPoint(zeropoint, 0.0)

    opacityTransferFunction.AddPoint(2**ColorDepth-1, maxpositiveintensity)
    opacityTransferFunction.ClampingOn()

    ###########################
    # Create color transfer map

    colorTransferFunction = vtk.vtkColorTransferFunction()

    r1, g1, b1 = NegativeColor
    r2, g2, b2 = PositiveColor
    r0, g0, b0 = BackgroundColor

    if zeropoint < 0:
        colorTransferFunction.AddRGBPoint(          0, r1, g1, b1)
    else:
        colorTransferFunction.AddRGBPoint(          0, r1, g1, b1)
        colorTransferFunction.AddRGBPoint(zeropoint-1, r1, g1, b1)
        colorTransferFunction.AddRGBPoint(zeropoint  , r0, g0, b0)
        colorTransferFunction.AddRGBPoint(zeropoint+1, r2, g2, b2)
    colorTransferFunction.AddRGBPoint(2**ColorDepth-1, r2, g2, b2)

    ########################
    # Now apply the scalings

    ScaledData = vtk.vtkImageShiftScale()
    ScaledData.SetInput(CubeData.GetGridOutput())
    ScaledData.SetShift(-mindatum)
    ScaledData.SetScale((2**ColorDepth-1)/(maxdatum-mindatum))

    if ColorDepth == 16:
        ScaledData.SetOutputScalarTypeToUnsignedShort()
    elif ColorDepth == 8:
        ScaledData.SetOutputScalarTypeToUnsignedChar()
    else:
        print
        print "Error! Unsupported color depth given"
        print
        print "valid values are 8 or 16"
        print
        raise ValueError

    ###############################
    # Form combined coloring scheme

    volumeProperty = vtk.vtkVolumeProperty()
    volumeProperty.SetColor(colorTransferFunction)
    volumeProperty.SetScalarOpacity(opacityTransferFunction)
    volumeProperty.SetInterpolationTypeToLinear()
    volumeProperty.ShadeOn()

    # The mapper / ray cast function know how to render the data
    compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()

    volumeMapper = vtk.vtkVolumeRayCastMapper()
    volumeMapper.SetVolumeRayCastFunction(compositeFunction)
    volumeMapper.SetInput(ScaledData.GetOutput())

    #Create a coarse representation
    #Actually a fake - won't display anything

    compositeFunction2 = vtk.vtkVolumeRayCastIsosurfaceFunction()
    compositeFunction2.SetIsoValue(2**ColorDepth-1)

    volumeMapperCoarse = vtk.vtkVolumeRayCastMapper()
    volumeMapperCoarse.SetVolumeRayCastFunction(compositeFunction2)
    volumeMapperCoarse.SetInput(ScaledData.GetOutput())

    # Create volumetric object to be rendered
    # Use level of detail prop so that it won't take forever to look around

    volume = vtk.vtkLODProp3D()
    id1 = volume.AddLOD(volumeMapper, volumeProperty, 0.)
    volume.SetLODProperty(id1, volumeProperty)
    id2 = volume.AddLOD(volumeMapperCoarse, volumeProperty, 0.)
    volume.SetLODProperty(id2, volumeProperty)

    # At this point, we can position and orient the volume

    #################################
    # End of volumetric data pipeline
    #################################

    #########
    #Contours
    #########

    contour = vtk.vtkContourFilter()
    contour.SetInput(CubeData.GetGridOutput())
    contour.SetNumberOfContours(1)
    contour.SetValue(0, 0.0)

    contourMapper = vtk.vtkPolyDataMapper()
    contourMapper.SetInput(contour.GetOutput())
    contourMapper.SetScalarRange(0,0)
    contourMapper.GetLookupTable().SetNumberOfTableValues(1)
    r0, g0, b0 = NodeColor
    contourMapper.GetLookupTable().SetTableValue(0, r0, g0, b0, NodeAlpha)

    contourActor = vtk.vtkLODActor()
    contourActor.SetMapper(contourMapper)
    contourActor.GetProperty().SetOpacity(NodeAlpha)

    ##########################################
    # Create a wireframe outline of the volume
    ##########################################

    frame = vtk.vtkOutlineFilter()
    frame.SetInput(CubeData.GetGridOutput())

    frameMapper = vtk.vtkPolyDataMapper()
    frameMapper.SetInput(frame.GetOutput())

    frameActor = vtk.vtkLODActor()
    frameActor.SetMapper(frameMapper)
    frameActor.GetProperty().SetColor(FrameColor)
    frameActor.GetProperty().SetOpacity(FrameAlpha)

    ######################
    # Draw balls for atoms
    ######################

    Sphere = vtk.vtkSphereSource()
    Sphere.SetThetaResolution(16)
    Sphere.SetPhiResolution(16)
    Sphere.SetRadius(0.4)

    Glyph = vtk.vtkGlyph3D()
    Glyph.SetInput(CubeData.GetOutput())
    Glyph.SetColorMode(1)
    Glyph.SetColorModeToColorByScalar()
    Glyph.SetScaleModeToScaleByVectorComponents()
    Glyph.SetSource(Sphere.GetOutput())

    AtomsMapper = vtk.vtkPolyDataMapper()
    AtomsMapper.SetInput(Glyph.GetOutput())
    AtomsMapper.SetImmediateModeRendering(1)
    AtomsMapper.UseLookupTableScalarRangeOff()
    AtomsMapper.SetScalarVisibility(1)
    AtomsMapper.SetScalarModeToDefault()

    Atoms = vtk.vtkLODActor()
    Atoms.SetMapper(AtomsMapper)

    ############
    # Draw bonds
    ############

    Tube = vtk.vtkTubeFilter()
    Tube.SetInput(CubeData.GetOutput())

    BondsMapper = vtk.vtkPolyDataMapper()
    BondsMapper.SetInput(Tube.GetOutput())

    Bonds = vtk.vtkLODActor()
    Bonds.SetMapper(BondsMapper)

    #######################
    # Now compose the image
    #######################
    if DrawVolume:
        ren.AddVolume(volume)

    if DrawNodes:
        ren.AddActor(contourActor)

    if DrawFrame:
        ren.AddActor(frameActor)

    if DrawAtoms:
        ren.AddActor(Atoms)

    if DrawBonds:
        ren.AddActor(Bonds)

    ren.SetBackground(BackgroundColor)
    renWin.SetSize(OutputHeight, OutputWidth)


    ######################################
    # Let VTK do its magic and render away
    ######################################

    renWin.Render()

    ###################################
    # Now allow user to play with image
    ###################################

    def Keypress(obj, event):
        #This function handles keyboard interaction

        key = obj.GetKeySym()

        if key == 'd' or key == 'F13':
            WriteToPNG()
        elif key == 'h' or key == 'question' or key =='?':
            PrintHelp()
        elif key == 'c':
            camera = ren.GetActiveCamera()
            print "Camera info:"
            print "------------"
            print "Position is: ", camera.GetPosition()
            print "Focal point is:", camera.GetFocalPoint()
            print "Orientation is:", ren.GetActiveCamera().GetOrientation()
            print "WXYZ", ren.GetActiveCamera().GetOrientationWXYZ()
            print "View up direction is:", camera.GetViewUp()
            print "Direction of projection is:", camera.GetDirectionOfProjection()

        else:
            if DEBUG:
                print 'User pressed key:', key 

    if Interactive:
        iren.SetDesiredUpdateRate(25.0) #25 fps when camera is moving around
        iren.SetStillUpdateRate(0.0) #0 fps when camera is not moving

        iren.Initialize()

        #The default interaction style is joystick, which seems unnatural
        style = vtk.vtkInteractorStyleTrackballCamera()
        iren.SetInteractorStyle(style)

        iren.AddObserver("KeyPressEvent", Keypress)
        iren.Start()
    else:
        WriteToPNG()


def PrintHelp():
    print """Help:

Mouse actions:

Button 1: rotates view
Button 2: pans (translates) view (Note: with 2-button mice, pan is defined as <Shift>-Button 1.)
Button 3: zoom view

Keys

c: prints information about camera position

d, PrtScn: dumps the current rendered window to screenshot.png

?, h: print this help dialog

---

3: toggle the render window into and out of stereo mode.

a: toggle between camera and actor modes. In camera mode, mouse events affect the camera position and focal point. In actor mode, mouse events affect the actor that is under the mouse pointer.

j, t: toggle between joystick (position sensitive) and trackball (motion sensitive) styles. In joystick style, motion occurs continuously as long as a mouse button is pressed. In trackball style, motion occurs when the mouse button is pressed and the mouse pointer moves.

e, q: exit the application.

f: fly to the picked point

p: perform a pick operation. The render window interactor has an internal instance of vtkCellPicker that it uses to pick.

r: reset the camera view along the current view direction. Centers the actors and moves the camera so that all actors are visible.

s: modify the representation of all actors so that they are surfaces. (Default)

w: modify the representation of all actors so that they are wireframe.
    """



def WriteToPNG(filename = 'screenshot.png'):
    global renWin
    winImg = vtk.vtkWindowToImageFilter()
    winImg.SetInput(renWin)
    winImg.SetInputBufferTypeToRGBA()
    winImg.Update()

    #Exports
    screenshot= vtk.vtkPNGWriter()
    screenshot.SetInputConnection(winImg.GetOutputPort())
    screenshot.SetFilePrefix(filename)
    screenshot.Write()

    #The following line fixes annoying behavior with vtkPNGWriter
    import os
    os.rename(filename+'.0',filename)
    if DEBUG:
        print "Wrote screenshot to ",filename 



def ScanCubeForRange(filename = 'test.cube'):
    CubeData = vtk.vtkGaussianCubeReader()
    CubeData.SetFileName(filename)
    CubeData.Update()
    return CubeData.GetGridOutput().GetPointData().GetScalars().GetRange()



if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        cubefilename = sys.argv[1]
        RenderCubeInVTK(cubefilename)

    else:
        #Nothing specified, goes into auto-scanning mode

        print """I will now go ahead and render all Gaussian cube files in the
current directory to PNG files. To avoid this behavior, specify the following command-line
options:"""
        print "Usage: "+sys.argv[0]+" [Gaussian Cube file name] [minimum (optional)] [maximum (optional)] [VTK data file name (optional)]"
        import glob

        print "I will now scan all cube files for global minimum and maximum"
        globalmax = -9999999.9
        globalmin =  9999999.9
        for cubefilename in glob.glob("*.cube"):
            print cubefilename
            thismin, thismax = ScanCubeForRange (cubefilename)
            globalmax = max(globalmax, thismax)
            globalmin = min(globalmin, thismin)
            print "The global range is [%f,%f]" % (globalmin, globalmax)

        print "I will now go ahead and render all cube files using the global range"

        Interactive = False

        for cubefilename in glob.glob("*.cube"):
            RenderCubeInVTK(cubefilename, globalmin, globalmax)


