import vtk
from numpy import asarray

def CombineGaussians(bounds, coords, weights, sizes):
    """Forms a linear combination of Gaussians"""

    assert len(bounds) == 6, 'Invalid bounds'
    assert len(coords) == len(weights) == len(sizes), 'Different numbers of coordinates and weights'
    combo = vtk.vtkImageWeightedSum()
    for i in range(len(coords)):
        primitive = vtk.vtkImageGaussianSource()
        primitive.SetWholeExtent(*bounds)
        primitive.SetCenter(coords[i])
        primitive.SetStandardDeviation(sizes[i])
        primitive.SetMaximum(weights[i]) #WRONG normalization
        combo.AddInputConnection(primitive.GetOutputPort())
        combo.SetWeight(i, 1.0)

    return combo



def pGaussian(bounds, cart, pos, size, h = 0.2):
    """Makes a p-type Cartesian Gaussian using finite difference
    cart = Cartesian direction: 0=x, 1=y, 2=z
    """
    pos = asarray(pos)

    if cart == 0:
        vec = (h, 0, 0)
    elif cart == 1:
        vec = (0, h, 0)
    elif cart == 2:
        vec = (0, 0, h)
    else:
        assert False, 'Invalid cart'
    vec = asarray(vec)

    coords = [pos - vec, pos + vec]
    weights = [1, -1]
    sizes = (size,) * 2
    return CombineGaussians(bounds, coords, weights, sizes)



def Rescale(combo, ZeroPoint):
    """Rescale to legal pixel values [-1, 1] --> [0, MaxScaledValue]"""
    scaled = vtk.vtkImageShiftScale()
    scaled.SetInputConnection(combo.GetOutputPort())
    scaled.SetShift(1.)
    scaled.SetScale(ZeroPoint + 1)
    scaled.SetOutputScalarTypeToUnsignedShort()

    return scaled



def MakeVTKColorScheme(ZeroPoint, MaxScaledVal, poscolor = (0.95, 0.5, 0.14), negcolor = (0, 0.23, 0.48)):
    #Opacity transfer function
    opacitytf = vtk.vtkPiecewiseFunction()
    opacitytf.AddPoint(0, 1.0)
    opacitytf.AddPoint(ZeroPoint - 1, 0.0)
    opacitytf.AddPoint(ZeroPoint + 1, 0.0)
    opacitytf.AddPoint(MaxScaledVal, 1.0)

    #Color transfer function
    colortf = vtk.vtkColorTransferFunction()
    colortf.AddRGBPoint(0, *negcolor)
    colortf.AddRGBPoint(ZeroPoint - 1, *negcolor)
    colortf.AddRGBPoint(ZeroPoint + 1, *poscolor)
    colortf.AddRGBPoint(MaxScaledVal, *poscolor)

    volprp = vtk.vtkVolumeProperty()
    volprp.SetColor(colortf)
    volprp.SetScalarOpacity(opacitytf)
    volprp.ShadeOn()
    volprp.SetInterpolationTypeToLinear()

    return volprp



def RenderVTKVolume(image, volprops):
    volmap = vtk.vtkVolumeRayCastMapper()
    volmap.SetVolumeRayCastFunction(vtk.vtkVolumeRayCastCompositeFunction())
    volmap.SetInputConnection(image.GetOutputPort())

    vol = vtk.vtkVolume()
    vol.SetMapper(volmap)
    vol.SetProperty(volprops)

    #Standard VTK stuff
    ren = vtk.vtkRenderer()
    ren.AddVolume(vol)
    ren.SetBackground((1, 1, 1))

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(ren)

    istyle = vtk.vtkInteractorStyleSwitch()
    istyle.SetCurrentStyleToTrackballCamera()

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.SetInteractorStyle(istyle)

    renwin.Render()
    iren.Start()


def main():
    bounds = (-5, 5, -5, 5, -5, 5)
    coords = [(-1, 0, 0), (1, 0, 0)]
    weights = (1, -1)
    sizes = (1, 1)

    MaxScaledVal = 65536
    ZeroPoint = int(MaxScaledVal / 2)

    #combo = CombineGaussians(bounds, coords, weights, sizes, ZeroPoint)
    combo = pGaussian(bounds, 0, (0, 0, 0), 1)
    volobj = Rescale(combo, ZeroPoint)
    volprp = MakeVTKColorScheme(ZeroPoint, MaxScaledVal)
    RenderVTKVolume(volobj, volprp)


if __name__ == '__main__':
    main()

