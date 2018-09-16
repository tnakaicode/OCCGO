import os

from OCC.Graphic3d import Graphic3d_NOM_ALUMINIUM
from OCC.TopoDS import TopoDS_Shape
from OCC.StlAPI import StlAPI_Writer
from OCC.AIS import AIS_Shape
from OCC.gp import gp_Pnt
from OCCUtils.Topology import dumpTopology
import airconics as act

from .base import Base
from .export import export_STEPFile

def print_xy_click(shp, *kwargs):
    shape = shp[0]
    print("Shape selected: %s" % (shape))
    dumpTopology (shape)
    print ("dumpTopology COMPLETE")
                                                                        
class Part (Base):

    def __init__(self, components={}, *args, **kwargs):
        # Set the components dictionary (default empty)
        self._Components = {}

        for name, component in components.items():
            self.__setitem__(name, component)

        for key, value in kwargs.items():
            self.__setattr__(key, value)

    def __getitem__(self, name):
        return self._Components[name]

    def __setitem__(self, name, component):
        self._Components[name] = component

    def __delitem__(self, name):
        del self._Components[name]

    def __iter__(self):
        return iter(self._Components)

    def __len__(self):
        return len(self._Components)

    def __str__(self):
        output = str(self.keys())
        return output

    def AddComponent(self, component, name=None):
        if name is None:
            name = 'untitled_component_' + str(len(self)+1)
        self.__setitem__(name, component)
    
    def AddComponents(self, components, name=None):
        if name is None:
            name = 'untitled_components_' + str(len(self)+1)
        for i, comp in enumerate(components):
            self.__setitem__(name + '-' + str(i+1), comp)
    
    def RemoveComponent(self, name):
        self.__delitem__(name)
    
    def Build (self):
        print ("ok")
    
    def Display (self, context, material=Graphic3d_NOM_ALUMINIUM, color=None):
        print ("display-Part")
        for name, component in self.items():
            ais = AIS_Shape(component)
            ais.SetMaterial(material)
            if color:
                try:
                    from OCC.Display.OCCViewer import get_color_from_name
                    color = get_color_from_name(color)
                except:
                    pass
                ais.SetColor(color)
            try:
                context.Context.Display(ais.GetHandle())
                #context.register_select_callback(print_xy_click)
            except:
                context.DisplayShape(component)
                #context.register_select_callback(print_xy_click)
    
    def Write (self, filename, single_export=True):
        print ("write-Part")
        path, ext = os.path.splitext(filename)

        # Default to a step writer if no extension type was provided:
        if not ext:
            ext = '.stp'

        status = []
        if ext == '.stl':
            stl_ascii_format = False
            if single_export:
                stl_writer = StlAPI_Writer()
                for name, component in self.items():
                    shape = component
                    status.append(stl_writer.Write(shape, filename, stl_ascii_format))
            else:
                for name, component in self.items():
                    stl_writer = StlAPI_Writer()
                    f = path + '_' + name + ext
                    shape = component
                    status.append(stl_writer.Write(shape, f, stl_ascii_format))

        elif ext in ['.stp', '.step']:
            if single_export:
                status.append(export_STEPFile(list(self.values()), filename))
            else:
                for name, component in self.items():
                    f = path + '_' + name + ext
                    status.append(export_STEPFile([component], f))
        else:
            raise ValueError('Unexpected file extension {}'.format(ext))

        return status

class Product (Part):

    def __init__(self, parts={}, *args, **kwargs):
        # Set the components dictionary (default empty)
        self._Parts = {}
        for name, part in parts.items():
            self.__setitem__(name, part)

        # Set all kwargs as attributes
        for key, value in kwargs.items():
            self.__setattr__(key, value)

    def __getitem__(self, name):
        return self._Parts[name]

    def __setitem__(self, name, part):
        self._Parts[name] = part

    def __delitem__(self, name):
        del self._Parts[name]

    def __iter__(self):
        return iter(self._Parts)

    def __len__(self):
        return len(self._Parts)

    def __str__(self):
        output = str(self.keys())
        return output

    def Build (self):
        print ("ok")

    def AddPart (self, part, name=None):
        if name is None:
            # set a default name:
            name = 'untitled_part_' + str(len(self))
        self.__setitem__(name, part)

    def AddProduct (self, product, name=None):
        if name is None:
            # set a default name:
            name = 'untitled_product_' + str(len(self))
        self.__setitem__(name, product)

    def Write(self, filename, single_export=True):
        print ("write-Product")
        path, ext = os.path.splitext(filename)

        status = []

        if ext == '.stl':
            stl_ascii_format = False

            if single_export:
                stl_writer = StlAPI_Writer()
                for partname, part in self.items():
                    for name, component in part.items():
                        status.append(stl_writer.Write(component, filename,
                                                       stl_ascii_format))
            else:
                for partname, part in self.items():
                    f = path + '_' + name + ext
                    status.append(path.Write(f, single_export=True))

        elif ext in ['.stp', '.step']:
            if single_export:
                shapes = []
                for partname, part in self.items():
                    print (part)
                    shapes.append (part)
                status.append(export_STEPFile(shapes, filename))
            else:
                for name, part in self.items():
                    f = path + '_' + name + ext
                    # Writes one file per part
                    p = Part ()
                    p.AddComponent(part, name)
                    status.append(p.Write(f, single_export=True))
        
        """
        Aim Output STEP file (single_export==True)
        Product1 -- Part1 -- Comp1
                  |        - Comp2
                  - Part2 -- Comp1
                  |        - Comp2
                  - Comp1
        """

        return status

    def Display(self, context, material=Graphic3d_NOM_ALUMINIUM, color=None):
        print ("display-Product")
        for name, component in self.items():
            try:
                component.Display(context, material, color)
            except AttributeError:
                # We are probably dealing with a core pythonocc shape:
                try:
                    context.DisplayShape(component)
                except:
                    print("Could not display shape type {}: skipping".format(
                        type(component)))

class Assembly (Base):

    def __init__(self, components={}, *args, **kwargs):
        # Set the components dictionary (default empty)
        self._Components = {}

        for name, component in components.items():
            self.__setitem__(name, component)

        for key, value in kwargs.items():
            self.__setattr__(key, value)
    
    def __str__(self):
        print ("str")

    def Display(self):
        print ("display")

    def Write(self):
        print ("write")