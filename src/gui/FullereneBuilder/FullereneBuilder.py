import numpy
from PyQt4 import QtCore, QtGui
from NL.CommonConcepts import PeriodicTable
from NL.CommonConcepts import PhysicalQuantity as Units
from NL.CommonConcepts.Configurations.BravaisLattice import UnitCell
from NL.CommonConcepts.Configurations.BulkConfiguration import BulkConfiguration
from NL.CommonConcepts.Configurations.MoleculeConfiguration import MoleculeConfiguration
from NL.CommonConcepts.PeriodicTable import Carbon
from NL.CommonConcepts.Configurations.NanoTube import NanoTube
from NL.CommonConcepts.ElementProperties import ElementProperties
from NL.CommonConcepts.PhysicalQuantity import Angstrom
from NL.GUI.Presenter.Presenter import Presenter
from NL.GUI.ScienceWidgets.PeriodicTable import PeriodicTable as PeriodicTableWidget
from NL.GUI.Tools.Builder.Stash.StashItem import StashItem
from NL.GUI.UI.CustomWidgets.Widgets.FloatLineEdit import FloatLineEdit


from FullereneWidget import FullereneWidget
from API import BuilderStashAddPlugin

class FullereneBuilder(BuilderStashAddPlugin):
    """ Class for building Carbon Nanotube builder """

    def setupWidget(self):
        """
        Create the widget.
        """
        if self._widget is None:
            self._widget = FullereneWidget(parent=self.parent())
            self._widget.setWindowFlags(QtCore.Qt.Tool)

            # Setup the initial element dialog.
            self._element_dialog = PeriodicTableWidget(multi_selectable=False)
            self._element_dialog.setWindowTitle("Select element")
            self._element_dialog.setWindowIcon(QtGui.QIcon(':/vnl_icon.png'))

            # Disable the spin boxes.
#            self.checkState()

            # Make sure preview is removed when the plugin is closed
            self._widget.closeEvent = self.close
            self._widget.connect(self._widget.button_preview, QtCore.SIGNAL('clicked()'), self.clicked_preview)
            self._widget.connect(self._widget.button_build, QtCore.SIGNAL('clicked()'), self.clicked_build)
        return self._widget


    def clicked_preview(self):
        positions = self._widget.current_isomer_geometry()*Angstrom;
        molecule = MoleculeConfiguration(elements=[Carbon for X in positions],cartesian_coordinates=positions)
        self.stash().setPreview(Presenter(molecule),self)

    def clicked_build(self):
        positions = self._widget.current_isomer_geometry()*Angstrom;
        molecule = MoleculeConfiguration(elements=[Carbon for X in positions],cartesian_coordinates=positions);
        S = self.stash()
        S.addConfiguration(molecule,set_active=True)
        print self._widget.current_isomer_title();
        S.activeItem().setTitle(self._widget.current_isomer_title())

    def close(self, event):
        """ Called when window is being closed - remove the preview """
        self.stash().setPreview(None)
        event.accept()

    def initialize(self):
        """
        Initialize the CNT Builder.
        """
        self._widget = None



    # def checkState(self):
    #     """
    #     Set the check state of the spin boxes.
    #     """
    #     # Construct a widget if needed.
    #     widget = self.setupWidget();
    #     # Disable the spin boxes if needed
    #     widget._wall1_n.setEnabled(widget._wall1.isChecked())
    #     widget._wall1_m.setEnabled(widget._wall1.isChecked())
    #     widget._wall2_n.setEnabled(widget._wall2.isChecked())
    #     widget._wall2_m.setEnabled(widget._wall2.isChecked())
    #     widget._wall3_n.setEnabled(widget._wall3.isChecked())
    #     widget._wall3_m.setEnabled(widget._wall3.isChecked())
    #     # Toggle the buttons.
    #     widget._add_button.setEnabled(widget._wall1.isChecked() or
    #                                   widget._wall2.isChecked() or
    #                                   widget._wall3.isChecked())
    #     widget._preview_button.setEnabled(widget._wall1.isChecked() or
    #                                       widget._wall2.isChecked() or
    #                                       widget._wall3.isChecked())

    # def _elementButtonClickedCallback(self):
    #     """
    #     Callback connected to the two element buttons
    #     """
    #     button = self._widget.sender()

    #     if button is self._widget._element1:
    #         self.changeElement(1)
    #     elif button is self._widget._element2:
    #         self.changeElement(2)
    #     else:
    #         pass

    # def changeElement(self, button_id):
    #     """
    #     Change the element in one of the dialogs.
    #     """
    #     # Make sure that we have setup the widget.
    #     self.setupWidget()

    #     # Prior to running the dialog removed all checked elements.
    #     self._element_dialog.deselectAll()

    #     # Run the dialog.
    #     if self._element_dialog.exec_():
    #         # Fetch the element class.
    #         element = self._element_dialog.selectedElements()[0]
    #         # Update the button.
    #         if button_id == 1:
    #             self._widget._element1.setText(ElementProperties[element]['name'])
    #         else:
    #             self._widget._element2.setText(ElementProperties[element]['name'])

    def menuText(self):
        """
        Set the menu text.
        """
        return "Fullerene Builder"

    def use(self):
        """
        Called when the menue item is clicked.
        """
        # Construct the widget if needed.
        widget = self.setupWidget()
        # Show and activate the widget.
        widget.show()
        widget.activateWindow()

    # def setupStructure(self):
    #     """
    #     Build the structure.
    #     """
    #     # # Construct the widget if needed.
    #     # widget = self.setupWidget()

    #     # # Find the elements.
    #     # self._element1 = [p for p in PeriodicTable.allElements() if p.name() == widget._element1.text()][0]
    #     # self._element2 = [p for p in PeriodicTable.allElements() if p.name() == widget._element2.text()][0]

    #     # # Build the first wall.
    #     # if (widget._wall1.isChecked()):
    #     #     n = max(widget._wall1_n.value(), widget._wall1_m.value())
    #     #     m = min(widget._wall1_n.value(), widget._wall1_m.value())
    #     #     wall1 = NanoTube(n=n,m=m,
    #     #                      atom_1=self._element1, atom_2=self._element2,
    #     #                      bond_length=widget._bond_distance.value()*Angstrom)

    #     # # Build the second wall.
    #     # if (widget._wall2.isChecked()):
    #     #     n = max(widget._wall2_n.value(), widget._wall2_m.value())
    #     #     m = min(widget._wall2_n.value(), widget._wall2_m.value())
    #     #     wall2 = NanoTube(n=n,m=m,
    #     #                      atom_1=self._element1, atom_2=self._element2,
    #     #                      bond_length=widget._bond_distance.value()*Angstrom)

    #     # # Build the third wall.
    #     # if (widget._wall3.isChecked()):
    #     #     n = max(widget._wall3_n.value(), widget._wall3_m.value())
    #     #     m = min(widget._wall3_n.value(), widget._wall3_m.value())
    #     #     wall3 = NanoTube(n=n,m=m,
    #     #                      atom_1=self._element1, atom_2=self._element2,
    #     #                      bond_length=widget._bond_distance.value()*Angstrom)

    #     # # First find the largest cell in the X,Y axises.
    #     # L = 0.0
    #     # if (widget._wall1.isChecked()):
    #     #     L = max(float(wall1.bravaisLattice().primitiveVectors()[0][0].inUnitsOf(Angstrom)),L)
    #     # if (widget._wall2.isChecked()):
    #     #     L = max(float(wall2.bravaisLattice().primitiveVectors()[0][0].inUnitsOf(Angstrom)),L)
    #     # if (widget._wall3.isChecked()):
    #     #     L = max(float(wall3.bravaisLattice().primitiveVectors()[0][0].inUnitsOf(Angstrom)),L)

    #     # # Calculate the longest
    #     # Z = [0,0,0]

    #     # # Center the first wall.
    #     # if (widget._wall1.isChecked()):
    #     #     Z[0] = float(wall1.bravaisLattice().primitiveVectors()[2][2].inUnitsOf(Angstrom))
    #     #     wall1 = BulkConfiguration(UnitCell([L,0,0]*Angstrom, [0, L, 0]*Angstrom, [0,0,Z[0]]*Angstrom),
    #     #                               elements=wall1.elements(),
    #     #                               cartesian_coordinates=wall1.cartesianCoordinates())
    #     #     wall1 = wall1.center(True, True, False)

    #     # # Center the second wall.
    #     # if (widget._wall2.isChecked()):
    #     #     Z[1] = float(wall2.bravaisLattice().primitiveVectors()[2][2].inUnitsOf(Angstrom))
    #     #     wall2 = BulkConfiguration(UnitCell([L,0,0]*Angstrom, [0, L, 0]*Angstrom, [0,0,Z[1]]*Angstrom),
    #     #                               elements=wall2.elements(),
    #     #                               cartesian_coordinates=wall2.cartesianCoordinates())
    #     #     wall2 = wall2.center(True, True, False)

    #     # # Center the third wall.
    #     # if (widget._wall3.isChecked()):
    #     #     Z[2] = float(wall3.bravaisLattice().primitiveVectors()[2][2].inUnitsOf(Angstrom))
    #     #     wall3 = BulkConfiguration(UnitCell([L,0,0]*Angstrom, [0, L, 0]*Angstrom, [0,0,Z[2]]*Angstrom),
    #     #                               elements=wall3.elements(),
    #     #                               cartesian_coordinates=wall3.cartesianCoordinates())
    #     #     wall3 = wall3.center(True, True, False)


    #     # # Determine if any of the systems could gain from being repeated,
    #     # # before the merge.
    #     # if (widget._wall1.isChecked()):
    #     #     wall1 = wall1.repeat(1,1,int(numpy.round(max(Z)/Z[0])))
    #     # if (widget._wall2.isChecked()):
    #     #     wall2 = wall2.repeat(1,1,int(numpy.round(max(Z)/Z[1])))
    #     # if (widget._wall3.isChecked()):
    #     #     wall3 = wall3.repeat(1,1,int(numpy.round(max(Z)/Z[2])))


    #     # # Calculate the new Z-max.
    #     # Zmax = 0.0
    #     # if widget._wall1.isChecked():
    #     #     Zmax = max(float(wall1.bravaisLattice().primitiveVectors()[2][2].inUnitsOf(Angstrom)), Zmax)
    #     # if widget._wall2.isChecked():
    #     #     Zmax = max(float(wall2.bravaisLattice().primitiveVectors()[2][2].inUnitsOf(Angstrom)), Zmax)
    #     # if widget._wall3.isChecked():
    #     #     Zmax = max(float(wall3.bravaisLattice().primitiveVectors()[2][2].inUnitsOf(Angstrom)), Zmax)

    #     # # Build the values need for the super structure.
    #     # elements = []
    #     # fractional_coordinates = []
    #     # unitcell = UnitCell([L,0,0]*Angstrom, [0,L,0]*Angstrom, [0,0,Zmax]*Angstrom)
    #     # if widget._wall1.isChecked():
    #     #     elements += numpy.array(wall1.elements()).tolist()
    #     #     fractional_coordinates += numpy.array(wall1.fractionalCoordinates()).tolist()
    #     # if widget._wall2.isChecked():
    #     #     elements += numpy.array(wall2.elements()).tolist()
    #     #     fractional_coordinates += numpy.array(wall2.fractionalCoordinates()).tolist()
    #     # if widget._wall3.isChecked():
    #     #     elements += numpy.array(wall3.elements()).tolist()
    #     #     fractional_coordinates += numpy.array(wall3.fractionalCoordinates()).tolist()

    #     # # Build the super structure.
    #     # return BulkConfiguration(unitcell, elements=elements,
    #     #                          fractional_coordinates=fractional_coordinates)

    # def build(self):
    #     """
    #     Build the structure.
    #     """
    #     # Get the widget.
    #     widget = self.setupWidget()
    #     # Build the super structure.
    #     wall = self.setupStructure()
    #     # Added as the active stash.
    #     self.stash().addConfiguration(wall, set_active=True)

    #     title = ""
    #     if widget._wall1.isChecked():
    #         title += "(%d,%d) "%(widget._wall1_n.value(), widget._wall1_m.value())
    #     if widget._wall2.isChecked():
    #         title += "(%d,%d) "%(widget._wall2_n.value(), widget._wall2_m.value())
    #     if widget._wall3.isChecked():
    #         title += "(%d,%d) "%(widget._wall3_n.value(), widget._wall3_m.value())

    #     # Handle the special case of both carbon separately.
    #     if self._element1 in [PeriodicTable.Carbon] and self._element2 in [PeriodicTable.Carbon]:
    #         title += "Carbon"
    #     else:
    #         title += self._element1.symbol() + "-" + self._element2.symbol()
    #     title = title + " Nanotube"

    #     # Set the title.
    #     itemtitle = title
    #     self.stash().activeItem().setTitle(itemtitle)
    #     # Hide the widget.
    #     widget.hide()

    # def showPreview(self):
    #     presenter = Presenter(self.setupStructure())
    #     self.stash().setPreview(presenter, self)

    # def preview(self):
    #     """
    #     Build the structure.
    #     """
    #     # Build the super structure.
    #     wall = self.setupStructure()
    #     return Presenter(wall)



