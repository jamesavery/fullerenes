import numpy;
import sys;
from PyQt4 import QtCore, QtGui;

from fullerenegui import FullereneGUI;
from fullwrap import *;

class FullereneWidget(FullereneGUI):
    def __init__(self,parent=None):
        FullereneGUI.__init__(self,parent);
        print "Shalom;"

        self.table_select.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows);
        self.table_select.setSelectionMode(QtGui.QAbstractItemView.SingleSelection);

        self.connect(self.button_show,       QtCore.SIGNAL('clicked()'),                 self.populate_table);
        self.connect(self.spinbox_N,         QtCore.SIGNAL('editingFinished()'),         self.N_changed);
        self.connect(self.combobox_symmetry, QtCore.SIGNAL('currentIndexChanged(int)'),  self.filter_changed);
        self.connect(self.table_select,      QtCore.SIGNAL('itemSelectionChanged()'), self.selection_changed);


    def set_label_Nisomers(self):
        N   = self.spinbox_N.value();
        IPR = self.checkbox_IPR.isChecked();
        filter_symmetry = str(self.combobox_symmetry.currentText());

        Nisomers     = IsomerDB.number_isomers(N,"Any",IPR);        
        Nisomers_sym = IsomerDB.number_isomers(N,filter_symmetry,IPR);

        self.label_Nisomers.setText(str(Nisomers)+" C"+str(N)+" isomers in total, "+str(Nisomers_sym)+" matching filter.");

    def N_changed(self):
        N = self.spinbox_N.value();
        IPR = self.checkbox_IPR.isChecked();
        Nisomers = IsomerDB.number_isomers(N,"Any",IPR);
        syms     = list(IsomerDB.symmetries(N,IPR));

        self.spinbox_iso_to.setValue(Nisomers);

        dropdown = self.combobox_symmetry;
        dropdown.clear();
        dropdown.insertItem(0,"Any");

        for i in range(len(syms)):
            dropdown.insertItem(i+1,syms[i]);

        dropdown.setCurrentIndex(0);
        self.set_label_Nisomers();
        
    def filter_changed(self,sym_index):
        self.set_label_Nisomers();
 
    def selection_changed(self):
#        print [str(I.text()) for I in self.table_select.selectedItems()]
        return;

    def populate_table(self):
        N = self.spinbox_N.value();
        IPR = self.checkbox_IPR.isChecked();
        iso_from = self.spinbox_iso_from.value();
        iso_to   = self.spinbox_iso_to.value();
        sym_filter = str(self.combobox_symmetry.currentText());

        state = self.fullstate = FullereneSelect(N,IPR,"" if(sym_filter == "Any" or sym_filter == " C1") else "-nontrivial");
        if(state.has_db): 
            Nisomers = state.db.Nisomers;
        else:
            Nisomers = iso_to-iso_from+1;

        entries = self.entries = state.get_fullerenes(Nisomers,iso_from,iso_to,sym_filter);

        self.table_select.clearContents();
        self.table_select.setRowCount(len(entries));
        for row in range(len(entries)):
            e = entries[row];
            v = [iso,sym,oc,gap,nmr,spiral] = [QtGui.QTableWidgetItem() for x in range(6)];
            iso.setText(str(e.isomer));
            sym.setText(e.group);
            oc.setText(str(e.NeHOMO));
            gap.setText(str(e.HLgap));
            nmr.setText(str(list(e.INMR))[1:-1]);
            spiral.setText(str(list(e.RSPI))[1:-1]);
            for i in range(6):
                self.table_select.setItem(row,i,v[i]);

        self.table_select.resizeColumnsToContents()

    def current_isomer_geometry(self):
#        rspi = [int(x) for x in str(self.table_select.selectedItems()[5].text()).split(',')];
#        print rspi;
        isomer = int(str(self.table_select.selectedItems()[0].text()));
        positions = self.fullstate.get_coordinates(isomer,3,1e-10);
        positions = numpy.array(positions).reshape(len(positions)/3,3);
        return positions;

    def current_isomer_title(self):
        row = [x.text() for x in list(self.table_select.selectedItems())];
        print "Selected row: ", row
        N = self.fullstate.db.N;
        (isomer,group) = (int(row[0]),str(row[1]));
        return "C%d-%s, isomer number %d" % (N,group,isomer);




if __name__ == "__main__":
    print "Constructing QtGui"
    app = QtGui.QApplication(sys.argv)

    print "Constructing widget"
    widget = FullereneWidget()
    print "Showing widget"
    widget.show()
    print "Waiting"

#    sys.exit(app.exec_())

