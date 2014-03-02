import QtQuick 2.0
import QtQuick.Controls 1.1


Rectangle {
    id: select_body
    width: 640
    height: 360

    Column {
        id: column1
        width: parent.width
        spacing: 20
        Row {
            id: row3
            spacing: 10

            Row {
                id: row2
                spacing: 2
                Label {
                    id: label_N
                    text: qsTr("N: ")
                    anchors.verticalCenter: parent.verticalCenter
                }
                SpinBox {
                    id: spinbox_N
                    minimumValue: 20
                    stepSize: 2
                }
            }

            CheckBox {
                id: checkBox1
                text: qsTr("IPR")
                anchors.verticalCenter: parent.verticalCenter
            }

            Row {
                id: row1
                spacing: 2
                Label {
                    id: label1
                    text: qsTr("Symmetry:")
                    anchors.verticalCenter: parent.verticalCenter
                }
                ComboBox {
                    id: combobox_symmetry
                    activeFocusOnPress: true
                    model: ["Any","C1","Cs","Ci","C2","C2v","C3","C3v","etc."]
                    onCurrentIndexChanged: console.debug(currentIndex)
                }
            }


            Row {
                spacing: 2

                Label { text:"Isomers"; anchors.verticalCenter: parent.verticalCenter}

                SpinBox {
                    id: spinbox_iso_min
                    minimumValue: 0
                    maximumValue: spinbox_iso_max.value
                }
                SpinBox {
                    id: spinbox_iso_max
                    minimumValue: 0
                    maximumValue: 9999
                }

            }

        }

        Row {
            width: parent.width
            TableView {
                width: parent.width
                TableViewColumn { role: "iso"; title: "Isomer no."    }
                TableViewColumn { role: "sym"; title: "Symmetry"      }
                TableViewColumn { role: "oc";  title: "Valence shell" }
                TableViewColumn { role: "gap"; title: "HOMO-LUMO gap" }
                TableViewColumn { role: "mnr"; title: "NMR pattern"   }
            }
        }
   }
}
