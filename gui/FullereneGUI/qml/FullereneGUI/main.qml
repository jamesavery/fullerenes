import QtQuick 2.0
import QtQuick.Controls 1.1

Rectangle {
    id: rectangle2
    width: 640
    height: 480

    Rectangle {
        id: rect_header
        height: 54
        color: "#ffffff"
        anchors.left: parent.left
        anchors.leftMargin: 0
        anchors.right: parent.right
        anchors.rightMargin: 0
        anchors.top: parent.top
        anchors.topMargin: 0

        Text {
            id: text_header
            x: 95
            y: 63
            color: "#172454"
            text: qsTr("Fullerene Test")
            font.bold: true
            horizontalAlignment: Text.AlignHCenter
            anchors.horizontalCenter: parent.horizontalCenter
            anchors.verticalCenter: parent.verticalCenter
            font.pixelSize: 24
        }
    }

    Rectangle {
        id: body
        y: 60
        height: 358
        color: "#ffffff"
        anchors.bottom: parent.bottom
        anchors.bottomMargin: 70
        anchors.left: parent.left
        anchors.leftMargin: 0
        anchors.right: parent.right
        anchors.rightMargin: 0

        TabView {
            height: 360
            anchors.fill: parent
            Component.onCompleted: {
                addTab("Select Fullerene", Qt.createComponent("fullerene-select.qml"));
                addTab("Transforms", tab_transform);
                addTab("2D Graph Layout",tab_layout2d)
                addTab("3D Geometry",tab_geometry)
                addTab("Analysis", tab_analysis)
                addTab("Output", tab_output)
            }

            Component {
                id: tab_transform
                Rectangle { color:"blue" }
            }
            Component {
                id: tab_analysis
                Rectangle { color: "purple" }
            }
            Component {
                id: tab_output
                Rectangle { color: "green" }
            }
            Component {
                id: tab_geometry
                Rectangle { }

            }
            Component {
                id: tab_layout2d
                Rectangle {}
            }
        }
    }
}
