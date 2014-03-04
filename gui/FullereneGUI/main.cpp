#include <QtQml>
#include <QtDebug>
#include <QtGui/QGuiApplication>
#include "qtquick2applicationviewer.h"
#include "fullereneselect.h"

#include <libgraph/isomerdb.hh>
#include <contrib/buckygen-wrapper.hh>

int main(int argc, char *argv[])
{
    QGuiApplication app(argc, argv);

    QtQuick2ApplicationViewer viewer;
    viewer.setMainQmlFile(QStringLiteral("qml/FullereneGUI/main.qml"));
    viewer.showExpanded();


    QObject *root = (QObject *)viewer.rootObject();
    QObject *button_OK = root->findChild<QObject*>("button_OK");


 //   button_OK->setProperty("text","BING!");

//    QObject::connect(button_OK, SIGNAL(clicked()), this, SLOT(loginClick()));

    FullereneSelect object(root);
    object.setObjectName("hellofriend");
    qDebug() << "Object name is: " << object.objectName() << endl;

    QQmlContext *ctx = viewer.rootContext();
    //ctx->setContextObject(&object);
    ctx->setContextProperty(QString("hellofriend"),&object);
    return app.exec();
}
