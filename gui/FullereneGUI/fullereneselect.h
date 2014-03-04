#ifndef FULLERENESELECT_H
#define FULLERENESELECT_H

#include <QObject>

class FullereneSelect : public QObject
{
    Q_OBJECT
public:
    explicit FullereneSelect(QObject *parent = 0);

    Q_INVOKABLE QString test_invoke() const { return QStringLiteral("Yo, world!"); }


signals:

public slots:
    void load_fullerene_table();
};

#endif // FULLERENESELECT_H
