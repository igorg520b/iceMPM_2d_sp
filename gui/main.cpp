#include "mainwindow.h"
#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <iostream>


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QApplication::setApplicationName("iceMPM");
    QApplication::setApplicationVersion("1.2");

    QCommandLineParser parser;
    parser.setApplicationDescription("MPM simulation of ice with GUI");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("parameters", QCoreApplication::translate("main", "JSON parameter file"));

    parser.process(a);

    const QStringList args = parser.positionalArguments();
    MainWindow w;

    if(args.size() == 1)
    {
        QString parameters_file = args[0];
        w.LoadParameterFile(parameters_file);
    }

    w.resize(1800,1000);
    w.move(0,0);
    w.show();
//    w.showMaximized();
    return a.exec();
}
