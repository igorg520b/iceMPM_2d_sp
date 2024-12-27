#include "pp_mainwindow.h"
#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <iostream>


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QApplication::setApplicationName("MPM Postprocessor");
    QApplication::setApplicationVersion("1.0");

    QCommandLineParser parser;
    parser.setApplicationDescription("Post-processing the HDF5 simulation output");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("frame", QCoreApplication::translate("main", "Data frame to visualize"));

    parser.process(a);

    const QStringList args = parser.positionalArguments();
    PPMainWindow w;

/*    if(args.size() == 1)
    {
        QString parameters_file = args[0];
        w.LoadParameterFile(parameters_file);
    }
*/
    w.resize(1400,900);
//    w.show();
    w.showMaximized();
    return a.exec();
}
