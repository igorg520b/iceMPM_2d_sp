#include "pp_mainwindow.h"
#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <iostream>
#include <omp.h>


int main(int argc, char *argv[])
{
    std::cout << "num_threads " << omp_get_max_threads() << std::endl;
    std::cout << "testing threads" << std::endl;
    int nthreads, tid;
#pragma omp parallel
    { std::cout << omp_get_thread_num(); }
    std::cout << std::endl;


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
