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

    QCommandLineOption resumeOption(
        QStringList() << "r" << "resume", // Option names: -r, --resume
        QCoreApplication::translate("main", "Resume simulation from the specified state file <filename>."), // Description
        QCoreApplication::translate("main", "filename") // Value name expected by the option
        );
    parser.addOption(resumeOption); // Add the option to the parser
    parser.process(a);

    const QStringList args = parser.positionalArguments();
    MainWindow w;

    if(args.size() >= 1)
    {
        QString resumeFilename;
        if (parser.isSet(resumeOption))
        {
            resumeFilename = parser.value(resumeOption);
        }

        QString parameters_file = args[0];
        w.LoadParameterFile(parameters_file, resumeFilename);
    }

    w.resize(1800,1000);
    w.move(0,0);
    w.show();
//    w.showMaximized();
    return a.exec();
}
