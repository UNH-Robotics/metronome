#ifdef GRAPHICS
#include <cairo.h>
#endif
#include "easylogging++.h"
#include "rapidjson/document.h"
#include "experiment/ConfigurationExecutor.hpp"
#include "rapidjson/filereadstream.h"
#include "rapidjson/istreamwrapper.h"
#include "utils/File.hpp"

//#include <iostream>
#include <cstdio>
#include <utils/Statistic.hpp>

#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/transport/TSocket.h>
#include <thrift/transport/TTransportUtils.h>

#ifdef GRAPHICS
#include <X11/Xlib.h>
#include <cairo/cairo-xlib.h>

#include <X11/Xatom.h>
#include <X11/Xutil.h>
#include <cairo.h>
#endif

#define NDEBUG

INITIALIZE_EASYLOGGINGPP

void printSplashScreen();

#ifdef GRAPHICS
cairo_surface_t* cairo_create_x11_surface(int x, int y) {
    Display* dsp;
    Drawable da;
    int screen;
    cairo_surface_t* sfc;

    if ((dsp = XOpenDisplay(NULL)) == NULL)
        exit(1);
    screen = DefaultScreen(dsp);
    da = XCreateSimpleWindow(dsp, DefaultRootWindow(dsp),
                             0, 0, x, y, 0, 0, 0);
    XSelectInput(dsp, da, ButtonPressMask | KeyPressMask);
    XMapWindow(dsp, da);

    sfc = cairo_xlib_surface_create(dsp, da,
                                    DefaultVisual(dsp, screen), x, y);
    cairo_xlib_surface_set_size(sfc, x, y);

    return sfc;
}
#endif

void thriftTest();

int main(int argc, char** argv) {
    using namespace metronome;
    printSplashScreen();
    thriftTest();

#ifdef GRAPHICS
    cairo_create_x11_surface(800,800);
#endif

    Statistic::initialize();
    if (argc == 1) {
        std::cerr << "Resource path is not provided. :: arg: " << argv[0] << std::endl;
        return 1;
    }

    std::string resourceDir{argv[1]};

    rapidjson::Document document;

    if (argc == 2) {
        std::stringstream jsonStream;

        for (std::string line; std::getline(std::cin, line);) {
            if (line.find_first_not_of(" \t\n\v\f\r") == std::string::npos) {
                break; // Terminate paring on empty line
            }

            LOG(INFO) << line;

            jsonStream << line;
        }

        rapidjson::IStreamWrapper streamWrapper{jsonStream};
        document.ParseStream(streamWrapper);
    } else {
        std::string configurationPath{argv[2]};

        if (!fileExists(configurationPath)) {
            std::cerr << "Invalid configuration file: " << configurationPath << std::endl;
        }

        std::ifstream configurationFile{configurationPath};
        rapidjson::IStreamWrapper streamWrapper{configurationFile};
        document.ParseStream(streamWrapper);
    }


    auto configuration = Configuration(std::move(document));
    LOG(INFO) << configuration.toString(); 
    const Result result = ConfigurationExecutor::executeConfiguration(configuration, resourceDir);

    LOG(INFO) << "Execution completed in " << result.planningTime / 1000000 << "ms";
    LOG(INFO) << "Path length: " << result.pathLength;
    LOG(INFO) << "Nodes :: expanded: " << result.expandedNodes << " generated: " << result.generatedNodes;

    //                for (auto action : result.actions) {
    //                    LOG(INFO) << action;
    //                }

    std::cout << "\n\nResult:" << std::endl;
    std::cout << result.getJsonString();
    std::cout << std::flush;

    return 0;
}

void thriftTest() {
    using namespace apache::thrift;
    using namespace apache::thrift::protocol;
    using namespace apache::thrift::transport;

    stdcxx::shared_ptr<TTransport> socket(new TSocket("localhost", 9090));
    stdcxx::shared_ptr<TTransport> transport(new TBufferedTransport(socket));
    stdcxx::shared_ptr<TProtocol> protocol(new TBinaryProtocol(transport));
// try {
//    transport->open();
//
//    client.ping();
//    cout << "ping()" << endl;
//
//    cout << "1 + 1 = " << client.add(1, 1) << endl;
//
//    Work work;
//    work.op = Operation::DIVIDE;
//    work.num1 = 1;
//    work.num2 = 0;
//
//    try {
//      client.calculate(1, work);
//      cout << "Whoa? We can divide by zero!" << endl;
//    } catch (InvalidOperation& io) {
//      cout << "InvalidOperation: " << io.why << endl;
       or using generated operator<<: cout << io << endl;
       or by using std::exception native method what(): cout << io.what() << endl;
//    }
//
//    work.op = Operation::SUBTRACT;
//    work.num1 = 15;
//    work.num2 = 10;
//    int32_t diff = client.calculate(1, work);
//    cout << "15 - 10 = " << diff << endl;
//
     Note that C++ uses return by reference for complex types to avoid
     costly copy construction
//    SharedStruct ss;
//    client.getStruct(ss, 1);
//    cout << "Received log: " << ss << endl;
//
//    transport->close();
//  } catch (TException& tx) {
//    cout << "ERROR: " << tx.what() << endl;
//  }   
}


void printSplashScreen() {
    std::cout << std::endl;
    std::cout << R"( ___            ___    )" << std::endl;
    std::cout << R"(|###\  ______  /###|   )" << std::endl;
    std::cout << R"(|#|\#\ \    / /#/|#|   )" << std::endl;
    std::cout << R"(|#| \#\ \  / /#/ |#|   )" << std::endl;
    std::cout << R"(|#|  \#\ \/ /#/  |#|   )" << std::endl;
    std::cout << R"(|#|      /\      |#|   )" << std::endl;
    std::cout << R"(|#|     /  \     |#|   )" << std::endl;
    std::cout << R"(|#|    /____\    |#|   )" << std::endl;
    std::cout << "---- Metronome  ----" << std::endl;
    std::cout << " When time matters!" << std::endl << std::endl;
}
