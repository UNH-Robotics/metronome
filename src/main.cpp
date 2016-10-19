#ifdef GRAPHICS
#include <cairo.h>
#endif
#include "easylogging++.h"
#include "experiment/ConfigurationExecutor.hpp"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/istreamwrapper.h"
#include "utils/File.hpp"

#include <cstdio>
#include <utils/Statistic.hpp>
#include "isolation/gridnav/main_gridnav.hpp"

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
    da = XCreateSimpleWindow(dsp, DefaultRootWindow(dsp), 0, 0, x, y, 0, 0, 0);
    XSelectInput(dsp, da, ButtonPressMask | KeyPressMask);
    XMapWindow(dsp, da);

    sfc = cairo_xlib_surface_create(dsp, da, DefaultVisual(dsp, screen), x, y);
    cairo_xlib_surface_set_size(sfc, x, y);

    return sfc;
}
#endif

int main(int argc, const char* argv[]) {
    using namespace metronome;
    printSplashScreen();

    if (argc > 3) {
        xmain(argc, argv, true);
        return 0;
    }

#ifdef GRAPHICS
    cairo_create_x11_surface(800, 800);
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

//        getchar(); // Wait for keypress

    const auto result = ConfigurationExecutor::executeConfiguration(Configuration(std::move(document)), resourceDir);

    LOG(INFO) << "Execution completed in " << result.planningTime / 1000000 << "ms";
    LOG(INFO) << "Path length: " << result.pathLength;
    LOG(INFO) << "Nodes :: expanded: " << result.expandedNodes << " generated: " << result.generatedNodes;

                for (auto action : result.actions) {
                    LOG(INFO) << action;
                }

    std::cout << "\n\nResult:" << std::endl;
    std::cout << result.getJsonString();
    std::cout << std::flush;

    return 0;
}

void printSplashScreen() {
    std::cout << std::endl;
    std::cout << " ___            ___     " << std::endl;
    std::cout << "|###\\  ______  /###|   " << std::endl;
    std::cout << "|#|\\#\\ \\    / /#/|#|   " << std::endl;
    std::cout << "|#| \\#\\ \\  / /#/ |#|   " << std::endl;
    std::cout << "|#|  \\#\\ \\/ /#/  |#|   " << std::endl;
    std::cout << "|#|      /\\      |#|   " << std::endl;
    std::cout << "|#|     /  \\     |#|   " << std::endl;
    std::cout << "|#|    /____\\    |#|   " << std::endl;
    std::cout << "---- Metronome  ----" << std::endl;
    std::cout << " When time matters!" << std::endl << std::endl;
}
