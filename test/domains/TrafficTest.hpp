#include "catch.hpp"
#include "domains/Traffic.hpp"
#include "easylogging++.h"
#include "experiment/ConfigurationExecutor.hpp"

namespace {
const char* json = "{\"timeLimit\" : 150000000000,\n"
                   "\"domainPath\" : "
                   "\"/input/vacuum/slalom.vw\",\n"
                   "\"domainInstanceName\" : \"Manual test instance\",\n"
                   "\"actionDuration\" : 6000000,\n"
                   "\"domainName\" : \"GRID_WORLD\",\n"
                   "\"terminationType\" : \"time\",\n"
                   "\"list\" : [1, 2],\n"
                   "\"objects\" : {\"a\": [], \"b\": {}},\n"
                   "\"algorithmName\" : \"A_STAR\"}";

const std::string resourceDir = "/home/aifs2/doylew/Public/metronome/resources";
metronome::Traffic testGrid =
        metronome::ConfigurationExecutor::extractDomain<metronome::Traffic>(metronome::Configuration(json),
                resourceDir);

bool metronome::Traffic::randomSeedFlag = false;
unsigned int metronome::Traffic::randomSeed = 1;
std::vector<std::vector<metronome::Traffic::Obstacle>> metronome::Traffic::generatedObstacles;

TEST_CASE("Traffic basic creation test", "[Traffic]") {
    metronome::Traffic traffic = testGrid;
    REQUIRE(traffic.getStartLocation() == metronome::Traffic::State(6, 0));
}
TEST_CASE("Traffic object movements test", "[Traffic]") {
    metronome::Traffic traffic = testGrid;
}
}