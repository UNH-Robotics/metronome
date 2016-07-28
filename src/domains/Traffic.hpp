#ifndef METRONOME_TRAFFIC_HPP
#define METRONOME_TRAFFIC_HPP

#include <boost/assert.hpp>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <unordered_map>
#include <utils/Location2D.hpp>
#include <vector>
#include "GridWorld.hpp"
#include "SuccessorBundle.hpp"

namespace metronome {

class Traffic {
public:
    typedef long long int Cost;
    static constexpr Cost COST_MAX = std::numeric_limits<Cost>::max();

    class Action {
    public:
        Action() : value(0) {}
        Action(unsigned int actionDuration) : value(actionDuration) {}
        constexpr char toChar() const {
            if (value == 1) {
                return 'N';
            } else if (value == 2) {
                return 'E';
            } else if (value == 3) {
                return 'S';
            } else if (value == 4) {
                return 'W';
            } else if (value == 5) {
                return '0';
            } else {
                return '~';
            }
        }
        const std::string toString() const {
            std::string s;
            s.push_back(toChar());
            return s;
        }

    private:
        unsigned int value;
    };

    class Obstacle {
    public:
        Obstacle() : x(0), y(0), xVelocity(0), yVelocity(0) {}
        Obstacle(const unsigned x, const unsigned y) : x(x), y(y) {
            if (metronome::Traffic::randomSeedFlag) {
                std::srand(metronome::Traffic::randomSeed);
            } else {
                std::srand(std::time(0));
            }

            if (metronome::Traffic::generatedObstacles[x][y] != nullptr) {
                xVelocity = std::rand() % 3;
                yVelocity = std::rand() % 3;
                metronome::Traffic::generatedObstacles[x][y] = this;
            } else {
                xVelocity = metronome::Traffic::generatedObstacles[x][y]->getXVelocity();
                yVelocity = metronome::Traffic::generatedObstacles[x][y]->getYVelocity();
            }
        }
        Obstacle& operator=(Obstacle toCopy) {
            swap(*this, toCopy);
            return *this;
        }
        unsigned int getX() const { return x; }
        unsigned int getY() const { return y; }
        std::size_t hash() const { return x ^ y << 16 ^ y >> 16; }
        int getXVelocity() const { return xVelocity; }
        int getYVelocity() const { return yVelocity; }
        bool operator==(const Obstacle& obstacle) const {
            return x == obstacle.x && y == obstacle.y && xVelocity == obstacle.xVelocity &&
                    yVelocity == obstacle.yVelocity;
        }

    private:
        unsigned int x;
        unsigned int y;
        int xVelocity;
        int yVelocity;

        friend void swap(Obstacle& first, Obstacle& second) {
            using std::swap;
            swap(first.x, second.x);
            swap(first.y, second.y);
        }
    };

    class State {
    public:
        State() : x(0), y(0) {}
        State(const unsigned int x, const unsigned int y) : x(x), y(y) {}
        State& operator=(State toCopy) {
            swap(*this, toCopy);
            return *this;
        }

        unsigned int getX() const { return x; }
        unsigned int getY() const { return y; }
        std::size_t hash() const { return x ^ y << 16 ^ y >> 16; }
        bool operator==(const State& state) const { return x == state.x && y == state.y; }
        const std::string toString() const {
            std::string string("x: ");
            return string + std::to_string(x) + " y: " + std::to_string(y);
        }

    private:
        unsigned int x;
        unsigned int y;

        friend void swap(State& first, State& second) {
            using std::swap;
            swap(first.x, second.x);
            swap(first.y, second.y);
        }
    };

    Traffic(const Configuration& configuration, std::istream& input) {
        if (!configuration.hasMember(ACTION_DURATION)) {
            throw MetronomeException("No value provided.");
        }
        actionDuration = configuration.getLong(ACTION_DURATION);
        deadCost = actionDuration * 2;
        unsigned int currentHeight = 0;
        unsigned int currentWidth = 0;
        unsigned int currentIndex = 0;
        std::string line;
        char* end;
        getline(input, line);
        std::stringstream convertWidth(line);
        if (std::strtol(line.c_str(), &end, 10) == 0) {
            throw MetronomeException("Traffic first line must be a number.");
        }
        convertWidth >> width;
        getline(input, line);
        if (std::strtol(line.c_str(), &end, 10) == 0) {
            throw MetronomeException("Traffic second line must be a number.");
        }
        std::stringstream convertHeight(line);
        convertHeight >> height;

        for (auto i = 0; i < width; ++i) {
            for (auto j = 0; j < height; ++j) {
                obstacles[i][j] = false;
                bunkers[i][j] = false;
            }
        }

        boost::optional<State> tempStartState;
        boost::optional<State> tempGoalState;

        while (getline(input, line)) {
            for (auto it = line.cbegin(); it != line.cend(); ++it) {
                // do something for IO
                if (*it == '@') { // find the start location
                    tempStartState = State(currentWidth, currentHeight);
                } else if (*it == '*') { // find the goal location
                    tempGoalState = State(currentWidth, currentHeight);
                } else if (*it == '#') { // store the objects
                    Obstacle newObstacle = Obstacle(currentWidth, currentHeight);
                    metronome::Traffic::generatedObstacles[currentWidth][currentHeight] = &newObstacle;
                    obstacleIndices.push_back(metronome::Location2D(currentWidth, currentHeight));
                    obstacles[currentWidth][currentHeight] = true;
                } else if (*it == '$') {
                    bunkers[currentWidth][currentHeight] = true;
                } else {
                    // open cell!
                }
                ++currentWidth;
                ++currentIndex;
            }
            if (currentWidth != width) {
                throw MetronomeException("Traffic is not complete. Width doesn't match the input configuration.");
            }
            currentWidth = 0;
            ++currentHeight;
        }

        if (currentHeight != height) {
            throw MetronomeException("Traffic is not complete. Height doesn't match input configuration.");
        }

        if (!tempStartState.is_initialized() || !tempGoalState.is_initialized()) {
            throw MetronomeException("Traffic unknown start or goal location. Start or goal location is not defined.");
        }

        startLocation = tempStartState.get();
        startLocation = tempGoalState.get();
    }
    const State transition(const State& state, const Action& action) {
        moveObstacles();

        if (action.toChar() == 'N') {
            State newState = State(state.getX(), state.getY() - 1);
            if (isLegalLocation(newState)) {
                return newState;
            }
        } else if (action.toChar() == 'E') {
            State newState = State(state.getX() + 1, state.getY());
            if (isLegalLocation(newState)) {
                return newState;
            }
        } else if (action.toChar() == 'S') {
            State newState = State(state.getX(), state.getY() + 1);
            if (isLegalLocation(newState)) {
                return newState;
            }
        } else if (action.toChar() == 'W') {
            State newState = State(state.getX() - 1, state.getY());
            if (isLegalLocation(newState)) {
                return newState;
            }
        } else if (action.toChar() == '0' && bunkers[state.getX()][state.getY()]) {
            State newState = State(state.getX(), state.getY());
            if (isLegalLocation(newState)) {
                return newState;
            }
        }

        return state;
    }
    const bool isObstacle(const State& location) const { return obstacles[location.getX()][location.getY()]; }
    const bool isLegalLocation(const State& location) const {
        return location.getX() < width && location.getY() < height && !isObstacle(location);
    }
    /*
     * this needs to be fixed....
     */
    std::vector<SuccessorBundle<Traffic>> successors(State state) {
        std::vector<SuccessorBundle<Traffic>> successors;

        unsigned int actions[] = {5, 4, 3, 2, 1};

        for (auto a : actions) {
            State newState = transition(state, Action(a));
            for (auto obstacleIndex : obstacleIndices) {
                if (obstacleIndex.x == newState.getX() && obstacleIndex.y == newState.getY()) {
                    successors.push_back(SuccessorBundle<Traffic>{newState, a, deadCost});
                } else {
                    successors.push_back(SuccessorBundle<Traffic>{newState, a, actionDuration});
                }
            }
        }

        return successors;
    }

    const State getStartLocation() const { return startLocation; }

    Cost heuristic(const State& state) { return heuristic(state); }

    const bool isGoal(const State& location) {
        return goalLocation.getX() == location.getX() && goalLocation.getY() == location.getY();
    }

    const bool isStart(const State& location) {
        return startLocation.getX() == location.getX() && startLocation.getY() == location.getY();
    }

private:
    void moveObstacles() {
        for (auto obstacleIndex : obstacleIndices) {
            Obstacle* curObstacle = generatedObstacles[obstacleIndex.x][obstacleIndex.y];
            int xVelocity = curObstacle->getXVelocity();
            int yVelocity = curObstacle->getYVelocity();

            unsigned int newXLocation = curObstacle->getX() + xVelocity;
            unsigned int newYLocation = curObstacle->getY() + yVelocity;

            // make sure new location is on the grid otherwise bounce
            if (newXLocation > width) {
                newXLocation = curObstacle->getX() + (xVelocity * -1);
            }
            if (newYLocation > height) {
                newYLocation = curObstacle->getY() + (yVelocity * -1);
            }
            // if the new location is a bunker bounce or in another obstacle
            if (bunkers[newXLocation][newYLocation] || obstacles[newXLocation][newYLocation]) {
                newXLocation = curObstacle->getX() + (xVelocity * -1);
                newYLocation = curObstacle->getY() + (yVelocity * -1);
            }
            // update the obstacleIndices
            obstacleIndex.x = newXLocation;
            obstacleIndex.y = newYLocation;
            // update the obstacle bit
            obstacles[obstacleIndex.x][obstacleIndex.y] = false;
            obstacles[newXLocation][newYLocation] = true;
            // update the generatedObstacles
            generatedObstacles[obstacleIndex.x][obstacleIndex.y] = nullptr;
            Obstacle newObstacle = Obstacle(newXLocation, newYLocation);
            generatedObstacles[newXLocation][newYLocation] = &newObstacle;
        }
    }

    /*
     * randomSeedFlag <- notifies the ctor if we are using internal seed or user-defined
     * randomSeed <- the seed we use for generating the velocities of the obstacles
     * width <- how wide the world is
     * height <- how tall the world is
     * obstacleIndices <- where the obstacles are using direct addressing *NOTE: vector is a container for a
     * bunkerIndices <- where the bunkers are using direct addressing      *NOTE: dynamic array!!!!!!!!!!!!!!
     * obstacles <- bit vector using direct addressing of the obstacles
     * bunkers <- bit vector using direct addressing of the bunkers ** NOTE: are these redundant?***
     * startLocation <- where the agent starts
     * goalLocation <- where the agent needs to go
     * deadCost <- to calculate if we are in a dead state if the cost is twice the actionDuration then we pronounce
     * it
     * dead
     * generatedObstacles <- our cache trick to insure velocity consistency when generating obstacles
     */
    static bool randomSeedFlag;
    static unsigned int randomSeed;
    unsigned int width;
    unsigned int height;
    std::vector<metronome::Location2D> obstacleIndices;
    std::vector<std::vector<bool>> obstacles;
    std::vector<std::vector<bool>> bunkers;
    State startLocation = State();
    State goalLocation = State();
    Cost actionDuration;
    Cost deadCost;
    static std::vector<std::vector<Obstacle*>> generatedObstacles;
};
}
#endif // METRONOME_TRAFFIC_HPP