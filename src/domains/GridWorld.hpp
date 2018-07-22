#ifndef VACUUM_WORLD_HPP
#define VACUUM_WORLD_HPP

#include <boost/assert.hpp>
#include <boost/optional.hpp>
#include <cstdlib>
#include <experiment/Configuration.hpp>
#include <functional>
#include <limits>
#include <ostream>
#include <unordered_set>
#include <utils/Hasher.hpp>
#include <vector>
#include "MetronomeException.hpp"
#include "SuccessorBundle.hpp"

namespace metronome {
class GridWorld {
public:
    typedef long long int Cost;
    static constexpr Cost COST_MAX = std::numeric_limits<Cost>::max();

    class Action {
    public:
        Action() : label{'~'} {};

        static const std::vector<Action>& getActions() {
            static const std::vector<Action> actions{
                Action('N'),
                Action('S'),
                Action('W'),
                Action('E'),
                Action('A'),
                Action('B'),
                Action('C'),
                Action('D'),
            };
            return actions;
        }

        bool operator==(const Action& rhs) const { return label == rhs.label; }

        bool operator!=(const Action& rhs) const { return !(rhs == *this); }

        char toChar() const { return label; }

        static Action getIdentity() { return Action('0'); }

        std::string toString() const { return std::string{label}; }

        friend std::ostream&
        operator<<(std::ostream& os, const Action& action) {
            os << action.label;
            return os;
        }

    private:
        Action(char label) : label{label} {}
        char label;
    };

    class State {
    public:
        State()
            : x(0),
              y(0) {}
        State(unsigned int x, unsigned int y) : x{x},
                                                y{y} {}
        /*Standard getters for the State(x,y)*/
        unsigned int getX() const { return x; }
        unsigned int getY() const { return y; }
        std::size_t hash() const { return x ^ y << 16 ^ y >> 16; }
        bool operator==(const State& state) const {
            return x == state.x && y == state.y;
        }
        bool operator!=(const State& state) const {
            return x != state.x || y != state.y;
        }
        State& operator=(State toCopy) {
            swap(*this, toCopy);
            return *this;
        }
        std::string toString() {
            std::string string{"x: "};
            return string + std::to_string(x) + " y: " + std::to_string(y);
        }

        friend std::ostream&
        operator<<(std::ostream& stream, const GridWorld::State& state) {
            stream << "x: " << state.getX() << " y: " << state.getY();
            return stream;
        }

    private:
        /*State(x,y) representation*/
        unsigned int x;
        unsigned int y;
        /*Function facilitating operator==*/
        friend void swap(State& first, State& second) {
            using std::swap;
            swap(first.x, second.x);
            swap(first.y, second.y);
        }
    };

    /*Entry point for using this Domain*/
    GridWorld(const Configuration& configuration, std::istream& input)
        : actionDuration{configuration.getLongOrThrow(ACTION_DURATION)},
          octileMovement{
              configuration.hasMember(OCTILE_MOVEMENT) ? configuration
                  .getBool(OCTILE_MOVEMENT) : false},
          octileActionDuration{(const Cost) sqrt(2.0) * actionDuration} {
        obstacles =
            std::unordered_set<State, typename metronome::Hasher<State>>{};
        unsigned int currentHeight = 0;
        unsigned int currentWidth = 0;
        std::string line;
        char* end;
        getline(input, line); // get the width
        std::stringstream convertWidth(line);
        if (std::strtol(line.c_str(), &end, 10) == 0) {
            throw MetronomeException("GridWorld first line must be a number.");
        }
        convertWidth >> width;
        getline(input, line); // get the height
        if (std::strtol(line.c_str(), &end, 10) == 0) {
            throw MetronomeException("GridWorld second line must be a number.");
        }
        std::stringstream convertHeight(line);
        convertHeight >> height;

        boost::optional<State> tempStarState;
        boost::optional<State> tempGoalState;

        while (getline(input, line) && currentHeight < height) {
            for (auto it = line.cbegin(); it != line.cend(); ++it) {
                if (*it == '@') { // find the start location
                    tempStarState = State(currentWidth, currentHeight);
                } else if (*it == '*') { // find the goal location
                    tempGoalState = State(currentWidth, currentHeight);
                } else if (*it == '#') { // store the objects
                    State object = State(currentWidth, currentHeight);
                    obstacles.insert(object);
                } else if (*it == '_') {
                    // its an open cell nothing needs to be done
                } else if (*it == 'v') {
                    State object = State(currentWidth, currentHeight);
                    directions[object] = 'S';
                } else if (*it == '^') {
                    State object = State(currentWidth, currentHeight);
                    directions[object] = 'N';
                } else if (*it == '<') {
                    State object = State(currentWidth, currentHeight);
                    directions[object] = 'W';
                } else if (*it == '>') {
                    State object = State(currentWidth, currentHeight);
                    directions[object] = 'E';
                }
                ++currentWidth; // at the end of the character parse move along
            }
            if (currentWidth != width) {
                throw MetronomeException(
                    "GridWorld is not complete. Width doesn't match input configuration.");
            }

            currentWidth = 0; // restart character parse at beginning of line
            ++currentHeight; // move down one line in charadter parse
        }

        if (currentHeight != height) {
            throw MetronomeException(
                "GridWorld is not complete. Height doesn't match input configuration.");
        }

        if (!tempStarState.is_initialized()
            || !tempGoalState.is_initialized()) {
            throw MetronomeException(
                "Traffic unknown start or goal location. Start or goal location is not defined.");
        }

        startLocation = tempStarState.get();
        goalLocation = tempGoalState.get();
    }

    /**
     * Transition to a new state from the given state applying the given action.
     *
     * @return the original state if the transition is not possible
     */
    boost::optional<State>
    transition(const State& sourceState, const Action& action) const {
        boost::optional<State> targetState;

        if (action.toChar() == 'N') {
            State newState = State(sourceState.getX(), sourceState.getY() - 1);
            if (isLegalLocation(newState)) {
                targetState = newState;
            }
        } else if (action.toChar() == 'E') {
            State newState = State(sourceState.getX() + 1, sourceState.getY());
            if (isLegalLocation(newState)) {
                targetState = newState;
            }
        } else if (action.toChar() == 'S') {
            State newState = State(sourceState.getX(), sourceState.getY() + 1);
            if (isLegalLocation(newState)) {
                targetState = newState;
            }
        } else if (action.toChar() == 'W') {
            State newState = State(sourceState.getX() - 1, sourceState.getY());
            if (isLegalLocation(newState)) {
                targetState = newState;
            }
        } else if (action.toChar() == 'A') {
            State newState =
                State(sourceState.getX() + 1, sourceState.getY() + 1);
            if (isLegalLocation(newState)) {
                targetState = newState;
            }
        } else if (action.toChar() == 'B') {
            State newState =
                State(sourceState.getX() + 1, sourceState.getY() - 1);
            if (isLegalLocation(newState)) {
                targetState = newState;
            }
        } else if (action.toChar() == 'C') {
            State newState =
                State(sourceState.getX() - 1, sourceState.getY() - 1);
            if (isLegalLocation(newState)) {
                targetState = newState;
            }
        } else if (action.toChar() == 'D') {
            State newState =
                State(sourceState.getX() - 1, sourceState.getY() + 1);
            if (isLegalLocation(newState)) {
                targetState = newState;
            }
        } else if (action.toChar() == '0') {
            return boost::make_optional(sourceState);
        }

        return targetState;
    }

    /*Validating a goal state*/
    bool isGoal(const State& location) const {
        return location.getX() == goalLocation.getX()
            && location.getY() == goalLocation.getY();
    }

    /*Validating an obstacle state*/
    bool isObstacle(const State& location) const {
        return obstacles.find(location) != obstacles.end();
    }

    /*Validating the agent can visit the state*/
    bool isLegalLocation(const State& location) const {
        return location.getX() < width && location.getY() < height
            && !isObstacle(location);
    }

    /*Standard getters for the (width,height) of the domain*/
    unsigned int getWidth() { return width; }
    unsigned int getHeight() { return height; }
    /*Adding an obstacle to the domain*/
    bool addObstacle(const State& toAdd) {
        if (isLegalLocation(toAdd)) {
            obstacles.insert(toAdd);
            return true;
        }
        return false;
    }

    std::vector<State>::size_type getNumberObstacles() {
        return obstacles.size();
    }

    State getStartState() const { return startLocation; }

    bool isStart(const State& state) const {
        return state.getX() == startLocation.getX()
            && state.getY() == startLocation.getY();
    }

    Cost distance(const State& state) const {
        return distance(state, goalLocation);
    }

    Cost distance(const State& sourceState, const State& targetState) const {
        const auto xDistance = std::abs(static_cast<int>(sourceState.getX()) -
            static_cast<int>(targetState.getX()));
        const auto yDistance = std::abs(static_cast<int>(sourceState.getY()) -
            static_cast<int>(targetState.getY()));

        const auto totalDistance = xDistance + yDistance;
        return static_cast<Cost>(totalDistance);
    }

    Cost heuristic(const State& state) const {
        return distance(state) * actionDuration;
    }

    Cost heuristic(const State& sourceState, const State& targetState) const {
        return distance(sourceState, targetState) * actionDuration;
    }

    bool safetyPredicate(const State& state) const { return true; }

    std::vector<SuccessorBundle<GridWorld>>
    successors(const State& state) const {
        std::vector<SuccessorBundle<GridWorld>> successors;

        auto iterator = directions.find(state);
        if (iterator != directions.end()) {
            const auto direction = iterator->second;
            if (direction == 'S') {
                addValidSuccessor(successors,
                                  state,
                                  0,
                                  1,
                                  Action::getActions()[1],
                                  actionDuration); // South
            } else if (direction == 'N') {
                addValidSuccessor(successors,
                                  state,
                                  0,
                                  -1,
                                  Action::getActions()[0],
                                  actionDuration); // North
            } else if (direction == 'W') {
                addValidSuccessor(successors,
                                  state,
                                  -1,
                                  0,
                                  Action::getActions()[2],
                                  actionDuration); // West
            } else if (direction == 'E') {
                addValidSuccessor(successors,
                                  state,
                                  1,
                                  0,
                                  Action::getActions()[3],
                                  actionDuration); // East
            }
            return successors;
        }

        addValidSuccessor(successors,
                          state,
                          0,
                          1,
                          Action::getActions()[1],
                          actionDuration); // South
        addValidSuccessor(successors,
                          state,
                          0,
                          -1,
                          Action::getActions()[0],
                          actionDuration); // North
        addValidSuccessor(successors,
                          state,
                          -1,
                          0,
                          Action::getActions()[2],
                          actionDuration); // West
        addValidSuccessor(successors,
                          state,
                          1,
                          0,
                          Action::getActions()[3],
                          actionDuration); // East

        if (octileMovement) {
            addValidSuccessor(successors,
                              state,
                              1,
                              1,
                              Action::getActions()[4],
                              octileActionDuration);
            addValidSuccessor(successors,
                              state,
                              1,
                              -1,
                              Action::getActions()[5],
                              octileActionDuration);
            addValidSuccessor(successors,
                              state,
                              -1,
                              -1,
                              Action::getActions()[6],
                              octileActionDuration);
            addValidSuccessor(successors,
                              state,
                              -1,
                              1,
                              Action::getActions()[7],
                              octileActionDuration);
        }

        return successors;
    }

    void visualize(std::ostream& display) const {
        for (unsigned int i = 0; i < height; ++i) {
            for (unsigned int j = 0; j < width; ++j) {
                if (startLocation.getX() == j && startLocation.getY() == i) {
                    display << '@';
                } else if (goalLocation.getX() == j
                    && goalLocation.getY() == i) {
                    display << '*';
                } else if (isObstacle(State(j, i))) {
                    display << '#';
                } else {
                    display << '_';
                }
            }
            display << "\n";
        }
        display << "\n";
    }

    Cost getActionDuration() const { return actionDuration; }

private:
    void addValidSuccessor(std::vector<SuccessorBundle<GridWorld>>& successors,
                           const State& sourceState,
                           const signed char relativeX,
                           const signed char relativeY,
                           const Action& action,
                           const Cost cost) const {
        const boost::optional<State>
            & successor = getSuccessor(sourceState, relativeX, relativeY);
        if (successor.is_initialized()) {
            successors.emplace_back(successor.get(), action, actionDuration);
        }
    }

    boost::optional<State> getSuccessor(const State& sourceState,
                                        signed char relativeX,
                                        signed char relativeY) const {
        State newState = State(sourceState.getX() + relativeX,
                               sourceState.getY() + relativeY);
        if (isLegalLocation(newState)) {
            return boost::make_optional(newState);
        }

        return boost::none;
    }

    /*
     * maxActions <- maximum number of actions
     * width/height <- internal size representation of world
     * obstacles <- stores locations of the objects in world
     * dirtyCells <- stores locations of dirt in the world
     * startLocation <- where the agent begins
     * goalLocation <- where the agent needs to end up
     * initalAmountDirty <- how many cells are dirty
     * initialCost <- constant cost value
     * obstacles <- stores references to obstacles
     */
    unsigned int width;
    unsigned int height;
    std::unordered_set<State, typename metronome::Hasher<State>> obstacles;
    std::unordered_map<State, char, typename metronome::Hasher<State>>
        directions;
    State startLocation{};
    State goalLocation{};
    const Cost actionDuration;
    const bool octileMovement;
    const Cost octileActionDuration;
};
}
#endif
