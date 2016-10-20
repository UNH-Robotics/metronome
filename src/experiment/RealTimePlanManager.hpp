#ifndef METRONOME_REALTIMEPLANMANAGER_HPP
#define METRONOME_REALTIMEPLANMANAGER_HPP

#include <boost/optional.hpp>
#include <chrono>
#include "MetronomeException.hpp"
#include "PlanManager.hpp"
#include "utils/TimeMeasurement.hpp"

namespace metronome {
template <typename Domain, typename Planner, typename TerminationChecker>
class RealTimePlanManager : PlanManager<Domain, Planner> {
public:
    Result plan(const Configuration& configuration, const Domain& domain, Planner& planner) {
        using namespace std::chrono;

        if (!configuration.hasMember(LOOKAHEAD_TYPE)) {
            LOG(ERROR) << "Lookahead type not found." << std::endl;
            return Result(configuration, "Missing: lookaheadType");
        }

        std::string lookaheadType{configuration.getString(LOOKAHEAD_TYPE)};

        bool dynamicLookahead;
        if (lookaheadType == LOOKAHEAD_STATIC) {
            dynamicLookahead = false;
        } else if (lookaheadType == LOOKAHEAD_DYNAMIC) {
            dynamicLookahead = true;
        } else {
            LOG(ERROR) << "Unknown lookahead type: " << lookaheadType << std::endl;
            return Result(configuration, "Unknown lookaheadType: " + lookaheadType);
        }

        std::vector<typename Domain::Action> actions;
        long long int planningTime{0};
        const auto actionDuration = configuration.getLong(ACTION_DURATION);
        auto currentState = domain.getStartState();

        TerminationChecker terminationChecker;

        long long int timeBound = actionDuration;
        long long int previousTimeBound;

        while (!domain.isGoal(currentState)) {
            auto planningIterationTime = measureNanoTime([&] {
                if (dynamicLookahead) {
                    terminationChecker.resetTo(timeBound);
                } else {
                    terminationChecker.resetTo(actionDuration);
                }

                auto actionBundles{planner.selectActions(currentState, terminationChecker)};

                previousTimeBound = timeBound;
                timeBound = 0;
                for (auto& actionBundle : actionBundles) {
                    boost::optional<typename Domain::State> nextState =
                            domain.transition(currentState, actionBundle.action);
                    if (!nextState.is_initialized()) {
                        throw MetronomeException("Invalid action. Partial plan is corrupt.");
                    }
                    currentState = nextState.get();
//                    LOG(INFO) << actionBundle.action << std::endl;
                    actions.emplace_back(actionBundle.action);
                    timeBound += actionBundle.actionDuration;
                }

            });

            //            LOG(INFO) << planningIterationTime - previousTimeBound;

            //            if (planningIterationTime > previousTimeBound + 100000) {
            //                LOG(INFO) << "Time bount violated: " << planningIterationTime - previousTimeBound << "
            //                bound: " << previousTimeBound;
            //            } else {
            //                LOG(INFO) << "Fine";
            //            }

            // Increase the total planning time after each planning iteration
            planningTime += planningIterationTime;
        }

        LOG(INFO) << "Planning: Done";

        auto pathLength = actions.size();

        std::vector<std::string> actionStrings;

        for (auto& action : actions) {
            actionStrings.push_back(action.toString());
        }

        return Result(configuration,
                planner.getExpandedNodeCount(),
                planner.getGeneratedNodeCount(),
                planningTime, // Planning time
                pathLength * configuration.getLong("actionDuration"), // Execution time
                domain.getActionDuration() + pathLength * configuration.getLong("actionDuration"), // GAT
                domain.getActionDuration(), // Idle planning time
                pathLength, // Path length
                actionStrings,
                planner.getIdentityActionCount());
    }
};
}

#endif // METRONOME_REALTIMEPLANMANAGER_HPP
