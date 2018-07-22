#ifndef METRONOME_ENVELOPESEARCH_HPP
#define METRONOME_ENVELOPESEARCH_HPP

#include <DynamicPriorityQueue.hpp>
#include <MemoryConfiguration.hpp>
#include <domains/SuccessorBundle.hpp>
#include <experiment/Configuration.hpp>
#include <unordered_map>
#include <utils/Hasher.hpp>
#include <utils/PriorityQueue.hpp>
#include <utils/StaticVector.hpp>
#include "OnlinePlanner.hpp"

namespace metronome {

template <typename Domain, typename TerminationChecker>
class EnvelopeSearch final : public OnlinePlanner<Domain, TerminationChecker> {
public:
    typedef typename Domain::State State;
    typedef typename Domain::Action Action;
    typedef typename Domain::Cost Cost;
    typedef typename OnlinePlanner<Domain, TerminationChecker>::ActionBundle ActionBundle;

    EnvelopeSearch(const Domain& domain, const Configuration&)
            : domain{domain}, pseudoFOpenList{PseudoFComparator{*this}} {
        // Initialize hash table
        nodes.max_load_factor(1);
        nodes.reserve(Memory::NODE_LIMIT);
    }

    std::vector<ActionBundle> selectActions(const State& sourceState, TerminationChecker& terminationChecker) override {
        initialize(sourceState);

        if (domain.isGoal(sourceState)) {
            return {};
        }

        forwardSearch(terminationChecker.limit(greedyResourceRatio + pseudoFResourceRatio));
        backup(terminationChecker.limit(backupLimit));

        //        const auto bestNode = explore(startState, terminationChecker);
        //
        //        return extractPath(bestNode, nodes[startState]);
    }

private:
    class Edge;

    class Node {
    public:
        Node(Node* parent, const State& state, Action action, Cost g, Cost h, bool open, unsigned int iteration = 0)
                : parent{parent},
                  state{state},
                  action{std::move(action)},
                  g{g},
                  h{h},
                  open{open},
                  iteration{iteration} {}

        Cost f() const { return g + h; }

        unsigned long hash() const { return state.hash(); }

        bool operator==(const Node& node) const { return state == node.state; }

        std::string toString() const {
            std::ostringstream stream;
            stream << "s: " << state << " g: " << g << " h: " << h << " f: " << f() << " a: " << action << " p: ";
            if (parent == nullptr) {
                stream << "None";
            } else {
                stream << parent->state;
            }
            stream << (open ? " Open" : " Not Open");
            return stream.str();
        }

        /** Index used by the priority queue */
        mutable unsigned int index;
        /** Parent node */
        Node* parent;
        /** Internal state */
        const State state;
        /** Action that led to the current node from the parent node */
        Action action;
        /** Cost from the root node */
        Cost g;
        /** Heuristic cost of the node */
        Cost h;
        /** True if the node is in the open list */
        bool open;
        /** Last iteration when the node was updated */
        unsigned int iteration;
        /** List of all the predecessors that were discovered in the current exploration phase. */
        std::vector<Edge> predecessors;

        // Envelope Search specific members

        Cost waveHeuristic;

        std::size_t waveCounter;

        std::size_t pseudoFOpenIndex = std::numeric_limits<std::size_t>::max();
        std::size_t heuristicOpenIndex = std::numeric_limits<std::size_t>::max();
        std::size_t backupFrontierIndex = std::numeric_limits<std::size_t>::max();

        std::size_t expanded = 0;
        std::size_t generated = 0;
    };

    struct PseudoFIndexFunction {
        static std::size_t& operator()(Node* node) { return node->pseudoFOpenIndex; }
    };

    struct HeuristicIndexFunction {
        static std::size_t& operator()(Node* node) { return node->heuristicOpenIndex; }
    };

    struct BackupIndexFunction {
        static std::size_t& operator()(Node* node) { return node->backupFrontierIndex; }
    };

    struct HeuristicComparator {
        static int operator()(Node* lhs, Node* rhs) {
            if (lhs->h < rhs->h)
                return -1;
            if (lhs->h > rhs->h)
                return 1;

            return 0;
        }
    };

    struct PseudoFComparator {
        const EnvelopeSearch& envelopeSearch;

        PseudoFComparator(const EnvelopeSearch& envelopeSearch) : envelopeSearch(envelopeSearch) {}

        int operator()(Node* lhs, Node* rhs) {
            Cost lhsPseudoG = envelopeSearch.domain.heuristic(envelopeSearch.agentState, lhs->state) * pseudoGWeight;

            Cost rhsPseudoG = envelopeSearch.domain.heuristic(envelopeSearch.agentState, rhs->state) * pseudoGWeight;

            Cost lhsPseudoF = lhs->h + lhsPseudoG;
            Cost rhsPseudoF = rhs->h + rhsPseudoG;

            if (lhsPseudoF < rhsPseudoF)
                return -1;
            if (lhsPseudoF > rhsPseudoF)
                return 1;
            if (lhs->h < rhs->h)
                return -1;
            if (lhs->h > rhs->h)
                return 1;
            if (lhs->generated < rhs->generated)
                return -1;
            if (lhs->generated > rhs->generated)
                return 1;

            return 0;
        }
    };

    struct WaveFComparator {
        const EnvelopeSearch& envelopeSearch;

        WaveFComparator(const EnvelopeSearch& envelopeSearch) : envelopeSearch(envelopeSearch) {}

        int operator()(Node* lhs, Node* rhs) {
            Cost lhsWaveF = envelopeSearch.domain.heuristic(envelopeSearch.agentState, lhs->state) + lhs->waveHeuristic;

            Cost rhsWaveF = envelopeSearch.domain.heuristic(envelopeSearch.agentState, rhs->state) + rhs->waveHeuristic;

            if (lhsWaveF < rhsWaveF)
                return -1;
            if (lhsWaveF > rhsWaveF)
                return 1;
            if (lhs->waveHeuristic < rhs->waveHeuristic)
                return -1;
            if (lhs->waveHeuristic > rhs->waveHeuristic)
                return 1;

            return 0;
        }
    };

    struct WaveComparator {
        const EnvelopeSearch& envelopeSearch;

        WaveComparator(const EnvelopeSearch& envelopeSearch) : envelopeSearch(envelopeSearch) {}

        int operator()(Node* lhs, Node* rhs) {
            Cost lhsWaveF = envelopeSearch.domain.heuristic(envelopeSearch.agentState, lhs->state) + lhs->waveHeuristic;

            Cost rhsWaveF = envelopeSearch.domain.heuristic(envelopeSearch.agentState, rhs->state) + rhs->waveHeuristic;

            if (lhs->waveCounter < rhs->waveCounter)
                return -1;
            if (lhs->waveCounter > rhs->waveCounter)
                return 1;
            if (lhs->h < rhs->h)
                return -1;
            if (lhs->h > rhs->h)
                return 1;

            return 0;
        }
    };

    class Edge {
    public:
        Edge(Node* predecessor, Action action, Cost actionCost)
                : predecessor{predecessor}, action{action}, actionCost{actionCost} {}

        Node* predecessor;
        const Action action;
        const Cost actionCost;
    };

    enum class Phase { GOAL_SEARCH, GOAL_BACKUP, PATH_IMPROVEMENT };
    enum class QueueSelector { HEURISTIC, PSEUDO_F };

    void initialize(const State& sourceState) {
        Node*& node = nodes[sourceState];
        if (node == nullptr) {
            node = nodePool.construct(nullptr, sourceState, Action(), 0, domain.heuristic(sourceState), true);

            addToOpenList(node);
        }

        agentState = sourceState;
    }

    void forwardSearch(TerminationChecker& terminationChecker) {
        activeQueue = QueueSelector::HEURISTIC;
        explore(terminationChecker.limit(greedyResourceRatio));
        activeQueue = QueueSelector::PSEUDO_F;
        explore(terminationChecker.limit(greedyResourceRatio + pseudoFResourceRatio));
    }

    void learn(TerminationChecker& terminationChecker) {
        ++iterationCounter;

        // Reorder the open list based on the heuristic values
        openList.reorder(hComparator);

        while (!terminationChecker.reachedTermination() && openList.isNotEmpty()) {
            auto currentNode = popOpenList();
            currentNode->iteration = iterationCounter;

            Cost currentHeuristicValue = currentNode->h;

            // update heuristic actionDuration of each predecessor
            for (auto predecessor : currentNode->predecessors) {
                Node* predecessorNode = predecessor.predecessor;

                if (predecessorNode->iteration == iterationCounter && !predecessorNode->open) {
                    // This node was already learned and closed in the current iteration
                    continue;
                    // TODO Review this. This could be incorrect if the action costs are not uniform
                }

                if (!predecessorNode->open) {
                    // This node is not open yet, because it was not visited in the current planning iteration

                    predecessorNode->h = currentHeuristicValue + predecessor.actionCost;
                    assert(predecessorNode->iteration == iterationCounter - 1);
                    predecessorNode->iteration = iterationCounter;

                    addToOpenList(*predecessorNode);
                } else if (predecessorNode->h > currentHeuristicValue + predecessor.actionCost) {
                    // This node was visited in this learning phase, but the current path is better then the previous
                    predecessorNode->h = currentHeuristicValue + predecessor.actionCost;
                    openList.update(*predecessorNode);
                }
            }
        }
    }

    const Node* explore(TerminationChecker& terminationChecker) {
        ++iterationCounter;

        Planner::incrementGeneratedNodeCount();

        while (!terminationChecker.reachedTermination() && openList.isNotEmpty()) {
            const auto topNode = topOfOpen();
            if (topNode != nullptr and domain.isGoal(topNode->state)) {
                return topNode;
            }

            Node* const currentNode = popOpenList();

            terminationChecker.notifyExpansion();
            expandNode(currentNode);
        }

        return openList.top();
    }

    void expandNode(Node* sourceNode) {
        Planner::incrementExpandedNodeCount();

        for (auto successor : domain.successors(sourceNode->state)) {
            auto successorState = successor.state;

            Node*& successorNode = nodes[successorState];

            if (successorNode == nullptr) {
                successorNode = createNode(sourceNode, successor);
            }

            // If the node is outdated it should be updated.
            if (successorNode->iteration != iterationCounter) {
                successorNode->iteration = iterationCounter;
                successorNode->predecessors.clear();
                successorNode->g = Domain::COST_MAX;
                successorNode->open = false; // It is not on the open list yet, but it will be
                // parent, action, and actionCost is outdated too, but not relevant.
            }

            // Add the current state as the predecessor of the child state
            successorNode->predecessors.emplace_back(sourceNode, successor.action, successor.actionCost);

            // Set iteration

            // only generate those state that are not visited yet or whose cost value are lower than this path
            Cost successorGValueFromCurrent{sourceNode->g + successor.actionCost};
            if (successorNode->g > successorGValueFromCurrent) {
                successorNode->g = successorGValueFromCurrent;
                successorNode->parent = sourceNode;
                successorNode->action = successor.action;

                if (!successorNode->open) {
                    addToOpenList(*successorNode);
                } else {
                    openList.update(*successorNode);
                }
            }
        }
    }

    Node* createNode(Node* sourceNode, SuccessorBundle<Domain> successor) {
        Planner::incrementGeneratedNodeCount();
        return nodePool->construct(Node{sourceNode,
                successor.state,
                successor.action,
                domain.COST_MAX,
                domain.heuristic(successor.state),
                true});
    }

    void clearOpenList() {
        openList.forEach([](Node* node) { node->open = false; });
        openList.clear();
    }

    Node* popOpenList() {
        // TODO
        if (openList.isEmpty()) {
            throw MetronomeException("Open list was empty, goal not reachable.");
        }

        Node* node = openList.pop();
        node->open = false;
        return node;
    }

    bool isOpenEmpty() {
        return false; // TODO
    }

    Node* topOfOpen() {
        return nullptr; // TODO
    }

    void addToOpenList(Node* node) {
        // TODO
        //        node.open = true;
        //        openList.push(node);
    }

    std::vector<ActionBundle> extractPath(const Node* targetNode, const Node* sourceNode) const {
        if (targetNode == sourceNode || targetNode == nullptr) {
            //                        LOG(INFO) << "We didn't move:" << sourceNode->toString();
            return std::vector<ActionBundle>();
        }

        std::vector<ActionBundle> actionBundles;
        auto currentNode = targetNode;

        while (currentNode != sourceNode) {
            // The g difference of the child and the parent gives the action cost from the parent
            actionBundles.emplace_back(currentNode->action, currentNode->g - currentNode->parent->g);
            currentNode = currentNode->parent;
        }

        std::reverse(actionBundles.begin(), actionBundles.end());
        return actionBundles;
    }

    static int fComparator(const Node& lhs, const Node& rhs) {
        if (lhs.f() < rhs.f())
            return -1;
        if (lhs.f() > rhs.f())
            return 1;
        if (lhs.g > rhs.g)
            return -1;
        if (lhs.g < rhs.g)
            return 1;
        return 0;
    }

    static int hComparator(const Node& lhs, const Node& rhs) {
        if (lhs.h < rhs.h)
            return -1;
        if (lhs.h > rhs.h)
            return 1;
        return 0;
    }

    const Domain& domain;
    PriorityQueue<Node> openList{Memory::OPEN_LIST_SIZE, fComparator};
    std::unordered_map<State, Node*, typename metronome::Hasher<State>> nodes{};
    StaticVector<Node, Memory::NODE_LIMIT> nodePool;
    unsigned int iterationCounter{0};

    // Envelope Search Specific members

    Phase phase = Phase::GOAL_SEARCH;
    QueueSelector activeQueue = QueueSelector::HEURISTIC;

    Node* agentState = nullptr;
    static constexpr double pseudoGWeight = 2;

    static constexpr double greedyResourceRatio = 7.0 / 9.0;
    static constexpr double pseudoFResourceRatio = 1.0 / 9.0;
    static constexpr double explorationLimit = greedyResourceRatio + pseudoFResourceRatio;
    static constexpr double backupLimit = 1.0;

    cserna::DynamicPriorityQueue<Node*,
            PseudoFIndexFunction,
            PseudoFComparator,
            Memory::OPEN_LIST_SIZE,
            Memory::OPEN_LIST_SIZE>
            pseudoFOpenList;

    cserna::DynamicPriorityQueue<Node*,
            HeuristicIndexFunction,
            HeuristicComparator,
            Memory::OPEN_LIST_SIZE,
            Memory::OPEN_LIST_SIZE>
            heuristicOpenList;

    cserna::DynamicPriorityQueue<Node*,
            BackupIndexFunction,
            WaveFComparator,
            Memory::OPEN_LIST_SIZE,
            Memory::OPEN_LIST_SIZE>
            backupFrontier;
};
} // namespace metronome
#endif // METRONOME_ENVELOPESEARCH_HPP
