#ifndef METRONOME_ASTAR_HPP
#define METRONOME_ASTAR_HPP
#define BOOST_POOL_NO_MT

#include "util/Hasher.hpp"
#include <boost/pool/object_pool.hpp>
#include <unordered_map>
#include <util/PriorityQueue.hpp>
#include <vector>

namespace metronome {

template <typename Domain>
class AStar {
    typedef typename Domain::State State;
    typedef typename Domain::Action Action;
    typedef typename Domain::Cost Cost;

public:
    AStar(Domain domain) : domain(domain), openList(10000000, fValueComparator) {
    }

    std::vector<Action> plan(State startState) {
        std::vector<Action> constructedPlan;

        const Node localStartNode =
                Node(nullptr, std::move(startState), Action(), 0, domain.heuristic(startState), true);
        auto startNode = nodePool.construct(localStartNode);

        nodes[localStartNode.state] = startNode;

        openList.push(localStartNode);

        //        while (!openList.isEmpty()) {
        //            Node* q = openList.pop();
        //
        //            constructedPlan.push_back(q->action);
        //
        //            std::vector<State> successors = domain.successors(q->state);
        //
        //            for (State successor : successors) {
        //                if (domain.isGoal(successor)) {
        //                    openList.clear();
        //                    return constructedPlan;
        //                }
        //
        //                Node n = Node(q, successor, successor.getAction(), q->g + successor.getCost(),
        //                        q->g + successor.getCost() + domain.heuristic(successor), true);
        //
        //                auto node = nodePool.construct(n);
        //
        //                hasher(n.state);
        //
        //                nodes[n.state] = node;
        //                openList.push(n);
        //            }
        //        }
        //
        //        return constructedPlan;

        while (!openList.isEmpty()) {
            // TODO increment expanded counter
            Node* currentNode = openList.pop();

            if (!currentNode->open) {
                continue; // This node was disabled
            }

            if (domain.isGoal(currentNode->state)) {
                // Goal is reached TODO extract plan
            }

            for (auto successor : domain.getSuccessors(currentNode->state)) {
                if (successor.state == currentNode->state) {
                    continue; // Skip parent TODO this might be unnecessary
                }

                // TODO increment generated node count

                auto& existingSuccessorNode = nodes[successor.state];
                auto newCost = successor.actionCost + currentNode->g;

                if (existingSuccessorNode == nullptr) {
                    const Node successorNode(currentNode, successor.state);
                    existingSuccessorNode =
                }
            }
        }

        return std::vector<Action>();
    }

private:
    class Node {
    public:
        Node(Node* parent, State state, Action action, Cost g, Cost f, bool open)
                : parent(parent), state(state), action(std::move(action)), g(g), f(f), open(open) {
        }

        unsigned long hash() {
            return state->hash();
        }

        bool operator==(const Node& node) const {
            return state == node.state;
        }

        mutable unsigned int index;
        Node* parent;
        State state;
        Action action;
        Cost g;
        Cost f;
        /** True if the node is in the open list. */
        bool open;
    };

    static int fValueComparator(const Node& lhs, const Node& rhs) {
        if (lhs.f < rhs.f)
            return -1;
        if (lhs.f > rhs.f)
            return 1;
        if (lhs.f > rhs.f)
            return -1;
        if (lhs.f < rhs.f)
            return 1;
        return 0;
    }

    Domain domain;
    PriorityQueue<Node> openList;
    std::unordered_map<State, Node*, typename metronome::Hasher<State>> nodes;
    boost::object_pool<Node> nodePool{100000000, 100000000};
};
}

#endif // METRONOME_ASTAR_HPP
