#pragma once
#include "search.hpp"
#include "../utils/pool.hpp"

template <class D> struct Frit_Offline : public SearchAlgorithm<D> {

	typedef typename D::State State;
	typedef typename D::PackedState PackedState;
	typedef typename D::Cost Cost;
	typedef typename D::Oper Oper;
	typedef typename D::Operators Operators;

	struct Node {
		ClosedEntry<Node, D> closedent;
		int openind;
		Node *parent;
		PackedState state;
		Oper op, pop;
		Cost h;
		int bfs_d;
		Node* idealTreeParent;
		unsigned int color;

		Node() : openind(-1), bfs_d(-1), idealTreeParent(NULL), color(0) {
		}

		static ClosedEntry<Node, D> &closedentry(Node *n) {
			return n->closedent;
		}

		static PackedState &key(Node *n) {
			return n->state;
		}

		static void setind(Node *n, int i) {
			n->openind = i;
		}

		static int getind(const Node *n) {
			return n->openind;
		}

		static bool pred(Node *a, Node *b) {
			if(a->bfs_d == b->bfs_d)
				return a->h < b->h;
			return a->bfs_d < b->bfs_d;
		}
	};

	Frit_Offline(int argc, const char *argv[]) :
		SearchAlgorithm<D>(argc, argv), all(30000001), seen(30000001), currentColor(0) {
		nodes = new Pool<Node>();
	}

	~Frit_Offline() {
		delete nodes;
	}

	void search(D &d, typename D::State &s0) {
		this->start();

		seen.init(d);
		Node *n0 = init(d, s0);
		Node *current = n0;
		State buf, &state = d.unpack(buf, current->state);

		while(!d.isgoal(state)) {
			std::vector<Node*> path = doBFS(d,current,state);

			//I used to loop over all the nodes returned from BFS
			//but sometimes you get impassable nodes because of the
			//free space assumption. Reevaluate at each step to be safe.

			current->idealTreeParent = path.front();

			current = path.front();
			this->res.ops.push_back(current->op);

			unsigned long hash = current->state.hash(&d);

			if(!seen.find(current->state,hash))
				seen.add(current, hash);

			state = d.unpack(buf, current->state);
		}
		this->finish();

		if (!d.isgoal(state)) {
			this->res.ops.clear();
			return;
		}

		// Rebuild the path from the operators to avoid storing very long
		// paths as we go.
		seen.clear();
		this->res.path.push_back(s0);

		for (auto it = this->res.ops.begin(); it != this->res.ops.end(); it++) {
			State copy = this->res.path.back();
			typename D::Edge e(d, copy, *it);
			this->res.path.push_back(e.state);
		}
		std::reverse(this->res.ops.begin(), this->res.ops.end());
		std::reverse(this->res.path.begin(), this->res.path.end());
	}

	virtual void reset() {
		SearchAlgorithm<D>::reset();
		all.clear();
		seen.clear();
		delete nodes;
		nodes = new Pool<Node>();
	}

	virtual void output(FILE *out) {
		SearchAlgorithm<D>::output(out);
		seen.prstats(stdout, "seen ");
		dfpair(stdout, "node size", "%u", sizeof(Node));
	}

private:

	std::vector<Node*> doBFS(D &d, Node *n, const State &s) {
		std::vector<Node*> path;
		OpenList<Node, Node, Cost> open;
	 	ClosedList<Node, Node, D> closed(1000);

		Node* current = n;
		State buf, &state = d.unpack(buf, current->state);

		assert(seen.find(current->state));

		currentColor++;

		current->bfs_d = 0;
		closed.add(current);
		expand(d, current, state, open, closed);

		bool foundGoal = false;

		while(!open.empty()) {
			current = open.pop();
			state = d.unpack(buf, current->state);

			if(inTreeGoalTest(d, current, currentColor)) {
				foundGoal = true;
				break;
			}
			
			if(seen.find(current->state))
				expand(d, current, state, open, closed);
			else
				naiveExpand(d, current, state, open, closed); 
		}

		if(foundGoal) {
			while(current != n) {
				path.push_back(current);
				current = current->parent;
			}
			std::reverse(path.begin(), path.end());
		}

		return path;
	}

	bool inTreeGoalTest(D &d, Node*n, unsigned int color) {

		Node* current = n;
		State buf, &state = d.unpack(buf, n->state);

		while(!d.isgoal(state)) {
			current->color = color;

			Node* next = followTree(d, current, state);

			current->idealTreeParent = next;

			if(next->color == color)
				return false;

			current = next;
			state = d.unpack(buf, current->state);
		}
		return true;
	}

	void expand(const D &d, Node *n, const State &s,
				OpenList<Node, Node, Cost>& open, ClosedList<Node, Node, D>& closed) {
		SearchAlgorithm<D>::res.expd++;

		Operators ops(d,s);

		for (unsigned int i = 0; i < ops.size(); i++) {
			SearchAlgorithm<D>::res.gend++;
			considerkid(d, n, s, ops[i], open, closed);
		}
	}

	void naiveExpand(const D &d, Node *n, const State &s,
				OpenList<Node, Node, Cost>& open, ClosedList<Node, Node, D>& closed) {
		SearchAlgorithm<D>::res.expd++;

		Operators inboundsOps = Operators::inboundsOps(d,s);

		for (unsigned int i = 0; i < inboundsOps.size(); i++) {
			SearchAlgorithm<D>::res.gend++;
			considerkid(d, n, s, inboundsOps[i], open, closed);
		}
	}

	void considerkid(const D &d, Node *parent, const State &s, const Oper op,
				OpenList<Node, Node, Cost>& open, ClosedList<Node, Node, D>& closed) {
		Node *kid = nodes->construct();
		assert (kid);
		typename D::Edge e(d, s, op);
		d.pack(kid->state, e.state);

		unsigned long hash = kid->state.hash(&d);

		Node *allDup = all.find(kid->state, hash);
		if(allDup) {
			nodes->destruct(kid);
			kid = allDup;
		}
		else {
			all.add(kid, hash);
		}

		Node *dup = closed.find(kid->state, hash);

		if (dup) {
			if(seen.find(parent->state) && kid->bfs_d > parent->bfs_d + 1) {
				kid->op = op;
				kid->pop = e.revop;
			}
			return;
		}
		else {
			kid->bfs_d = parent->bfs_d + 1;
			kid->parent = parent;
			kid->op = op;
			kid->h = d.h(e.state);
			kid->pop = e.revop;
			closed.add(kid, hash);
			open.push(kid);
		}
	}

	Node* followTree(D &d, Node *n, State &s) {
		SearchAlgorithm<D>::res.expd++;

		Node* minKid = NULL;
		Cost minCost;

		Operators ops = seen.find(n->state) ? Operators(d, s) : Operators::inboundsOps(d, s);

		assert(ops.size() > 0);

		for (unsigned int i = 0; i < ops.size(); i++) {
			SearchAlgorithm<D>::res.gend++;

			std::pair<Node*, Cost> nodeCostPair = generateSloppyKid(d, n, s, ops[i]);
			Node* kid = nodeCostPair.first;
			Cost kidValue = kid->h + nodeCostPair.second;

			if(minKid == NULL || minCost > kidValue) {
				minKid = kid;
				minCost = kidValue;
			}
		}

		return minKid;
	}

	std::pair<Node*, Cost> generateSloppyKid(const D &d, Node *parent, const State &s, const Oper op) {
		Node *kid = nodes->construct();
		assert (kid);
		typename D::Edge e(d, s, op);
		d.pack(kid->state, e.state);

		unsigned long hash = kid->state.hash(&d);

		Node *allDup = all.find(kid->state, hash);
		if(allDup) {
			nodes->destruct(kid);
			kid = allDup;
		}
		else {
			kid->bfs_d = 0;
			kid->parent = parent;
			kid->op = op;
			kid->h = d.h(e.state);
			kid->pop = e.revop;
			all.add(kid, hash);
		}

		return std::pair<Node*, Cost>(kid, e.cost);
	}

	Node *init(const D &d, State &s0) {
		Node *n0 = nodes->construct();
		d.pack(n0->state, s0);
		n0 ->h = d.h(s0);
		n0->pop = n0->op = D::Nop;
		n0->parent = NULL;
		seen.add(n0);
		all.add(n0);
		return n0;
	}

	ClosedList<Node, Node, D> all;
 	ClosedList<Node, Node, D> seen;
	Pool<Node> *nodes;
	unsigned int currentColor;
};
