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

		bool isDeadend() const {
			return this == idealTreeParent;
		}
		
		void setDeadend() { 
			idealTreeParent = this;
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

//		Operators allOps = Operators(d);

		seen.init(d);
		Node *n0 = init(d, s0);
		Node *current = n0;
		seen.add(current);
		all.add(current);

		State buf, &state = d.unpack(buf, current->state);
int cur = 0;
		while(!d.isgoal(state)) {
			std::vector<Node*> path;
fprintf(stderr, "%d: %d\n", cur++, state.getLoc()); d.dumpstate(stderr, state);
			path.push_back(getBestNeighbor(d, current, state));

			currentColor++;

			if(path[0] == NULL || !inTreeGoalTest(d, path[0], currentColor)) {
fprintf(stderr, "no best neighbor\n");
				displayCurrentTree(d, state);

				current->setDeadend();
				path = doBFS(d, current, state);

				if(path.size() == 0) {
					fprintf(stderr, "NO SOLUTION EXISTS!\n");
					break;
				}
			}
			else {
				fprintf(stderr, "GOT  best neighbor\n");
			}

			current->idealTreeParent = path[0];

			for(unsigned int i = 0; i < path.size(); i++) {
				state = d.unpack(buf, current->state);

//				if(!allOps.applicable(d, state, path[i]->op)) break;

				displayCurrentTree(d, state);

				current = path[i];

				this->res.ops.push_back(path[i]->op);
				
				if(i < path.size() - 1)
					path[i]->idealTreeParent = path[i+1];

				unsigned long hash = path[i]->state.hash(&d);
				if(!seen.find(path[i]->state,hash))
					seen.add(path[i], hash);
			}

			state = d.unpack(buf, current->state);
if(cur > 10) break;
		}

		this->finish();

		if (!d.isgoal(state)) {
			this->res.ops.clear();
			return;
		}
			displayCurrentTree(d, state);

		// Rebuild the path from the operators to avoid storing very long
		// paths as we go.
		seen.clear();
		this->res.path.push_back(s0);

		for (auto it = this->res.ops.begin(); it != this->res.ops.end(); it++) {
			State copy = this->res.path.back();
			typename D::Edge e(d, copy, *it);
			this->res.path.push_back(e.state);
			assert(d.map->ok(copy.getLoc(), d.map->mvs[*it]));
		}
		std::reverse(this->res.ops.begin(), this->res.ops.end());
		std::reverse(this->res.path.begin(), this->res.path.end());
	}

	virtual void reset() {
		SearchAlgorithm<D>::reset();
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
/* Local Examination

    These functions are meant to maintain feasibility of moves and update the seen list */

	/* Get the best neighbor that can actually be moved to */

	Node* getBestNeighbor(D &d, Node *n, const State &s) {
		SearchAlgorithm<D>::res.expd++;

		Node* minKid = NULL;
		Cost minCost;

		Operators inboundsOps(d,s,true);

		for (unsigned int i = 0; i < inboundsOps.size(); i++) {

			SearchAlgorithm<D>::res.gend++;

			std::pair<Node*, Cost> nodeCostPair = generateKid(d, n, s, inboundsOps[i]);
			Node* kid = nodeCostPair.first;
			Cost kidValue = kid->h + nodeCostPair.second;

			if(minKid == NULL || minCost >= kidValue) {
				minKid = kid;
				minCost = kidValue;
			}
		}

		return inboundsOps.applicable(d, s, minKid->op) ? minKid : NULL;
	}

	/* The agent is in a state that can "see" all these children so update seen list */
	std::pair<Node*, Cost>generateKid(const D &d, Node *parent, const State &state, const Oper op) {
		Node *kid = nodes->construct();
		assert(kid);
		typename D::Edge e(d, state, op);
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

		Node *dup = seen.find(kid->state, hash);
		if (dup) {
			dup->parent = parent;
			dup->op = op;
			dup->pop = e.revop;
			return std::pair<Node*, Cost>(dup, e.cost);
		} else {
			kid->h = d.h(e.state);
			kid->parent = parent;
			kid->op = op;
			kid->pop = e.revop;
			seen.add(kid, hash);
			return std::pair<Node*, Cost>(kid, e.cost);
		}
	}

/* END LOCAL EXAMINATION */

/* BFS STUFF

    These functions are meant to operate on the hypothetical tree and the pieces of the space
    that have been seen previously */

	/* Doing BFS but not updating the global DS because the agent isn't actually 
		experiencing these states */

	std::vector<Node*> doBFS(D &d, Node *n, const State &s) {
		std::vector<Node*> path;
		OpenList<Node, Node, Cost> open;
	 	ClosedList<Node, Node, D> closed(1000);

		Node* current = n;
		State buf, &state = d.unpack(buf, current->state);

		if(state.getLoc() == 22)
			fprintf(stderr, "GOTCHA\n");

		currentColor++;

		current->bfs_d = 0;
		open.push(current);
		closed.add(current);
		bool foundGoal = false;
		while(!open.empty()) {
			current = open.pop();
			state = d.unpack(buf, current->state);

			if(inTreeGoalTest(d, current, currentColor)) {
				foundGoal = true;
				break;
			}
			expand(d, current, state, open, closed);
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

	/* checking to see if the state is actually going to lead to the goal in the imaginary tree
	     and pieces of the tree we've already seen and know about */

	bool inTreeGoalTest(D &d, Node*n, unsigned int color) {
		Node* current = n;
		State buf, &state = d.unpack(buf, n->state);

		if(current->isDeadend()) return false;

		while(!d.isgoal(state)) {

			current->color = color;

			Node* next = current->idealTreeParent ? current->idealTreeParent : followTree(d, current, state);

			if(next->isDeadend() || next->color == color)
				return false;

			current = next;
			state = d.unpack(buf, current->state);
		}
		return true;
	}

	/* According to the tree, that we may not have seen states from, where should the agent go
	    (basically according to the heuristic */

	Node* followTree(D &d, Node *n, State &s) {
		SearchAlgorithm<D>::res.expd++;

		Node* minKid = NULL;
		Cost minCost;

		Operators inboundsOps(d,s,true);
		for (unsigned int i = 0; i < inboundsOps.size(); i++) {
			SearchAlgorithm<D>::res.gend++;

			std::pair<Node*, Cost> nodeCostPair = generateSloppyKid(d, n, s, inboundsOps[i]);
			Node* kid = nodeCostPair.first;
			if(kid == NULL) continue;
			Cost kidValue = kid->h + nodeCostPair.second;

			if(minKid == NULL || minCost > kidValue) {
				/* we might be able to avoid a single expansion BFS if we check
					ties on h in the case of a minHKid that is has already been
					removed from the ideal tree */
				minKid = kid;
				minCost = kidValue;
			}
		}

		return minKid;
	}

	std::pair<Node*, Cost>generateSloppyKid(const D &d, const Node *parent, 
							const State &state, const Oper op) {
		Node *kid = nodes->construct();
		assert(kid);

		typename D::Edge e(d, state, op);
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

		Node *dup = seen.find(kid->state, hash);
		if (dup) {
			if(d.map->ok(state.getLoc(), d.map->mvs[op]))
				return std::pair<Node*, Cost>(dup, e.cost);
			else
				return std::pair<Node*, Cost>(NULL, e.cost);
		} else {
			kid->h = d.h(e.state);
			kid->op = op;
			kid->pop = e.revop;
			return std::pair<Node*, Cost>(kid, e.cost);
		}
	}

	/* expand this state but don't update global DS */

	void expand(const D &d, Node *n, const State &s,
				OpenList<Node, Node, Cost>& open, ClosedList<Node, Node, D>& closed) {
		SearchAlgorithm<D>::res.expd++;

		Operators inboundsOps = Operators(d,s,true);

		for (unsigned int i = 0; i < inboundsOps.size(); i++) {
			SearchAlgorithm<D>::res.gend++;
			considerkid(d, n, s, inboundsOps[i], inboundsOps, open, closed);
		}
	}

	/* look at this kid and decide if you should keep it */

	void considerkid(const D &d, Node *parent, const State &s, const Oper op, const Operators &ops,
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

		Node *seenDup  = seen.find(kid->state, hash);

		if (dup || (seenDup && !ops.applicable(d,s,op))) {
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

/* DONE BFS STUFF */

	/* this really has nothing to do with anything except keeping the code looking decent */

	Node *init(const D &d, State &s0) {
		Node *n0 = nodes->construct();
		d.pack(n0->state, s0);
		n0 ->h = d.h(s0);
		n0->pop = n0->op = D::Nop;
		n0->parent = NULL;
		return n0;
	}

	void displayCurrentTree(const D &d, const State &s) {
		unsigned int loc = s.getLoc();
		unsigned int goal = d.finish;
		for(unsigned int y = 0; y < d.map->h; y++) {
			for(unsigned int x = 0; x < d.map->w; x++) {
				unsigned int cell = d.map->index(x,y);
				
				if(d.map->blkd(cell)) {
					fprintf(stderr, "#");
				}
				else if(cell == loc) {
					fprintf(stderr, "@");
				}
				else if(cell == goal) {
					fprintf(stderr, "*");
				}
				else {
					PackedState ps;
					ps.loc = cell;
					Node* n  = seen.find(ps);
					if(n && n->idealTreeParent) {
						unsigned int dir = n->idealTreeParent->state.loc;
						std::pair<int,int> xy = d.map->coord(dir);
						if(xy.first > (int)x) fprintf(stderr, ">");
						else if(xy.first < (int)x) fprintf(stderr, "<");
						else if(xy.second < (int)y) fprintf(stderr, "↑");
						else if(xy.second > (int)y) fprintf(stderr, "↓");
						else fprintf(stderr, "?");
					}
					else {
						fprintf(stderr, " ");
					}
				}
			}
			fprintf(stderr, "\n");
		}


//		fprintf(stderr, "waiting for keypress"); std::cin.ignore();
	}

 	ClosedList<Node, Node, D> all;
 	ClosedList<Node, Node, D> seen;
	Pool<Node> *nodes;
	unsigned int currentColor;
};
