#include "search.hpp"
#include "../utils/utils.hpp"
#include <algorithm>
#include <forward_list>

template <class D> struct Rtastar : public SearchAlgorithm<D> {
	typedef typename D::State State;
	typedef typename D::PackedState PackedState;
	typedef typename D::Cost Cost;
	typedef typename D::Oper Oper;

	struct Node {
		ClosedEntry<Node, D> closedent;
		PackedState packed;
		Cost h;

		static PackedState &key(Node *n) { return n->packed; }

		static ClosedEntry<Node, D> &closedentry(Node *n) { return n->closedent; }
	};

	Rtastar(int argc, const char *argv[]) :
			SearchAlgorithm<D>(argc, argv),
			seen(30000001) {
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "-learning") == 0)
			    learning = true;
		}
		lsslim = LookaheadLimit::fromArgs(argc, argv);
		nodes = new Pool<Node>();
	}

	~Rtastar() {
		delete nodes;
	}

	struct Current {
		Current(const State &s, Oper o, Oper p, Cost e, Cost fval) :
			state(s),
			op(o),
			pop(p),
			edgecost(e),
			f(fval) { }

		State state;
		Oper op, pop;
		Cost edgecost;
		Cost f;
	};

	void search(D &d, typename D::State &s0) {
		this->start();
		seen.init(d);
		Current cur(s0, D::Nop, D::Nop, typename D::Cost(0), d.h(s0));
		Cost *curh = storenode(d, cur.state, cur.f);

		while (!d.isgoal(cur.state) && !this->limit()) {
			std::vector<Current> bests = bestkids(d, cur, *curh);
			if (bests.size() == 0) {
				this->res.ops.clear();
				break;	// deadend;
			}
			storenode(d, cur.state, *curh);
			//std::random_shuffle(bests.begin(),bests.end());
			cur = bests.front();
			this->res.ops.push_back(cur.op);
		}

		this->finish();

		if (!d.isgoal(cur.state)) {
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
		seen.clear();
		delete nodes;
		nodes = new Pool<Node>();
	}

	virtual void output(FILE *out) {
		SearchAlgorithm<D>::output(out);
		// We output the 'closed ' prefix here so that they key
		// matches the closed keys for other search algorithms.
		// This is desriable because the seen list performs a
		// similar functionality as a closed list, and we will
		// probably want to compare with closed list stats.
		seen.prstats(stdout, "closed ");
		dfpair(stdout, "node size", "%u", sizeof(Node));
	}

private:

	// bestkids returns a vector of all of the successors that have the
	// best f value.
	std::vector<Current> bestkids(D &d, Current &cur, Cost &sndf) {
		sndf = Cost(-1);
		std::vector<Current> bests;

		this->res.expd++;
		typename D::Operators ops(d, cur.state);
		unsigned int left = lsslim->lookahead( cur.edgecost );
		unsigned int left_prev = 0;
		for (unsigned int n = 0; n < ops.size(); n++) {

			this->res.gend++;
			typename D::Edge e(d, cur.state, ops[n]);
			Cost h = heuristic(d, e.state, e.revop, left * ( n + 1 ) / ops.size( ) - left_prev );
			left_prev = left * ( n + 1 ) / ops.size( );
			Cost f = ((h == Cost(-1)) ? h : h + e.cost);

			if (bests.empty() || better(f, bests[0].f)) {
				if (!bests.empty())
					sndf = bests[0].f;
				bests.clear();
				bests.push_back(Current(e.state, ops[n], e.revop, e.cost, f));
			} else if (bests[0].f == f) {
				bests.push_back(Current(e.state, ops[n], e.revop, e.cost, f));
			} else if (better(f, sndf)) {
				sndf = f;
			}
		}

		if( learning )
		    sndf = bests[0].f;

		return bests;
	}

	// value returns the f value for the given state
	// either from the seen table or by performing
	// lookahead search.
	Cost heuristic(D &d, State &cur, Oper pop, unsigned int left) {
		Node *n = nodes->construct();
		d.pack(n->packed, cur);
		unsigned long hash = n->packed.hash(&d);

		Node *dup = seen.find(n->packed, hash);
		if (dup) {
			nodes->destruct(n);
			return dup->h;
		}
		Cost alpha = Cost(-1);
		nodes->destruct(n);
		return look(d, cur, alpha, pop, Cost(0), left);
	}

    // look returns the lookahead cost via a depth-limited,
    // alpha pruned, depth-first search.
    Cost look(D &d, State &state, Cost &alpha, Oper pop, Cost g, unsigned int left) {
        Cost f = d.h(state) + g;
        if (this->limit() || left == 0 || !better(f, alpha) || d.isgoal(state)) {
            if (better(f, alpha))
                alpha = f;
            return f;
        }

        this->res.expd++;
        Cost bestf = Cost(-1);
        typename D::Operators ops(d, state);
        for (unsigned int n = 0; n < ops.size(); n++) {
            if (ops[n] == pop)
                continue;

            this->res.gend++;
            typename D::Edge e(d, state, ops[n]);
            f = look(d, e.state, alpha, e.revop, g + e.cost, left-1);
            if (better(f, bestf))
                bestf = f;
        }

        return bestf;
    }

	// storenode marks the node as seen if it is not already
	// in the seen list and returns a pointer to the stored
	// heuristic value.  This pointer can be used to update
	// the value of this node without looking it up again.
	Cost *storenode(D &d, State &cur, Cost h) {
		Node *n = nodes->construct();
		d.pack(n->packed, cur);
		unsigned long hash = n->packed.hash(&d);

		Node *dup = seen.find(n->packed, hash);
		if (dup) {
			nodes->destruct(n);
			dup->h = h;
			return &dup->h;
		}

		n->h = h;
		seen.add(n, hash);
		return &n->h;
	}

	// better returns true if a is better than b.
	static bool better(Cost a, Cost b) {
		return b == Cost(-1) || (a != Cost(-1) && a < b);
	}

	bool learning;
	LookaheadLimit *lsslim;
	ClosedList<Node, Node, D> seen;
	Pool<Node> *nodes;
};
