#include "search.hpp"
#include "../utils/pool.hpp"

void fatal(const char*, ...);	// utils.hpp

template <class D> struct Wastar : public SearchAlgorithm<D> {

	typedef typename D::State State;
	typedef typename D::PackedState PackedState;
	typedef typename D::Cost Cost;
	typedef typename D::Oper Oper;

	struct Node {

		ClosedEntry<Node, D> closedent;
		int openind;
		Node *parent;
		PackedState state;
		Oper op, pop;
		Cost f, g;
		double fprime;

		Node() : openind(-1) {
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
			if (a->fprime == b->fprime) {
				if (a->f == b->f)
					return a->g > b->g;
				return a->f < b->f;
			}
			return a->fprime < b->fprime;
		}	
	};

	Wastar(int argc, const char *argv[]) :
			SearchAlgorithm<D>(argc, argv), dropdups(false),
			wt(-1.0), closed(30000001) {
		for (int i = 0; i < argc; i++) {
			if (i < argc - 1 && strcmp(argv[i], "-wt") == 0)
				wt = strtod(argv[++i], NULL);
			if (strcmp(argv[i], "-dropdups") == 0)
				dropdups = true;
		}

		if (wt < 1)
			fatal("Must specify a weight ≥1 weight using -wt");

		nodes = new Pool<Node>();
	}

	~Wastar() {
		delete nodes;
	}

	void search(D &d, typename D::State &s0) {
		this->start();
		closed.init(d);

		Node *n0 = init(d, s0);
		closed.add(n0);
		open.push(n0);

		while (!open.empty() && !SearchAlgorithm<D>::limit()) {
			Node *n = *open.pop();
			State buf, &state = d.unpack(buf, n->state);

			if (d.isgoal(state)) {
				solpath<D, Node>(d, n, this->res);
				break;
			}

			expand(d, n, state);
		}
		this->finish();
	}

	virtual void reset() {
		SearchAlgorithm<D>::reset();
		open.clear();
		closed.clear();
		delete nodes;
		nodes = new Pool<Node>();
	}

	virtual void output(FILE *out) {
		SearchAlgorithm<D>::output(out);
		closed.prstats(stdout, "closed ");
		dfpair(stdout, "open list type", "%s", "binary heap");
		dfpair(stdout, "node size", "%u", sizeof(Node));
		dfpair(stdout, "weight", "%lg", wt);
	}

private:

	void expand(D &d, Node *n, State &state) {
		SearchAlgorithm<D>::res.expd++;

		typename D::Operators ops(d, state);
		for (unsigned int i = 0; i < ops.size(); i++) {
			if (ops[i] == n->pop)
				continue;
			SearchAlgorithm<D>::res.gend++;
			considerkid(d, n, state, ops[i]);
		}
	}

	void considerkid(D &d, Node *parent, State &state, Oper op) {
		Node *kid = nodes->construct();
		typename D::Edge e(d, state, op);
		kid->g = parent->g + e.cost;
		d.pack(kid->state, e.state);

		unsigned long hash = kid->state.hash(&d);
		Node *dup = static_cast<Node*>(closed.find(kid->state, hash));
		if (dup) {
			this->res.dups++;
			if (dropdups || kid->g >= dup->g) {
				nodes->destruct(kid);
				return;
			}
			if (dup->openind < 0)
				this->res.reopnd++;
			dup->fprime = dup->fprime - dup->g + kid->g;
			dup->f = dup->f - dup->g + kid->g;
			dup->g = kid->g;
			dup->parent = parent;
			dup->op = op;
			dup->pop = e.revop;
			open.pushupdate(dup, dup->openind);
			nodes->destruct(kid);
		} else {
			typename D::Cost h = d.h(e.state);
			kid->f = kid->g + h;
			kid->fprime = kid->g + wt * h;
			kid->parent = parent;
			kid->op = op;
			kid->pop = e.revop;
			closed.add(kid, hash);
			open.push(kid);
		}
	}

	Node *init(D &d, State &s0) {
		Node *n0 = nodes->construct();
		d.pack(n0->state, s0);
		n0->g = Cost(0);
		typename D::Cost h = d.h(s0);
		n0->f = h;
		n0->fprime = wt * h;
		n0->op = n0->pop = D::Nop;
		n0->parent = NULL;
		return n0;
	}

	bool dropdups;
	double wt;
	BinHeap<Node, Node*> open;
 	ClosedList<Node, Node, D> closed;
	Pool<Node> *nodes;
};
