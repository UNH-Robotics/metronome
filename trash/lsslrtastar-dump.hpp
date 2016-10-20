#pragma once
#include "search.hpp"
#include "../utils/geom2d.hpp"
#include "../utils/pool.hpp"
#include <vector>
#include <queue>

template <class D>
class Lsslrtastar_dump : public SearchAlgorithm<D> {
private:

	typedef typename D::PackedState PackedState;
	typedef typename D::Operators Operators;
	typedef typename D::State State;
	typedef typename D::Cost Cost;
	typedef typename D::Oper Oper;
	typedef typename D::Edge Edge;

	class Node;
	class Nodes;

	struct inedge {
		inedge(Node *n, double c) : node(n), incost(c) {
		}

		Node *node;
		double incost;
	};

	struct outedge {
		outedge(Node *n, Oper o, Cost rc, Cost c) :
			node(n),
			op(o),
			revcost(rc == Cost(-1) ? geom2d::Infinity : rc),
			outcost(c == Cost(-1) ? geom2d::Infinity : c) {
		}

		Node *node;
		Oper op;
		double revcost;
		double outcost;
	};

	class Node {
	public:
		static ClosedEntry<Node, D> &closedentry(Node *n) {
			return n->closedent;
		}

		static PackedState &key(Node *n) {
			return n->state;
		}

		PackedState state;
		double h;
		bool expd;	// was this expanded before?
		bool goal;
		std::vector<outedge> succs;
		std::vector<inedge> preds;
		int g;

	private:
		friend class Nodes;
		friend class Pool<Node>;

		Node() : expd(false), g(-1) {
		}

		ClosedEntry<Node, D> closedent;
	};

	class Nodes {
	public:
		Nodes(unsigned int sz) : tbl(sz) {
			pool = new Pool<Node>();
		}

		~Nodes() {
			delete pool;
		}

		void clear() {
			tbl.clear();
			delete pool;
			pool = new Pool<Node>();
		}

		Node *get(D &d, State &s) {
			Node *n = pool->construct();
			d.pack(n->state, s);

			unsigned long hash = n->state.hash(&d);
			Node *found = tbl.find(n->state, hash);
			if (found) {
				pool->destruct(n);
				return found;
			}

			n->goal = d.isgoal(s);
			n->h = d.h(s);
			tbl.add(n, hash);
			return n;
		}

	private:
		ClosedList<Node, Node, D> tbl;
		Pool<Node> *pool;
	};

	class AstarNode {
	public:

		AstarNode() : openind(-1), updated(false) {
		}
	
		class Nodes {
		public:
			static ClosedEntry<AstarNode, D> &closedentry(AstarNode *n) {
				return n->nodesent;
			}
	
			static PackedState &key(AstarNode *n) {
				return n->node->state;
			}
		};
	
		class Closed {
		public:
			static ClosedEntry<AstarNode, D> &closedentry(AstarNode *n) {
				return n->closedent;
			}
	
			static PackedState &key(AstarNode *n) {
				return n->node->state;
			}
		};

		class FSort {
		public:
			static void setind(AstarNode *n, int i) {
				n->openind = i;
			}
		
			static bool pred(AstarNode *a, AstarNode *b) {	
				if (geom2d::doubleeq(a->f, b->f))
					return a->g > b->g;
				return a->f < b->f;
			}
		};

		class HSort {
		public:
			static void setind(AstarNode *n, int i) {
				n->openind = i;
			}
		
			static bool pred(AstarNode *a, AstarNode *b) {
				return a->node->h < b->node->h;
			}
		};

		Node *node;
		AstarNode *parent;
		double g, f;
		Oper op;
		long openind;
		bool updated;

	private:
		ClosedEntry<AstarNode, D> closedent, nodesent;
	};

public:

	Lsslrtastar_dump(int argc, const char *argv[]) :
		SearchAlgorithm<D>(argc, argv),
		astarClosed(1),
		astarNodes(1),
		nodes(30000001),
		lookahead(0) {

		for (int i = 0; i < argc; i++) {
			if (i < argc - 1 && strcmp(argv[i], "-lookahead") == 0)
				lookahead = strtoul(argv[++i], NULL, 10);
		}
		if (lookahead < 1)
			fatal("Must specify a lookahead ≥1 using -lookahead");

		astarPool = new Pool<AstarNode>();
		astarNodes.resize(lookahead*3);
		astarClosed.resize(lookahead*3);
	}

	~Lsslrtastar_dump() {
		delete astarPool;
	}

	void reset() {
		SearchAlgorithm<D>::reset();
		nodes.clear();
		astarOpen.clear();
		astarClosed.clear();
		astarNodes.clear();
		delete astarPool;
		astarPool = new Pool<AstarNode>();
	}

	void search(D &d, State &s0) {
		this->start();

		Node *cur = nodes.get(d, s0);

		dfrowhdr(stdout, "dumpdata", 13, "h0", "op0", "h1", "op1","h2", "op2","h3", "op3",
								"h4", "op4","h5", "op5","h6", "op6","act");

		while (!cur->goal && !this->limit()) {
			AstarNode *goal = expandLss(d, cur);
			if (!goal && !this->limit()) {
				hCostLearning(d);
			}

			lookahead += lookahead;
			if(goal && goal->node->goal)
				cur = move(cur, goal);

			times.push_back(cputime() - this->res.cpustart);
		}

		this->finish();

		if (!cur->goal) {
			this->res.ops.clear();
			return;
		}

		// Rebuild the path from the operators.
		nodes.clear();
		this->res.path.push_back(s0);
		for (auto it = this->res.ops.begin(); it != this->res.ops.end(); it++) {
			State copy = this->res.path.back();
			Edge e(d, copy, *it);
			this->res.path.push_back(e.state);
		}
		std::reverse(this->res.ops.begin(), this->res.ops.end());
		std::reverse(this->res.path.begin(), this->res.path.end());
	}

	virtual void output(FILE *out) {
		SearchAlgorithm<D>::output(out);
		dfpair(out, "num steps", "%lu", (unsigned long) times.size());
		assert (lengths.size() == times.size());
		if (times.size() != 0) {
			double min = times.front();
			double max = times.front();
			for (unsigned int i = 1; i < times.size(); i++) {
				double dt = times[i] - times[i-1];
				if (dt < min)
					min = dt;
				if (dt > max)
					max = dt;
			}
			dfpair(out, "first emit cpu time", "%f", times.front());
			dfpair(out, "min step cpu time", "%f", min);
			dfpair(out, "max step cpu time", "%f", max);
			dfpair(out, "mean step cpu time", "%f", (times.back()-times.front())/times.size());
		}
		if (lengths.size() != 0) {
			unsigned int min = lengths.front();
			unsigned int max = lengths.front();
			unsigned long sum = 0;
			for (auto l : lengths) {
				if (l < min)
					min = l;
				if (l > max)
					max = l;
				sum += l;
			}
			dfpair(out, "min step length", "%u", min);
			dfpair(out, "max step length", "%u", max);
			dfpair(out, "mean step length", "%g", sum / (double) lengths.size());
		}
	}

private:

	// ExpandLss returns the cheapest  goal node if one was generated
	// and NULL otherwise.
	AstarNode *expandLss(D &d, Node *rootNode) {
		astarOpen.clear();
		astarNodes.clear();
		astarClosed.clear();
		astarPool->releaseall();

		AstarNode *a = astarPool->construct();
		a->node = rootNode;
		a->parent = NULL;
		a->op = D::Nop;
		a->g = 0;
		a->f = rootNode->h;
		astarOpen.push(a);
		astarNodes.add(a);

		AstarNode *goal = NULL;

		for (unsigned int exp = 0; !astarOpen.empty() && exp < lookahead && !this->limit(); exp++) {
			AstarNode *s = *astarOpen.pop();

			if (!astarClosed.find(s->node->state))
				astarClosed.add(s);

			AstarNode *g = expandPropagate(d, s, true);

			if (g && (!goal || g->g < goal->g))
				goal = g;

			if (s->node->goal) {
				astarOpen.push(s);
				goal = s;
				break;
			}

			if(exp % 100 == 0) {
				dump(rootNode, false);
			}

		}

		dump(rootNode, true);

		return goal;
	}

	void dump(Node* s, bool act) {
		fprintf(stdout, "#altrow  \"dumpdata\"");
		fflush(stdout);

		int i = 0;
		for(auto succ : s->succs) {

			double h = findMinF(succ.node, 1);
			fprintf(stdout, "\t%g\t%d", h, (int)succ.op);
			fflush(stdout);

			i++;
		}
		for( ; i < 6; i++)
			fprintf(stdout, "\t-1\t-1");

		fprintf(stdout, "\t%d\n", act);
		fflush(stdout);
	}

	double findMinF(Node* s, unsigned int g) {
		for(AstarNode* n : astarNodes)
			n->node-> g = -1;

		std::queue<Node*> q;
		q.push(s);

		s->g = 0;

		double minF = geom2d::Infinity;

		while(q.size() > 0) {
			Node* n = q.front();
			q.pop();

			int nSuccs = 0;

			for(auto succ: n->succs) {
				if(!astarNodes.find(succ.node->state) || succ.node->g > 0) continue;

				succ.node->g = n->g + 1;

				q.push(succ.node);
				nSuccs++;
			}

			if(nSuccs == 0 && (n->g + n->h) < minF) //leaf
				minF = (n->g + n->h);

		}
		return minF;
	}

	// ExpandPropagate returns the cheapest goal node if one is generated
	// and NULL otherwise.
	AstarNode *expandPropagate(D &d, AstarNode *s, bool doExpand) {
		AstarNode *goal = NULL;

		if (this->limit())
			return NULL;

		auto kids = expand(d, s->node);
		for (unsigned int i = 0; i < kids.size(); i++) {
			auto e = kids[i];
			Node *k = e.node;
			AstarNode *kid = astarNodes.find(k->state);

			if (doExpand && astarNodes.find(s->node->state)) {
				if (!kid) {
					kid = astarPool->construct();
					kid->node = k;
					kid->g = geom2d::Infinity;
					kid->openind = -1;
					astarNodes.add(kid);
				}
				if (kid->g > s->g + e.outcost) {
					if (kid->parent)	// !NULL if a dup
						this->res.dups++;
					kid->parent = s;
					kid->g = s->g + e.outcost;
					kid->f = kid->g + kid->node->h;
					kid->op = e.op;
					astarOpen.pushupdate(kid, kid->openind);
				}
			}

			if (k->goal && kid && (!goal || kid->g < goal->g))
				goal = kid;
		}

		return goal;
	}

	void hCostLearning(D &d) {
		BinHeap<typename AstarNode::HSort, AstarNode*> open;

		open.append(astarOpen.data());

		unsigned long left = astarClosed.getFill();

		while (left > 0 && !open.empty()) {
			AstarNode *s = *open.pop();

			if (astarClosed.find(s->node->state))
				left--;

			for (auto e : s->node->preds) {
				Node *sprime = e.node;

				AstarNode *sp = astarClosed.find(sprime->state);
				if (!sp)
					continue;

				if (!sp->updated || sprime->h > e.incost + s->node->h) {
					sprime->h = e.incost + s->node->h;
					sp->updated = true;
					open.pushupdate(sp, sp->openind);
				}
			}
		}
	}

	Node *move(Node *cur, AstarNode *goal) {
		AstarNode *best = goal;
		if (!best) {
			for (auto n : astarOpen.data()) {
				if (n->node == cur)
					continue;
				if (best == NULL || AstarNode::FSort::pred(n, best))
					best = n;
			}
		}
	
		assert (best->node != cur);
		std::vector<Oper> ops;
//		Node* prev;
		for (AstarNode *p = best; p->node != cur; p = p->parent) {
			assert (p->parent != best);	// no cycles
			ops.push_back(p->op);
//			prev = p->node;
		}

//		this->res.ops.push_back(ops.back());
//		lengths.push_back(1);
//		return prev;

		this->res.ops.insert(this->res.ops.end(), ops.rbegin(), ops.rend());
		lengths.push_back(ops.size());
		return best->node;

	}

	// Expand returns the successor nodes of a state.
	std::vector<outedge> expand(D &d, Node *n) {
		if (n->expd)
			return n->succs;

		State buf, &s = d.unpack(buf, n->state);

		this->res.expd++;

		Operators ops(d, s);
		for (unsigned int i = 0; i < ops.size(); i++) {
			this->res.gend++;

			Edge e(d, s, ops[i]);
			Node *k = nodes.get(d, e.state);
			k->preds.emplace_back(n, e.cost);
			n->succs.emplace_back(k, ops[i], e.revcost, e.cost);
		}
		n->expd = true;

		return n->succs;
	}

	BinHeap<typename AstarNode::FSort, AstarNode*> astarOpen;
	ClosedList<typename AstarNode::Closed, AstarNode, D> astarClosed;
	ClosedList<typename AstarNode::Nodes, AstarNode, D> astarNodes;
	Pool<AstarNode> *astarPool;

	Nodes nodes;
	unsigned int lookahead;

	std::vector<double> times;
	std::vector<unsigned int> lengths;
};