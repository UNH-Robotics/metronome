#pragma once
#include "search.hpp"
#include "../utils/geom2d.hpp"
#include "../utils/pool.hpp"
#include "lsslrtastar2.hpp"
#include <vector>

template <class D>
class GStar : public SearchAlgorithm<D> {
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
		double h, horig, d;
		bool expd;	// was this expanded before?
		bool goal;
		std::vector<outedge> succs;
		std::vector<inedge> preds;

	private:
		friend class Nodes;
		friend class Pool<Node>;

		Node() : expd(false) {
		}

		ClosedEntry<Node, D> closedent;
	};

	class Nodes {
	public:
		Nodes(unsigned int sz) : tbl(sz) {
		}

		~Nodes() {
		}

		void clear() {
			for (auto n : tbl)
				pool.destruct(n);
			tbl.clear();
		}

		Node *get(D &d, State &s) {
			Node *n = pool.construct();
			d.pack(n->state, s);

			unsigned long hash = n->state.hash(&d);
			Node *found = tbl.find(n->state, hash);
			if (found) {
				pool.destruct(n);
				return found;
			}

			n->goal = d.isgoal(s);
			n->h = n->horig = 0.0;
			n->d = d.d(s);
			tbl.add(n, hash);
			return n;
		}

	private:
		ClosedList<Node, Node, D> tbl;
		Pool<Node> pool;
	};

	class LssNode {
	public:

		LssNode() : openind(-1), updated(false), closed(false) {
		}
	
		class Nodes {
		public:
			static ClosedEntry<LssNode, D> &closedentry(LssNode *n) {
				return n->nodesent;
			}
	
			static PackedState &key(LssNode *n) {
				return n->node->state;
			}
		};

        class FSort {
        public:
            static void setind(LssNode *n, int i) {
                n->openind = i;
            }

            static bool pred(LssNode *a, LssNode *b) {
                double a_f = a->g + a->node->h;
                double b_f = b->g + b->node->h;

                if (geom2d::doubleeq(a_f, b_f))
                    return a->g > b->g;
                return a_f < b_f;
            }
        };

		class GSort {
		public:
			static void setind(LssNode *n, int i) {
				n->openind = i;
			}
		
			static bool pred(LssNode *a, LssNode *b) {
			    return a->g < b->g;
			}
		};

		class HSort {
		public:
			static void setind(LssNode *n, int i) {
				n->openind = i;
			}
		
			static bool pred(LssNode *a, LssNode *b) {
				return a->node->h < b->node->h;
			}
		};

		Node *node;
		LssNode *parent;
		double g, f;
		Oper op;
		long openind;
		bool updated;
		bool closed;

	private:
		ClosedEntry<LssNode, D> nodesent;
	};

public:

	GStar(int argc, const char *argv[]) :
		SearchAlgorithm<D>(argc, argv),
		lssClosed(4051),
		nodes(30000001) {

		lsslim = LookaheadLimit::fromArgs(argc, argv);
	}

	~GStar() {
	}

	void reset() {
		SearchAlgorithm<D>::reset();
		nodes.clear();
		lssOpen.clear();
		for (auto n : lssClosed)
			lssPool.destruct(n);
		lssClosed.clear();
	}

	void search(D &d, State &s0) {
		this->start();

		Node *cur = nodes.get(d, s0);

		lsslim->start(0);

		newIt = true;
		while (!cur->goal && !this->limit()) {

			LssNode *goal = expandLss(d, cur);
			newIt = false;
			if (this->limit())
				break;
			if (!goal)
				hCostLearning(d);

			auto m = move(cur, goal);
			cur = m.first;
			//fprintf( stderr, "%lu\n", lsslim->lookahead(m.second) );
			lsslim->start(m.second);
			times.push_back(walltime() - this->res.wallstart);
			cputimes.push_back(cputime() - this->res.cpustart);
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
			double cpumin = cputimes.front();
			double cpumax = cputimes.front();
			for (unsigned int i = 1; i < times.size(); i++) {
				double dt = times[i] - times[i-1];
				min = std::min(min, dt);
				max = std::max(max, dt);

				dt = cputimes[i] - cputimes[i-1];
				cpumin = std::min(cpumin, dt);
				cpumax = std::max(cpumax, dt);
			}
			dfpair(out, "first emit wall time", "%f", times.front());
			dfpair(out, "min step wall time", "%f", min);
			dfpair(out, "max step wall time", "%f", max);
			dfpair(out, "mean step wall time", "%f", (times.back()-times.front())/times.size());
			dfpair(out, "first emit cpu time", "%f", cputimes.front());
			dfpair(out, "min step cpu time", "%f", cpumin);
			dfpair(out, "max step cpu time", "%f", cpumax);
			dfpair(out, "mean step cpu time", "%f", (cputimes.back()-cputimes.front())/cputimes.size());
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
		lsslim->output(out);
	}

private:

	// ExpandLss returns the cheapest  goal node if one was generated
	// and NULL otherwise.
	LssNode *expandLss(D &d, Node *rootNode) {

	    if( newIt )
	    {
	        lssOpen.clear();
	        for (auto n : lssClosed)
	            lssPool.destruct(n);
	        lssClosed.clear();
	        nclosed = 0;

	        LssNode *a = lssPool.construct();
	        a->node = rootNode;
	        a->parent = NULL;
	        a->op = D::Nop;
	        a->g = 0;
	        a->f = 0;
	        lssOpen.push(a);
	        lssClosed.add(a);
	    }

		LssNode *goal = NULL;

		while (!lssOpen.empty() && !lsslim->stop() && !this->limit()) {
			LssNode *s = *lssOpen.pop();
			nclosed += !s->closed;
			s->closed = true;

			for (auto e : expand(d, s->node)) {
				Node *k = e.node;
				if (s->parent && k == s->parent->node)
					continue;

				LssNode *kid = lssClosed.find(k->state);

				if (!kid) {
					kid = lssPool.construct();
					kid->node = k;
					kid->g = geom2d::Infinity;
					kid->openind = -1;
					lssClosed.add(kid);
				}
				if (kid->g > s->g + e.outcost) {
					if (kid->parent)	// !NULL if a dup
						this->res.dups++;
					kid->parent = s;
					kid->g = s->g + e.outcost;

					kid->op = e.op;
					lssOpen.pushupdate(kid, kid->openind);
				}
				if (k->goal && (!goal || kid->g < goal->g))
					goal = kid;
			}

			if (s->node->goal) {
				lssOpen.push(s);
				goal = s;
				break;
			}
		}

		return goal;
	}

	void hCostLearning(D &d) {
		BinHeap<typename LssNode::HSort, LssNode*> open;
		BinHeap<typename LssNode::GSort, LssNode*> open2;

		for( LssNode *s : lssClosed )
		    s->updated = false;

		minF = NULL;
		for( LssNode *s : lssOpen.data() )
		{
		    State state;
		    state = d.unpack( state, s->node->state );

		    double newH = d.h( state );
		    s->node->h = std::max( s->node->h, newH );
		    s->updated = true;
		    if( !minF || LssNode::FSort::pred( s, minF ) )
		        minF = s;
		}

		open2.append(lssOpen.data());
		open.append(lssOpen.data());

		std::vector<Node*> updated;

		unsigned int temp = nclosed;

		while (nclosed > 0 && !open.empty()) {
			LssNode *s = *open.pop();

			nclosed -= s->closed;

			for (auto e : s->node->preds) {
				Node *sprime = e.node;

				LssNode *sp = lssClosed.find(sprime->state);
				if (!sp || !sp->closed)
					continue;

				if (!sp->updated || sprime->h > e.incost + s->node->h) {
					sprime->h = e.incost + s->node->h;
					sprime->d = s->node->d + 1;
					if (!sp->updated)
						updated.push_back(sprime);
					sp->updated = true;
					open.pushupdate(sp, sp->openind);
				}
			}
		}

		for (auto n : updated)
			n->h = std::max(n->h, n->horig);

		lssOpen.clear( );
		lssOpen.append( open2.data( ) );

		if( !newIt )
		    nclosed = temp;
	}

	std::pair<Node*, double> move(Node *cur, LssNode *goal) {
		//LssNode *best = goal;
		//if (!best) {
		    Oper qop = D::Nop;
		    LssNode *q = NULL;
		    double cost = 0.0;

		    for( outedge e : cur->succs )
		    {
		        LssNode *s = lssClosed.find( e.node->state );

		        if( !q || s->node->h + e.outcost < cost )
		        {
		            q = s;
		            cost = s->node->h + e.outcost;
		            qop = e.op;
		        }
		    }

		    if( q->openind >= 0 || nclosed >= 10000 )
		        newIt = true;

		    //fprintf( stderr, "HIT\n" );

            this->res.ops.push_back(qop);
            lengths.push_back(1);
            return std::make_pair(q->node, q->g);
		/*}
		assert (best);
		assert (best->node != cur);
		std::vector<Oper> ops;

		for (LssNode *p = best; p->node != cur; p = p->parent) {
			assert (p->parent != best);	// no cycles
			ops.push_back(p->op);
		}

		assert (ops.size() >= 1);

		this->res.ops.insert(this->res.ops.end(), ops.rbegin(), ops.rend());
		lengths.push_back(ops.size());
		return std::make_pair(best->node, best->g);*/
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
			k->h = std::max(k->h, n->h - e.cost);
			n->succs.emplace_back(k, ops[i], e.revcost, e.cost);
		}
		n->expd = true;

		return n->succs;
	}

	BinHeap<typename LssNode::GSort, LssNode*> lssOpen;
	ClosedList<typename LssNode::Nodes, LssNode, D> lssClosed;
	Pool<LssNode> lssPool;
	unsigned int nclosed;

	LssNode *minF;
	bool newIt;

	LookaheadLimit *lsslim;
	Nodes nodes;

	std::vector<double> times;
	std::vector<double> cputimes;
	std::vector<unsigned int> lengths;
};
