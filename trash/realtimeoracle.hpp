#pragma once
#include "search.hpp"
#include "../utils/geom2d.hpp"
#include "../utils/pool.hpp"
#include "../utils/utils.hpp"
#include "../rdb/rdb.hpp"
#include <vector>

template <class D>
class Realtimeoracle : public SearchAlgorithm<D> {
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

		void updateH(double newH) {
			h = newH;
			hHistory.push_back(newH);
		}

		void swapH(double newH) {
			assert(hHistory.size() > 0);
			revertH();
			updateH(newH);
		}

		void revertH() {
			assert(hHistory.size() > 0);

			h = hHistory.back();
			hHistory.pop_back();
		}

		std::vector<double> hHistory;
		PackedState state;
		double h, horig;
		bool dead;
		bool expd;	// was this expanded before?
		bool goal;
		std::vector<outedge> succs;
		std::vector<inedge> preds;

	private:
		friend class Nodes;
		friend class Pool<Node>;

		Node() : dead(false), expd(false) {
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
			n->updateH(d.h(s));
			n->horig = n->h;
//			n->h = n->horig = d.h(s);
			tbl.add(n, hash);
			return n;
		}

		void output(FILE *out) {
			tbl.prstats(out, "nodes ");
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
				if (geom2d::doubleeq(a->f, b->f))
					return a->g > b->g;
				return a->f < b->f;
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

	class GATCost {
	public:
		GATCost() : gat(0), solutionCost(0), searchTime(0), depth(0) {}

		void finalize() {
			std::reverse(lssSizes.begin(), lssSizes.end());
			std::reverse(lssTimes.begin(), lssTimes.end());
			std::reverse(solutionPath.begin(), solutionPath.end());
		}

		double gat;
		double solutionCost;
		double searchTime;
		unsigned int depth;
		std::vector<unsigned int> lssSizes;
		std::vector<double> lssTimes;
		std::vector<Oper> solutionPath;
	};

public:

	Realtimeoracle(int argc, const char *argv[]) :
		SearchAlgorithm<D>(argc, argv),
		nodes(30000001), frameRate(0), maxLSS(0), deltaLSS(0) {

		for (int i = 0; i < argc; i++)
			if (strcmp(argv[i], "-framerate") == 0) 
				frameRate = strtod(argv[++i], NULL);
			else if (strcmp(argv[i], "-maxlss") == 0) 
				maxLSS = strtol(argv[++i], NULL, 10);
			else if (strcmp(argv[i], "-deltalss") == 0)
				deltaLSS = strtol(argv[++i], NULL, 10);

		if(frameRate == 0)
			fatal("please supply a \"frame rate\" with -framerate");
		else if (maxLSS == 0)
			fatal("please supply a max lss sizewith -maxlss");
		else if(deltaLSS == 0)
			fatal("please supply a delta lss size with -deltalss");
	}

	~Realtimeoracle() {
	}

	void reset() {
		SearchAlgorithm<D>::reset();
		nodes.clear();
	}

	void search(D &d, State &s0) {
		this->start();

		Node *cur = nodes.get(d, s0);

		GATCost zeroPath;
		GATCost infIncumbent;
		infIncumbent.gat = std::numeric_limits<double>::infinity();
		infIncumbent.solutionCost = std::numeric_limits<double>::infinity();
		infIncumbent.searchTime = std::numeric_limits<double>::infinity();

		GATCost cost = searchHelper(d, cur, zeroPath, infIncumbent, frameRate);

		cost.finalize();
		this->res.ops = cost.solutionPath;

		fprintf(stderr, "Solution Cost: %f\n" , cost.solutionCost);
		fprintf(stderr, "Search Time: %f\n", cost.searchTime);
		fprintf(stderr, "GAT: %f\n", cost.gat);

		fprintf(stderr, "LSS Sizes:\n");
		for(auto lss : cost.lssSizes)
			fprintf(stderr, "\t%u\n", lss);
		fprintf(stderr, "LSS Times\n");
		for(auto lss : cost.lssTimes)
			fprintf(stderr, "\t%f\n", lss);

		this->finish();

		if (std::isinf(cost.solutionCost)) {
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
		nodes.output(stdout);
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

	GATCost searchHelper(D &d, Node *root, const GATCost &path, const GATCost &incumbent, double planningTime) {
		GATCost cost;
		cost.gat = std::numeric_limits<double>::infinity();
		cost.solutionCost = std::numeric_limits<double>::infinity();
		cost.searchTime = std::numeric_limits<double>::infinity();

		BinHeap<typename LssNode::FSort, LssNode*> *lssOpen;
		ClosedList<typename LssNode::Nodes, LssNode, D> *lssNodes;
		Pool<LssNode> *lssPool;

		for(int effort = maxLSS; effort >= 1; effort -= deltaLSS) {
			lssOpen = new BinHeap<typename LssNode::FSort, LssNode*>();
			lssNodes = new ClosedList<typename LssNode::Nodes, LssNode, D>(4051);
			lssPool = new Pool<LssNode>();

			double startTime = cputime();

			LssNode* goal = doSearchIteration(d, root, effort, *lssOpen, *lssNodes, *lssPool);

			//If we didn't find the goal do some h learning
			if (!goal) {
				hCostLearning(d, *lssOpen, *lssNodes);
			}
			else fprintf(stderr, "\tfound goal\n");

			//retrieve the best node on the fringe and its g-cost from the lss
			LssNode* best = getBestMove(root, goal, *lssOpen);

			double gvalue = best->g;
			Node* bestNode = best->node;
			double hvalue = bestNode->h;

			std::vector<Oper> ops;
			for (LssNode *p = best; p->node != root; p = p->parent) {
				assert (p->parent != best);	// no cycles
				ops.push_back(p->op);
			}
			assert (ops.size() >= 1);

			//how long did this searching take us
			double deltaTime = cputime() - startTime;
			
			// The amount of execution time we'll get for our current leg
			double execTime = frameRate * best->g;
			
			// how much time is left over from the previous planning iteration
			double remainingTime = planningTime - deltaTime;

			double gat = execTime;

			if(remainingTime < 0) {
				gat += -remainingTime;
				remainingTime = frameRate - fmod(-remainingTime, frameRate);
			}

			//clean up the damn mess
			delete lssOpen;
			delete lssNodes;
			delete lssPool;

			GATCost deeperCost;
			if(!goal) {
				//keeping track of the cost so far
				GATCost morePath = path;
				morePath.solutionCost += gvalue;
				morePath.searchTime += deltaTime;
				morePath.gat += gat;

				//don't bother continuing if we're worse than the incumbent
				if(incumbent.gat < (morePath.gat+ hvalue)) {
					continue;
				}

				//dig deeper
				deeperCost = searchHelper(d, bestNode, morePath, cost, remainingTime);

				//however long it took, add this current iteration costs
				deeperCost.solutionCost += gvalue;
				deeperCost.searchTime += deltaTime;
				deeperCost.gat += gat;

			}
			else {
				//we found the goal so just return our current stats
				deeperCost.solutionCost = gvalue;
				deeperCost.searchTime = deltaTime;
				deeperCost.gat = gat;
			}

			if(deeperCost.gat < cost.gat) {
				cost.solutionCost = deeperCost.solutionCost;
				cost.searchTime = deeperCost.searchTime;
				cost.gat = deeperCost.gat;
				cost.depth = deeperCost.depth + 1;
				cost.lssSizes = deeperCost.lssSizes;
				cost.lssSizes.push_back(effort);
				cost.lssTimes = deeperCost.lssTimes;
				cost.lssTimes.push_back(deltaTime);

				cost.solutionPath = deeperCost.solutionPath;
				cost.solutionPath.insert(cost.solutionPath.end(), ops.begin(), ops.end());
			}

		}

		return cost;
	}

	LssNode *doSearchIteration(D &d, Node *rootNode, unsigned int expansions,
									BinHeap<typename LssNode::FSort, LssNode*> &lssOpen,
									ClosedList<typename LssNode::Nodes, LssNode, D> &lssNodes,
									Pool<LssNode> &lssPool) {
		LssNode *a = lssPool.construct();
		a->node = rootNode;
		a->parent = NULL;
		a->op = D::Nop;
		a->g = 0;
		a->f = rootNode->h;
		lssOpen.push(a);
		lssNodes.add(a);

		nclosed = 0;
		unsigned int nexpanded = 0;

		LssNode *goal = NULL;

		while (!lssOpen.empty() && !this->limit() && nexpanded < expansions) {

			LssNode *s = *lssOpen.pop();

			nclosed += !s->closed;
			s->closed = true;

			nexpanded++;

			for (auto e : expand(d, s->node)) {
				Node *k = e.node;
				if (s->parent && k == s->parent->node)
					continue;

				LssNode *kid = lssNodes.find(k->state);

				if (!kid) {
					kid = lssPool.construct();
					kid->node = k;
					kid->parent = NULL;
					kid->g = geom2d::Infinity;
					lssNodes.add(kid);
				}
				if (kid->g > s->g + e.outcost) {
					if (kid->parent)	// !NULL if a dup
						this->res.dups++;
					kid->parent = s;
					kid->g = s->g + e.outcost;
					kid->f = kid->g + kid->node->h;
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

	void hCostLearning(D &d, BinHeap<typename LssNode::FSort, LssNode*> &lssOpen,
					ClosedList<typename LssNode::Nodes, LssNode, D> &lssNodes) {
		BinHeap<typename LssNode::HSort, LssNode*> open;

		open.append(lssOpen.data());

		std::vector<Node*> updated;

		while (nclosed > 0 && !open.empty()) {
			LssNode *s = *open.pop();

			nclosed -= s->closed;

			for (auto e : s->node->preds) {
				Node *sprime = e.node;

				LssNode *sp = lssNodes.find(sprime->state);
				if (!sp || !sp->closed)
					continue;

				if (!sp->updated || sprime->h > e.incost + s->node->h) {
					sprime->updateH(e.incost + s->node->h);
					if (!sp->updated)
						updated.push_back(sprime);
					sp->updated = true;
					open.pushupdate(sp, sp->openind);
				}
			}
		}

		for (auto n : updated)
			n->swapH(std::max(n->h, n->horig));
//If we want the heuristic to be strictly increasing, this line might need to also update n->horig

	}

	LssNode* getBestMove(Node *cur, LssNode *goal,
									BinHeap<typename LssNode::FSort, LssNode*> &lssOpen) {
		LssNode *best = goal;
		if (!best) {
			for (auto n : lssOpen.data()) {
				if (n->node == cur)
					continue;
				if (best == NULL || LssNode::FSort::pred(n, best))
					best = n;
			}
		}
		assert (best);
		assert (best->node != cur);

		return best;
	}

	// Expand returns the successor nodes of a state.
	std::vector<outedge> expand(D &d, Node *n) {
		assert(!n->dead);

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
			k->swapH(std::max(k->h, n->h - e.cost));
			n->succs.emplace_back(k, ops[i], e.revcost, e.cost);
		}
		n->expd = true;

		return n->succs;
	}

	Nodes nodes;

	double frameRate;

	unsigned int nclosed;
	unsigned int maxLSS;
	unsigned int deltaLSS;

	std::vector<double> times;
	std::vector<unsigned int> lengths;
};
