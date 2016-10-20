// Best-first Utility-Guided Search Yes! From:
// "Best-first Utility-Guided Search," Wheeler Ruml and Minh B. Do,
// Proceedings of the Twentieth International Joint Conference on
// Artificial Intelligence (IJCAI-07), 2007

#include "search.hpp"
#include "../structs/binheap.hpp"
#include "../utils/pool.hpp"
#include <cstring>

template <class D> struct Bugsy : public SearchAlgorithm<D> {
	typedef typename D::State State;
	typedef typename D::PackedState PackedState;
	typedef typename D::Cost Cost;
	typedef typename D::Oper Oper;

	struct Node {
		ClosedEntry<Node, D> closedent;
		int openind;
		Node *parent;
		PackedState state;
		Cost f, g, h, d;
		Oper op, pop;
		double u, t;

		// The expansion count when this node was
		// generated.
		unsigned long expct;

		// path-based mean expansion delay.
		double avgdelay;
		unsigned long depth;

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
			if (a->u != b->u)
				return a->u > b->u;
			if (a->t != b->t)
				return a->t < b->t;
			if (a->f != b->f)
				return a->f < b->f;
			return a->g > b->g;
		}
	};

	enum {
		// Resort1 is the number of expansions to perform
		// before the 1st open list resort.
		Resort1 = 128,
	};

	Bugsy(int argc, const char *argv[]) :
			SearchAlgorithm<D>(argc, argv),
			resort(true),
			usehlms(false),
			usehhat(false),
			usedhat(false),
			navg(0),
			herror(0),
			derror(0),
			herrlast(0),
			derrlast(0),
			useexpdelay(false),
			avgdelay(0),
			delaylast(0),
			dropdups(false),
			lasttime(0.0),
			timeper(0.0),
			lastnodes(0),
			nextresort(Resort1),
			nresort(0),
			closed(30000001) {
		wf = wt = -1;
		for (int i = 0; i < argc; i++) {
			if (i < argc - 1 && strcmp(argv[i], "-wf") == 0)
				wf = strtod(argv[++i], NULL);
			else if (i < argc - 1 && strcmp(argv[i], "-wt") == 0)
				wt = strtod(argv[++i], NULL);
			else if (i < argc - 1 && strcmp(argv[i], "-noresort") == 0)
				resort = false;
			else if (i < argc - 1 && strcmp(argv[i], "-expdelay") == 0)
				useexpdelay = true;
			else if (i < argc - 1 && strcmp(argv[i], "-hhat") == 0)
				usehhat = true;
			else if (i < argc - 1 && strcmp(argv[i], "-dhat") == 0)
				usedhat = true;
			else if (i < argc - 1 && strcmp(argv[i], "-dropdups") == 0)
				dropdups = true;
			else if (i < argc - 1 && strcmp(argv[i], "-interph") == 0)
				interph = usehhat = true;
			else if (i < argc - 1 && strcmp(argv[i], "-hlms") == 0) {
				usehlms = true;
				initlms(sizeof(hcoeffs)/sizeof(hcoeffs[0]), argv[++i], hcoeffs);
			} else if (i < argc - 1 && strcmp(argv[i], "-dlms") == 0) {
				usedlms = true;
				initlms(sizeof(dcoeffs)/sizeof(dcoeffs[0]), argv[++i], dcoeffs);
			}
		}

		if (wf < 0)
			fatal("Must specify non-negative f-weight using -wf");
		if (wt < 0)
			fatal("Must specify non-negative t-weight using -wt");

		nodes = new Pool<Node>();
	}

	~Bugsy() {
		delete nodes;
	}

	void search(D &d, typename D::State &s0) {
		this->start();
		lasttime = walltime();
		closed.init(d);
		Node *n0 = init(d, s0);
		closed.add(n0);
		open.push(n0);

		while (!open.empty() && !SearchAlgorithm<D>::limit()) {
			Node* n = *open.pop();
			State buf, &state = d.unpack(buf, n->state);
			if (d.isgoal(state)) {
				solpath<D, Node>(d, n, this->res);
				break;
			}
			expand(d, n, state);
			updateopen();
		}

		this->finish();
	}

	virtual void reset() {
		SearchAlgorithm<D>::reset();
		open.clear();
		closed.clear();
		delete nodes;
		timeper = 0.0;
		lastnodes = 0;
		nresort = 0;
		nextresort = Resort1;
		navg = 0;
		herror = herrlast = derror = derrlast = 0;
		avgdelay = delaylast = 0;
		nodes = new Pool<Node>();
	}

	virtual void output(FILE *out) {
		SearchAlgorithm<D>::output(out);
		closed.prstats(stdout, "closed ");
		dfpair(stdout, "open list type", "%s", "binary heap");
		dfpair(stdout, "node size", "%u", sizeof(Node));
		dfpair(stdout, "wf", "%g", wf);
		dfpair(stdout, "wt", "%g", wt);
		dfpair(stdout, "final time per expand", "%g", timeper);
		dfpair(stdout, "number of resorts", "%lu", nresort);
		if (usehhat)
			dfpair(stdout, "mean single-step h error", "%g", herror);
		if (usedhat)
			dfpair(stdout, "mean single-step d error", "%g", derror);
		if (useexpdelay)
			dfpair(stdout, "mean expansion delay", "%g", avgdelay);
	}

private:

	// Kidinfo holds information about a node used for
	// correcting the heuristic estimates.
	struct Kidinfo {
		Kidinfo() : f(-1), h(-1), d(-1) { }

		Kidinfo(Cost gval, Cost hval, Cost dval) : f(gval + hval), h(hval), d(dval) { }

		Cost f, h, d;
	};

	// expand expands the node, adding its children to the
	// open and closed lists as appropriate.
	void expand(D &d, Node *n, State &state) {
		this->res.expd++;

		if (useexpdelay) {
			unsigned long delay = this->res.expd - n->expct;
			avgdelay = avgdelay + (delay - avgdelay)/this->res.expd;
		}

		Kidinfo bestinfo;
		typename D::Operators ops(d, state);
		for (unsigned int i = 0; i < ops.size(); i++) {
			if (ops[i] == n->pop)
				continue;

			this->res.gend++;

			Kidinfo kinfo = considerkid(d, n, state, ops[i]);
			if (bestinfo.f < Cost(0) || kinfo.f < bestinfo.f)
				bestinfo = kinfo;
		}

		if (bestinfo.f < Cost(0))
			return;

		navg++;
		if (usehhat) {
			double herr = bestinfo.f - n->f;
			if (herr < 0)	// floating point rounding
				herr = 0;
			herror = herror + (herr - herror)/navg;
		}

		if (usedhat) {
			double derr = bestinfo.d + 1 - n->d;
			if (derr < 0)	// floating point rounding
				derr = 0;
			derror = derror + (derr - derror)/navg;
		}
	}

	// considers adding the child to the open and closed lists.
	Kidinfo considerkid(D &d, Node *parent, State &state, Oper op) {
		Node *kid = nodes->construct();
		typename D::Edge e(d, state, op);
		kid->g = parent->g + e.cost;
		kid->depth = parent->depth + 1;
		kid->avgdelay = parent->avgdelay;

		// single step path-max on d
		kid->d = d.d(e.state);
		if (kid->d < parent->d - Cost(1))
			kid->d = parent->d - Cost(1);

		// single step path-max on h
		kid->h = d.h(e.state);
		if (kid->h < parent->h - e.cost)
			kid->h = parent->h - e.cost;

		kid->f = kid->g + kid->h;

		if (useexpdelay)
			kid->expct = this->res.expd;

		Kidinfo kinfo(kid->g, kid->h, kid->d);

		d.pack(kid->state, e.state);

		unsigned long hash = kid->state.hash(&d);
		Node *dup = static_cast<Node*>(closed.find(kid->state, hash));
		if (dup) {
			this->res.dups++;
			if (!dropdups && kid->g < dup->g) {
				if (dup->openind < 0)
					this->res.reopnd++;
				dup->f = dup->f - dup->g + kid->g;
				dup->g = kid->g;
				dup->parent = parent;
				dup->op = op;
				dup->pop = e.revop;
				computeutil(dup);
				open.pushupdate(dup, dup->openind);
			}
			nodes->destruct(kid);
		} else {
			kid->parent = parent;
			kid->op = op;
			kid->pop = e.revop;
			computeutil(kid);
			closed.add(kid, hash);
			open.push(kid);
		}
		return kinfo;
	}

	Node *init(D &d, State &s0) {
		Node *n0 = nodes->construct();
		d.pack(n0->state, s0);
		n0->g = Cost(0);
		n0->h = n0->f = d.h(s0);
		n0->d = d.d(s0);
		n0->depth = 0;
		n0->avgdelay = 0;
		computeutil(n0);
		n0->op = n0->pop = D::Nop;
		n0->parent = NULL;
		n0->expct = 0;
		return n0;
	}

	// compututil computes the utility value of the given node
	// using corrected estimates of d and h.
	void computeutil(Node *n) {
		if (!resort) {
			derrlast = derror;
			herrlast = herror;
			avgdelay = delaylast;

			// Update time for each expansion if we aren't resorting.
			if (this->res.expd > 0) {
				double t = walltime();
				double dt = t - lasttime;
				lastnodes++;
				timeper += (dt - timeper)/this->res.expd;
				lasttime = t;
			}
		}
		double d = usedlms ? evallms(n, dcoeffs) : n->d;
		if (usedhat)
			d /= (1 - derrlast);

		double h = usehlms ? evallms(n, hcoeffs) : n->h;
		if (usehhat) {
			if (interph) {
				double hhat = h + d*herrlast;
				double w = wf/(wf+wt);
				h = n->h*w + hhat*(1-w);
			} else {
				h += d * herrlast;
			}
		}

		if (useexpdelay && delaylast > 0)
			d *= delaylast;

		double f = h + n->g;
		n->t = timeper * d;
		n->u = -(wf * f + wt * n->t);
	}

	// updateopen updates the utilities of all nodes on open and
	// reinitializes the heap every 2^i expansions.
	void updateopen() {
		if (!resort || this->res.expd < nextresort)
			return;
		nextresort *= 2;
		nresort++;
		herrlast = herror;
		derrlast = derror;
		delaylast = avgdelay;

		double t = walltime();
		timeper = (t-lasttime) / (this->res.expd - lastnodes);
		lasttime = t;
		lastnodes = this->res.expd;

		reinitheap();
	}

	// Reinitheap reinitialize the heap property in O(n)
	// time while also updating the utilities.
	void reinitheap() {
		if (open.size() <= 0)
			return;

		for (long i = open.size()-1; i > (long) open.size()/2; i--) 
			computeutil(open.at(i));

		for (long i = (long) open.size()/2; i >= 0; i--) {
			computeutil(open.at(i));
			open.pushdown(i);
		}
	}

	void initlms(unsigned int n, const char *wts, double coeffs[]) {
		char *cpy = new char[strlen(wts)+1];
		memcpy(cpy, wts, strlen(wts)+1);

		char *tok = strtok(cpy, ",");
		for (unsigned int i = 0; i < n; i++) {
			if (!tok)
				fatal("Expected %u coefficients", n);
			char *end = NULL;
			coeffs[i] = strtod(tok, &end);
			if (end == tok)
				fatal("Coefficient %u (%s) is invalid", i, tok);		
			tok = strtok(NULL, ",");
		}

		delete []cpy;
	}

	double evallms(Node *n, double coeffs[]) const {
		return n->h*coeffs[0] +
			n->g*coeffs[1] +
			n->d*coeffs[2] +
			n->depth*coeffs[3];
	}

	// resort the open list sometimes?
	bool resort;

	bool usehlms;
	// h, g, d, and D coefficients for LMS-based heuristic;
	double hcoeffs[4];

	bool usedlms;
	// h, g, d, and D coefficients for LMS-based heuristic;
	double dcoeffs[4];

	// wf and wt are the cost and time weight respectively.
	double wf, wt;

	// heuristic correction
	bool usehhat, usedhat;
	unsigned long navg;
	double herror, derror;
	// Herrlast and derrlast are the h and d errors at
	// the last re-sort.  These are used to order nodes
	// until the next resort.
	double herrlast, derrlast;

	// interph — when using heuristic correction,
	// interpolate between h and hhat depending
	// on the ratio of wf to wt.
	bool interph;

	// expansion delay
	bool useexpdelay;
	double avgdelay;
	// Delaylast is the expansion delay at the last resort,
	// it is used to order nodes until the next resort.
	double delaylast;

	bool dropdups;

	// Time-per-expansion estimation.
	// Lasttime is the last time the estimate was made.
	// Lastnodes is the expansion count at the last estimate.
	// Timeper is the time-per-expansion at the last estimate.
	double lasttime, timeper;
	unsigned long lastnodes;

	// for resorting the open list
	unsigned long nextresort;
	unsigned int nresort;

	BinHeap<Node, Node*> open;
 	ClosedList<Node, Node, D> closed;
	Pool<Node> *nodes;
};
