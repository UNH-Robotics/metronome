#pragma once
#include <cassert>
#include <isolation/search/closedlist.hpp>
#include <vector>
#include "../utils/geom2d.hpp"
#include "../utils/pool.hpp"
#include "lsslrtastar2.hpp"
#include "search.hpp"

#define PI 3.14159265359

template <class D>
class MORTS : public SearchAlgorithm<D> {
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
        inedge(Node* n, double c) : node(n), incost(c) {}

        Node* node;
        double incost;
    };

    struct outedge {
        outedge(Node* n, Oper o, Cost rc, Cost c)
                : node(n),
                  op(o),
                  revcost(rc == Cost(-1) ? geom2d::Infinity : rc),
                  outcost(c == Cost(-1) ? geom2d::Infinity : c) {}

        Node* node;
        Oper op;
        double revcost;
        double outcost;
    };

    class Node {
    public:
        static ClosedEntry<Node, D>& closedentry(Node* n) { return n->closedent; }

        static PackedState& key(Node* n) { return n->state; }

        PackedState state;
        double h, horig, d, derr;
        bool expd; // Was this expanded before?
        bool goal;

        // Nop if node has no identity was not expanded, otherwise it's the identity action at this node.
        std::pair<Oper, Cost> ident;

        std::vector<outedge> succs;
        std::vector<inedge> preds;

    private:
        friend class Nodes;
        friend class Pool<Node>;

        Node() : expd(false) {}

        ClosedEntry<Node, D> closedent;
    };

    class Nodes {
    public:
        Nodes(unsigned int sz) : tbl(sz) {}

        ~Nodes() {}

        void clear() {
            for (auto n : tbl)
                pool.destruct(n);
            tbl.clear();
        }

        Node* get(D& d, State& s) {
            Node* n = pool.construct();
            d.pack(n->state, s);

            unsigned long hash = n->state.hash(&d);
            Node* found = tbl.find(n->state, hash);
            if (found) {
                pool.destruct(n);
                return found;
            }

            n->goal = d.isgoal(s);
            n->h = n->horig = d.h(s);
            n->d = n->derr = d.d(s);
            tbl.add(n, hash);
            return n;
        }

    private:
        ClosedList<Node, Node, D> tbl;
        Pool<Node> pool;
    };

    class LssNode {
    public:
        LssNode() : openind(-1), updated(false), closed(false) {}

        class Nodes {
        public:
            static ClosedEntry<LssNode, D>& closedentry(LssNode* n) { return n->nodesent; }

            static PackedState& key(LssNode* n) { return n->node->state; }
        };

        class FHatSort {
        public:
            static void setind(LssNode* n, int i) { n->openind = i; }

            static bool pred(LssNode* a, LssNode* b) {
                if (geom2d::doubleeq(a->fhat, b->fhat)) {
                    if (geom2d::doubleeq(a->f, b->f))
                        return a->g > b->g;
                    return a->f < b->f;
                }
                return a->fhat < b->fhat;
            }
        };

        class HSort {
        public:
            static void setind(LssNode* n, int i) { n->openind = i; }

            static bool pred(LssNode* a, LssNode* b) { return a->node->h < b->node->h; }
        };

        Node* node;
        LssNode *parent, *tla, *best;
        unsigned long genCount, expCount, depth; // Counters for node generation and expansion.
        unsigned long iteration;
        double g, f, fhat, herr;
        Oper op;
        long openind;
        bool updated, closed, path;

    private:
        ClosedEntry<LssNode, D> nodesent;
    };

public:
    MORTS(int argc, const char* argv[])
            : SearchAlgorithm<D>(argc, argv),
              lssClosed(4051),
              herror(0),
              derror(0),
              nodes(30000001),
              onestep(false),
              delaySum(0),
              delayCount(0),
              nshort(0),
              nident(0),
              identStates(0),
              identCurs(0) {
        argv[1] = "fhatlrtastar";
        lsslim = LookaheadLimit::fromArgs(argc, argv);
        P = std::vector<Oper>();

        for (int i = 0; i < argc; i++) {
            if (strcmp(argv[i], "-onestep") == 0)
                onestep = true;
        }
    }

    ~MORTS() {}

    void reset() {
        SearchAlgorithm<D>::reset();
        nodes.clear();
        lssOpen.clear();
        for (auto n : lssClosed)
            lssPool.destruct(n);
        lssClosed.clear();
        P.clear();
        herror = 0;
        derror = 0;
        nshort = 0;
        nident = 0;
        identStates = 0;
        identCurs = 0;
    }

    void search(D& d, State& s0) {
        this->start();

        Node* cur = nodes.get(d, s0);

        lsslim->start(0);

        bool ident = false;

        curIteration = 0;
        expCount = 0;

        while (!cur->goal && !this->limit()) {
            if (cur->ident.first != D::Nop)
                identCurs++;

            LssNode* goal = expandLss(d, cur, ident);
            if (this->limit())
                break;
            auto m = move(d, cur, goal);
            ident = (cur == m.first);
            cur = m.first;
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

    virtual void output(FILE* out) {
        SearchAlgorithm<D>::output(out);
        dfpair(out, "num steps", "%lu", (unsigned long)times.size());
        assert(lengths.size() == times.size());
        // dfpair(out, "h error last", "%g", herror);
        // dfpair(out, "d error last", "%g", derror);
        if (times.size() != 0) {
            double min = times.front();
            double max = times.front();
            double cpumin = cputimes.front();
            double cpumax = cputimes.front();
            double omit = 0.0;
            for (unsigned int i = 1; i < times.size(); i++) {
                double dt = times[i] - times[i - 1];
                min = std::min(min, dt);
                max = std::max(max, dt);

                dt = cputimes[i] - cputimes[i - 1];
                cpumin = std::min(cpumin, dt);
                cpumax = std::max(cpumax, dt);
            }

            for (unsigned int i = 1; i < omitTimes.size(); i += 2) {
                omit += omitTimes[i] - omitTimes[i - 1];
            }

            dfpair(out, "time to omit", "%f", omit);
            dfpair(out, "first emit wall time", "%f", times.front());
            dfpair(out, "min step wall time", "%f", min);
            dfpair(out, "max step wall time", "%f", max);
            dfpair(out, "mean step wall time", "%f", (times.back() - times.front()) / times.size());
            dfpair(out, "first emit cpu time", "%f", cputimes.front());
            dfpair(out, "min step cpu time", "%f", cpumin);
            dfpair(out, "max step cpu time", "%f", cpumax);
            dfpair(out, "mean step cpu time", "%f", (cputimes.back() - cputimes.front()) / cputimes.size());
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
            dfpair(out, "mean step length", "%g", sum / (double)lengths.size());
        }
        lsslim->output(out);
        dfpair(out, "number of short trajectories taken", "%u", nshort);
        dfpair(out, "number of identity actions taken", "%u", nident);
        dfpair(out, "number of expanded states with identity actions", "%lu", identStates);
        dfpair(out, "number of current states with identity actions", "%lu", identCurs);
    }

private:
    // ExpandLss returns the cheapest  goal node if one was generated
    // and NULL otherwise.
    LssNode* expandLss(D& d, Node* rootNode, bool ident) {
        static double herrnext = 0;
        static double derrnext = 0;

        ++curIteration;

        delaySum = 0.0;
        delayCount = 0;

        if (!ident) {
            lssOpen.clear();
            for (auto n : lssClosed)
                lssPool.destruct(n);
            lssClosed.clear();
            nclosed = 0;

            LssNode* a = lssPool.construct();
            a->node = rootNode;
            a->parent = NULL;
            a->tla = NULL;
            a->op = D::Nop;
            a->g = 0;
            a->depth = 0;
            a->iteration = curIteration;
            a->genCount = 0;
            a->f = rootNode->h;
            lssOpen.push(a);
            lssClosed.add(a);

            herrnext = 0;
            derrnext = 0;
            expCount = 0;
        }

        LssNode* goal = NULL;

        while (!lssOpen.empty() && !lsslim->stop() && !this->limit()) {
            LssNode* s = *lssOpen.pop();

            nclosed += !s->closed;
            s->closed = true;
            s->expCount = ++expCount;

            delaySum += s->expCount - s->genCount;
            ++delayCount;

            LssNode* bestkid = NULL;
            for (auto e : expand(d, s->node)) {
                Node* k = e.node;
                if (s->parent && k == s->parent->node)
                    continue;

                LssNode* kid = lssClosed.find(k->state);

                if (!kid) {
                    kid = lssPool.construct();
                    kid->node = k;
                    kid->g = geom2d::Infinity;
                    kid->openind = -1;
                    kid->genCount = expCount;
                    lssClosed.add(kid);
                }
                if (kid->g > s->g + e.outcost) {
                    kid->iteration = curIteration;

                    kid->depth = s->depth + 1;

                    if (kid->parent) // !NULL if a dup
                        this->res.dups++;

                    kid->parent = s;
                    kid->g = s->g + e.outcost;
                    kid->f = kid->g + kid->node->h;

                    kid->tla = kid->parent->tla;
                    if (!kid->tla)
                        kid->tla = kid;

                    assert(derror != 1);
                    double d = kid->node->derr / (1 - derror);
                    double h = kid->node->h + herror * d;
                    // kid->herr = std::max( 0.0, kid->f - s->f );
                    // kid->herr = s->herr + ( kid->herr - s->herr ) / kid->depth;
                    kid->fhat = kid->g + h; // kid->node->h + kid->herr * kid->node->d;

                    kid->op = e.op;
                    lssOpen.pushupdate(kid, kid->openind);
                }
                if (k->goal && (!goal || kid->g < goal->g))
                    goal = kid;

                if (!bestkid || kid->f < bestkid->f)
                    bestkid = kid;
            }

            if (bestkid) {
                double herr = bestkid->f - s->f;
                if (herr < 0)
                    herr = 0;
                herrnext = herrnext + (herr - herrnext) / (expCount + 1);

                double derr = bestkid->node->d + 1 - s->node->d;
                if (derr < 0)
                    derr = 0;
                if (derr >= 1)
                    derr = 1 - geom2d::Threshold;
                derrnext = derrnext + (derr - derrnext) / (expCount + 1);
            }

            if (s->node->goal) {
                lssOpen.push(s);
                goal = s;
                break;
            }
        }

        herror = herrnext;
        if (derrnext == 1)
            derrnext = 0;
        derror = derrnext;

        return goal;
    }

    void hCostLearning(D& d) {
        BinHeap<typename LssNode::HSort, LssNode*> open;

        open.append(lssOpen.data());

        std::vector<Node*> updated;

        while (nclosed > 0 && !open.empty()) {
            LssNode* s = *open.pop();

            nclosed -= s->closed;

            for (auto e : s->node->preds) {
                Node* sprime = e.node;

                LssNode* sp = lssClosed.find(sprime->state);
                if (!sp || !sp->closed)
                    continue;

                if (!sp->updated || sprime->h > e.incost + s->node->h) {
                    sprime->h = e.incost + s->node->h;
                    sprime->derr = s->node->derr;
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

        for (LssNode* s : lssClosed)
            s->path = false;

        LssNode* a = *lssOpen.pop();
        LssNode* s = a->parent;
        int count = 0;

        while (s) {
            ++count;
            s->path = true;
            s->best = NULL;
            s = s->parent;
        }

        while (count > 0 && !lssOpen.empty()) {
            s = *lssOpen.pop();

            if (s->tla != a->tla)
                continue;

            LssNode* sPrime = s->parent;

            while (!sPrime->path)
                sPrime = sPrime->parent;

            if (!sPrime->best) {
                sPrime->best = s;
                --count;
            }
        }
    }

    std::pair<Node*, double> move(D& d, Node* cur, LssNode* goal) {
        LssNode* best = goal;
        LssNode* secondBest = NULL;
        if (!best) {
            for (auto n : lssOpen.data()) {
                if (n->node == cur)
                    continue;
                if (best == NULL || LssNode::FHatSort::pred(n, best))
                    best = n;
            }
            for (auto n : lssOpen.data()) {
                if (n->node == cur || n->tla == best->tla)
                    continue;
                if (secondBest == NULL || LssNode::FHatSort::pred(n, secondBest))
                    secondBest = n;
            }
        }
        assert(best);
        assert(best->node != cur);
        std::vector<Oper> ops;

        LssNode* to = best;
        assert(to->parent);
        omitTimes.push_back(cputime());

        if (!goal) {
            if (cur->ident.first != D::Nop &&
                    computeBenefit(d, cur, best, secondBest, cur->ident.second) > cur->ident.second) {
                this->res.ops.push_back(cur->ident.first);
                lengths.push_back(1);
                nshort++;
                nident++;
                return std::make_pair(cur, cur->ident.second);

                // If we are using prefix commitment, available quanta are all >= time to state
                // We only have > options if we have an identity action at/before state
                // If we are using some other form of commitment we have quanta available < time to state
                // The best state is the one with highest Utility(state, quanta)

                // WHAT IS UTILITY?
                // Benefit - Cost
                // Benefit = average (currently) possible reduction in search time (currently reduction in f) *
                // probability of reduction
                // Cost = guaranteed loss (time spent not moving) + possible loss (bypassed possible reduction)

                // If the best state is utilizing an inherited identity action, commit to the state from which it
                // inherits,
                // but search at the state itself, and hold on to the remaining prefix
            }

            hCostLearning(d);

            LssNode *alpha = best, *beta;

            std::vector<LssNode*> path;

            for (LssNode* p = best; p->node != cur && p->iteration == curIteration; p = p->parent)
                path.push_back(p);

            std::reverse(path.begin(), path.end());

            // LssNode *prefix = NULL;
            // double maxBenefit = 0.0, benefit;

            // THIS IS STILL REAL TIME,
            // SINCE IT'S ONLY DONE FOR THE LAST ITERATION

            for (LssNode* p : path) {
                to = p;

                if (to == best)
                    break;

                beta = to->best;

                if (!beta)
                    continue;

                double area = (std::max(alpha->depth, beta->depth) - to->depth) * (delaySum / delayCount);

                if (computeBenefit(d, to->node, alpha, beta, to->g) > area / lsslim->lookahead(1))
                    break;
            }

            if (to != best)
                ++nshort;
        }

        omitTimes.push_back(cputime());

        for (LssNode* p = to; p->node != cur; p = p->parent) {
            assert(p->parent != to); // no cycles
            ops.push_back(p->op);
        }

        assert(ops.size() >= 1);

        this->res.ops.insert(this->res.ops.end(), ops.rbegin(), ops.rend());
        lengths.push_back(ops.size());
        return std::make_pair(to->node, to->g);
    }

    // Returns the utility of searching at a given state
    double computeBenefit(D& d, Node* n, LssNode* alpha, LssNode* beta, double g) {
        if (beta == NULL || alpha == beta)
            return 0.0;

        double x = lsslim->lookahead(g);

        int i = 1;
        double delay_alpha = 0;
        for (LssNode* p = alpha->parent; p->genCount > 0 && p->node != n && delay_alpha < x; p = p->parent) {
            ++i;
            delay_alpha += p->expCount - p->genCount;
        }
        delay_alpha += 1 + expCount - alpha->genCount;
        delay_alpha /= i;

        i = 1;
        double delay_beta = 0;
        for (LssNode* p = beta->parent; p->genCount > 0 && p->node != n && delay_beta < x; p = p->parent) {
            ++i;
            delay_beta += p->expCount - p->genCount;
        }
        delay_beta += 1 + expCount +
                (beta->fhat - alpha->fhat) /
                        (std::abs(alpha->fhat - lssClosed.find(n->state)->fhat) / alpha->genCount) -
                beta->genCount;

        delay_beta /= i;

        LssNode* cur = lssClosed.find(n->state);

        double bv_alpha = pow((cur->f - alpha->f) / (alpha->depth - cur->depth) * alpha->node->d, 2),
               bv_beta = pow((cur->f - beta->f) / (beta->depth - cur->depth) * beta->node->d, 2);

        double mu_alpha = alpha->f + sqrt(bv_alpha), mu_beta = beta->f + sqrt(bv_beta);

        double x_alpha = std::min(alpha->node->d, x / delay_alpha), x_beta = std::min(beta->node->d, x / delay_beta);

        // TODO if not only considering identity actions, fix this
        double v_alpha = bv_alpha * (x_alpha / alpha->node->d), v_beta = bv_beta * (x_beta / beta->node->d);

        double start = (mu_alpha - 2.0 * sqrt(v_alpha) < mu_beta - 2.0 * sqrt(v_beta)) ?
                mu_alpha - 2.0 * sqrt(v_alpha) :
                mu_beta - 2.0 * sqrt(v_beta);
        double end = (mu_alpha + 2.0 * sqrt(v_alpha) > mu_beta + 2.0 * sqrt(v_beta)) ? mu_alpha + 2.0 * sqrt(v_alpha) :
                                                                                       mu_beta + 2.0 * sqrt(v_beta);
        double a, b;
        double benefit = 0.0;

        for (int i = 0; i < 100; ++i) {
            a = start + ((end - start) * ((i + 0.5) / 100.0));

            int j;
            double sum = 0.0;
            for (j = 0; j < 100; ++j) {
                b = start + ((end - start) * ((j + 0.5) / 100.0));

                if (a < b)
                    break;

                sum += (a - b) * std::exp(-pow(b - mu_beta, 2) / (2 * v_beta)) / (sqrt(2 * PI * v_beta));
            }

            sum *= (end - start) / 100.0;

            benefit += sum * std::exp(-pow(a - mu_alpha, 2) / (2 * v_alpha)) / (sqrt(2 * PI * v_alpha)) *
                    ((end - start) / 100.0);
        }

        return benefit;
    }

    // Expand returns the successor nodes of a state.
    std::vector<outedge> expand(D& d, Node* n) {
        this->res.expd++;

        if (n->expd)
            return n->succs;

        State buf, &s = d.unpack(buf, n->state);

        n->ident = d.ident(s);
        if (n->ident.first != D::Nop)
            identStates++;

        Operators ops(d, s);
        for (unsigned int i = 0; i < ops.size(); i++) {
            this->res.gend++;

            Edge e(d, s, ops[i]);
            Node* k = nodes.get(d, e.state);
            k->preds.emplace_back(n, e.cost);
            k->h = std::max(k->h, n->h - e.cost);
            n->succs.emplace_back(k, ops[i], e.revcost, e.cost);
        }
        n->expd = true;

        return n->succs;
    }

    BinHeap<typename LssNode::FHatSort, LssNode*> lssOpen;
    ClosedList<typename LssNode::Nodes, LssNode, D> lssClosed;
    Pool<LssNode> lssPool;
    unsigned int nclosed;

    double herror;
    double derror;

    LookaheadLimit* lsslim;
    Nodes nodes;

    bool onestep;

    unsigned long expCount, delaySum, delayCount;

    unsigned long curIteration;

    unsigned int nshort, nident; // Number of short and identity trajectories taken.
    unsigned long identStates; // Number of expanded states with identity actions.
    unsigned long identCurs; // Number of 'current states' with identity actions.

    std::vector<Oper> P;

    std::vector<double> times;
    std::vector<double> cputimes;
    std::vector<double> omitTimes;
    std::vector<unsigned int> lengths;
};
