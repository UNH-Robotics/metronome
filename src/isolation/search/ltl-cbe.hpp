#pragma once
#include "search.hpp"
#include "../utils/geom2d.hpp"
#include "../utils/pool.hpp"
#include "lsslrtastar2.hpp"
#include <vector>

#define PI 3.14159265359

template <class D>
class LTL_CBE : public SearchAlgorithm<D> {
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
        double h, horig, d, derr;
        bool expd;  // Was this expanded before?
        bool goal;

        // Nop if node has no identity was not expanded, otherwise it's the identity action at this node.
        std::pair<Oper,Cost> ident;

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

        class FHatSort {
        public:
            static void setind(LssNode *n, int i) {
                n->openind = i;
            }

            static bool pred(LssNode *a, LssNode *b) {
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
            static void setind(LssNode *n, int i) {
                n->openind = i;
            }

            static bool pred(LssNode *a, LssNode *b) {
                return a->node->h < b->node->h;
            }
        };

        Node *node;
        LssNode *parent;
        unsigned long genCount, expCount; //Counters for node generation and expansion.
        double g, f, fhat;
        Oper op;
        long openind;
        bool updated;
        bool closed;

    private:
        ClosedEntry<LssNode, D> nodesent;
    };

public:

    LTL_CBE(int argc, const char *argv[]) :
        SearchAlgorithm<D>(argc, argv),
        lssClosed(4051),
        herror(0),
        derror(0),
        nodes(30000001),
        nshort(0),
        nident(0),
        identStates(0),
        identCurs(0) {

        argv[1] = "fhatlrtastar";
        lsslim = LookaheadLimit::fromArgs(argc, argv);
        P = std::vector<Oper>();
    }

    ~LTL_CBE() {
    }

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

    void search(D &d, State &s0) {
        this->start();

        Node *cur = nodes.get(d, s0);
        Node *root = cur;
        double leftOver = 0;

        lsslim->start(0);

        while (!cur->goal && !this->limit()) {

            if (cur->ident.first != D::Nop)
                identCurs++;

            LssNode *goal = expandLss(d, root);
            if (this->limit())
                break;
            if (!goal)
                hCostLearning(d);
            auto m = move(d, cur, root, goal, leftOver);
            cur = m.first.first;
            root = m.first.second;
            lsslim->start(m.second.first);
            leftOver = m.second.second;
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
        dfpair(out, "h error last", "%g", herror);
        dfpair(out, "d error last", "%g", derror);
        if (times.size() != 0) {
            double min = times.front();
            double max = times.front();
            double cpumin = cputimes.front();
            double cpumax = cputimes.front();
            double omit = 0.0;
            for (unsigned int i = 1; i < times.size(); i++) {
                double dt = times[i] - times[i-1];
                min = std::min(min, dt);
                max = std::max(max, dt);

                dt = cputimes[i] - cputimes[i-1];
                cpumin = std::min(cpumin, dt);
                cpumax = std::max(cpumax, dt);
            }

            for (unsigned int i = 1; i < omitTimes.size(); i += 2) {
                omit += omitTimes[ i ] - omitTimes[ i - 1 ];
            }

            dfpair(out, "time to omit", "%f", omit);
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
        dfpair(out, "number of short trajectories taken", "%u", nshort);
        dfpair(out, "number of identity actions taken", "%u", nident);
        dfpair(out, "number of expanded states with identity actions", "%lu", identStates);
        dfpair(out, "number of current states with identity actions", "%lu", identCurs);
    }

private:

    // ExpandLss returns the cheapest  goal node if one was generated
    // and NULL otherwise.
    LssNode *expandLss(D &d, Node *rootNode) {

        //Clear OPEN and CLOSED
        lssOpen.clear();
        for (auto n : lssClosed)
            lssPool.destruct(n);
        lssClosed.clear();
        nclosed = 0;

        //Push the agent onto each list
        LssNode *a = lssPool.construct();
        a->node = rootNode;
        a->parent = NULL;
        a->op = D::Nop;
        a->g = 0;
        a->f = rootNode->h;
        a->genCount = 0;
        lssOpen.push(a);
        lssClosed.add(a);

        double herrnext = 0;
        double derrnext = 0;

        LssNode *goal = NULL;

        unsigned int exp = 0;
        while (!lssOpen.empty() && !lsslim->stop() && !this->limit()) {
            LssNode *s = *lssOpen.pop();

            nclosed += !s->closed;
            s->closed = true;
            s->expCount = ++exp;

            delaySum += s->expCount - s->genCount;
            ++delayCount;

            LssNode *bestkid = NULL;
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
                    kid->genCount = exp;
                    lssClosed.add(kid);
                }
                if (kid->g > s->g + e.outcost) {
                    if (kid->parent)    // !NULL if a dup
                        this->res.dups++;
                    kid->parent = s;
                    kid->g = s->g + e.outcost;
                    kid->f = kid->g + kid->node->h;

                    assert (derror != 1);
                    double d = kid->node->derr / (1 - derror);
                    double h = kid->node->h + herror*d;
                    kid->fhat = kid->g + h;

                    kid->op = e.op;
                    lssOpen.pushupdate(kid, kid->openind);
                }
                if (k->goal && (!goal || kid->g < goal->g))
                    goal = kid;

                if (!bestkid || kid->f < bestkid->f)
                    bestkid = kid;
            }

            if (bestkid) {
                double herr =  bestkid->f - s->f;
                if (herr < 0)
                    herr = 0;
                herrnext = herrnext + (herr - herrnext)/(exp+1);

                double derr = bestkid->node->d + 1 - s->node->d;
                if (derr < 0)
                    derr = 0;
                if (derr >= 1)
                    derr = 1 - geom2d::Threshold;
                derrnext = derrnext + (derr - derrnext)/(exp+1);
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

    void hCostLearning(D &d) {
        BinHeap<typename LssNode::HSort, LssNode*> open;

        open.append(lssOpen.data());

        std::vector<Node*> updated;

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
    }

    std::pair<std::pair<Node*, Node*>, std::pair<double, double>> move(D &d, Node *cur, Node *root, LssNode *goal, double leftOver) {

        //Search along the path from s to s'
        LssNode *best = goal;
        if (!best) {
            for (auto n : lssOpen.data()) {
                if (n->node == cur)
                    continue;
                if (best == NULL || LssNode::FHatSort::pred(n, best))
                    best = n;
            }
        }
        assert (best);
        assert (best->node != cur);
        std::vector<Oper> ops;

        LssNode *s = lssClosed.find(cur->state);
        std::pair<std::pair<std::pair<LssNode*, LssNode*>, double>, double> maxUtility = std::make_pair( std::make_pair( std::make_pair( best, best ), lsslim->lookahead(best->g)), 0.0 );

        std::pair<LssNode *, LssNode *> to = std::make_pair(best, best);
        assert (to.first->parent);
        omitTimes.push_back(cputime());

        LssNode *parent = best;
        LssNode *next = best->parent;
        best->parent = NULL;

        //Reversing the tree pointers
        for (LssNode *p = next; next != NULL; p = next)
        {
            next = p->parent;
            p->parent = parent;
            parent = p;
        }

        if (!goal) //no goal, find destination node.
        {

            LssNode *ident = NULL;

            if( s->node->ident.first == D::Nop )
                s = s -> parent;

            //ACTUALLY Search along the path from s to s'
            for ( ; s != best; s = s->parent)
            {
                //Do we have an identity action?  Great, then all future states have one too.
                if( s->node->ident.first != D::Nop )
                    ident = s;

                //If we are using prefix commitment, available quanta are all >= time to state
                    //We only have > options if we have an identity action at/before state
                //If we are using some other form of commitment we have quanta available < time to state
                //The best state is the one with highest Utility(state, quanta)

                //WHAT IS UTILITY?
                    //Benefit - Cost
                        //Benefit = average (currently) possible reduction in search time (currently reduction in f) * probability of reduction
                        //Cost = guaranteed loss (time spent not moving) + possible loss (bypassed possible reduction)

                //If the best state is utilizing an inherited identity action, commit to the state from which it inherits,
                    //but search at the state itself, and hold on to the remaining prefix

                std::pair<std::pair<std::pair<LssNode*, LssNode*>, double>, double> utility = computeUtility( d, s, ident, leftOver );

                //TODO account for bypassed utility?
                double waitTime = std::max( 0.0, utility.first.second - s->g - leftOver );
                utility.second -= waitTime;



                if( utility.second > maxUtility.second )
                {
                    //fprintf( stderr, "%f\n", utility.second );
                    maxUtility = utility;
                }
            }
        }

        parent = lssClosed.find(cur->state);
        next = parent->parent;
        parent->parent = NULL;

        //Reversing the tree pointers back
        for (LssNode *p = next; next != NULL; p = next)
        {
            next = p->parent;
            p->parent = parent;
            parent = p;
        }

        omitTimes.push_back(cputime());

        to = maxUtility.first.first;

        if( to.first->node == cur )
        {
            this->res.ops.push_back(cur->ident.first);
            lengths.push_back(1);
            nshort++;
            nident++;

            std::vector<Oper> temp = std::vector<Oper>();

            for (LssNode *p = to.second; p->node != root; p = p->parent) {
                assert (p->parent != to.second);   // no cycles
                temp.push_back(p->op);
            }

            P.insert(P.begin(), temp.begin(), temp.end());

            return std::make_pair(std::make_pair(cur, to.second->node), std::make_pair( cur->ident.second, leftOver + to.second->g ));
        }

        if (to.first != best)
            nshort++;

        for (LssNode *p = to.first; p->node != cur; p = p->parent) {
            assert (p->parent != to.first);   // no cycles
            ops.push_back(p->op);
        }

        ops.insert(ops.end(), P.begin(), P.end());

        P.clear();

        for (LssNode *p = to.second; p != to.first; p = p->parent) {
            assert (p->parent != to.second);   // no cycles
            P.push_back(p->op);
        }

        assert (ops.size() >= 1);

        this->res.ops.insert(this->res.ops.end(), ops.rbegin(), ops.rend());
        lengths.push_back(ops.size());
        return std::make_pair(std::make_pair(to.first->node, to.second->node), std::make_pair( to.first->g + leftOver, to.second->g - to.first->g ));
    }

    /*double gaussrand()
    {
        static double V1, V2, S;
        static int phase = 0;
        double X;

        if(phase == 0) {
            do {
                double U1 = (double)rand() / RAND_MAX;
                double U2 = (double)rand() / RAND_MAX;

                V1 = 2 * U1 - 1;
                V2 = 2 * U2 - 1;
                S = V1 * V1 + V2 * V2;
                } while(S >= 1 || S == 0);

            X = V1 * sqrt(-2 * log(S) / S);
        } else
            X = V2 * sqrt(-2 * log(S) / S);

        phase = 1 - phase;

        return X;
    }*/

    // Returns the best quanta/utility combo
    std::pair<std::pair<std::pair<LssNode*, LssNode*>, double>, double> computeUtility(D &d, LssNode *sPrimePrime, LssNode *ident, double leftOver) {
        BinHeap<typename LssNode::FHatSort, LssNode*> kids;

        Node *n = sPrimePrime->node;
        for (auto k : expand(d, n)) {
            auto kid = lssClosed.find(k.node->state);
            assert(kid);    // Must have been generated  since n was expanded.
            kids.push(kid);
        }

        assert (kids.size() > 0);

        if (kids.size() == 1)
            return std::make_pair(std::make_pair(std::make_pair(sPrimePrime, sPrimePrime), 0), 0.0);

        LssNode *alpha = *kids.pop();
        LssNode *beta = *kids.pop();

        //double belVarA = pow((alpha->fhat - alpha->f), 2t );
        double belVarB = pow((beta->fhat - beta->f), 2 );

        double quanta = 0.0;
        double lookahead;

        double expansionDelay = ((double)delaySum) / delayCount;

        double maxQuanta = sPrimePrime->g + leftOver;

        if( ident )
        {
            maxQuanta = std::max(maxQuanta, n->d / (lsslim->lookahead(1) / expansionDelay));
        }

        while( ++quanta < maxQuanta )
        {
            lookahead = lsslim->lookahead( quanta );

            //double futureVarA = belVarA * ( n->d / (lookahead / expansionDelay) );
            double futureVarB = belVarB * ( n->d / (lookahead / expansionDelay) );

            //double stdDevA = sqrt( futureVarA );
            double stdDevB = sqrt( futureVarB );

            if( beta->fhat - stdDevB + std::max(0.0, quanta - sPrimePrime->g) < alpha->fhat )
            {
                break;
            }
        }

        quanta = std::min(quanta, maxQuanta);

        double start = beta->fhat - sqrt(belVarB);
        double end = alpha->fhat;

        double benefit = (sqrt(belVarB * PI/2) * (alpha->fhat - beta->fhat) * erf((end-beta->fhat)/sqrt(2 * belVarB))+belVarB*exp(-(pow(end-beta->fhat,2))/(2 * belVarB)))/sqrt(2 * PI * belVarB);
        benefit -= (sqrt(belVarB * PI/2) * (alpha->fhat - beta->fhat) * erf((start-beta->fhat)/sqrt(2 * belVarB))+belVarB*exp(-(pow(start-beta->fhat,2))/(2 * belVarB)))/sqrt(2 * PI * belVarB);

        benefit /= 0.5 * erfc((beta->fhat - alpha->fhat)/sqrt(2 * belVarB));

        lookahead = lsslim->lookahead( quanta );
        //fprintf( stderr, "%f\n", lookahead );

        double futureVarB = belVarB * ( n->d / (lookahead / expansionDelay) );

        benefit *= 0.5 * erfc((beta->fhat - alpha->fhat)/sqrt(2 * futureVarB));

        LssNode *target = sPrimePrime;
        if( quanta > sPrimePrime->g + leftOver )
            target = ident;

        return std::make_pair( std::make_pair( std::make_pair( target, sPrimePrime ), quanta), benefit );
    }

    // Expand returns the successor nodes of a state.
    std::vector<outedge> expand(D &d, Node *n) {
        if (n->expd)
            return n->succs;

        State buf, &s = d.unpack(buf, n->state);

        this->res.expd++;

        n->ident = d.ident(s);
        if (n->ident.first != D::Nop)
            identStates++;

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

    BinHeap<typename LssNode::FHatSort, LssNode*> lssOpen;
    ClosedList<typename LssNode::Nodes, LssNode, D> lssClosed;
    Pool<LssNode> lssPool;
    unsigned int nclosed;

    double herror;
    double derror;

    LookaheadLimit *lsslim;
    Nodes nodes;

    unsigned long delaySum, delayCount;

    unsigned int nshort, nident;    // Number of short and identity trajectories taken.
    unsigned long identStates;  // Number of expanded states with identity actions.
    unsigned long identCurs;    // Number of 'current states' with identity actions.

    std::vector<Oper> P;

    std::vector<double> times;
    std::vector<double> cputimes;
    std::vector<double> omitTimes;
    std::vector<unsigned int> lengths;
};
