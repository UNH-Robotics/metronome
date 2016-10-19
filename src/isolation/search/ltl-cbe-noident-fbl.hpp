#pragma once
#include "search.hpp"
#include "../utils/geom2d.hpp"
#include "../utils/pool.hpp"
#include "lsslrtastar2.hpp"
#include <vector>

#define PI 3.14159265359

template <class D>
class LTL_CBE_NOIDENT_FBL : public SearchAlgorithm<D> {
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
        bool expd; // Was this expanded before?
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
        LssNode *parent, *best;
        unsigned long genCount, expCount, depth; //Counters for node generation and expansion.
        double g, f, fhat, herr;
        Oper op;
        long openind;
        bool updated;
        bool closed;

    private:
        ClosedEntry<LssNode, D> nodesent;
    };

public:

    LTL_CBE_NOIDENT_FBL(int argc, const char *argv[]) :
        SearchAlgorithm<D>(argc, argv),
        lssClosed(4051),
        herror(0),
        derror(0),
        nodes(30000001),
        onestep(false),
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

    ~LTL_CBE_NOIDENT_FBL() {
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

        lsslim->start(0);

        while (!cur->goal && !this->limit()) {

            if (cur->ident.first != D::Nop)
                identCurs++;

            //fprintf( stderr, "%lu:%lu\n", this->res.expansions(), this->res.generations() );

            LssNode *goal = expandLss(d, cur);
            if (this->limit())
                break;
            if (!goal)
                hCostLearning(d);
            auto m = move(d, cur, goal);
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

    virtual void output(FILE *out) {
        SearchAlgorithm<D>::output(out);
        dfpair(out, "num steps", "%lu", (unsigned long) times.size());
        assert (lengths.size() == times.size());
        //dfpair(out, "h error last", "%g", herror);
        //dfpair(out, "d error last", "%g", derror);
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

        delaySum = 0.0;
        delayCount = 0;

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
        a->herr = 0.0;
        a->depth = 0;
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

                    kid->depth = s->depth + 1;

                    if (kid->parent)    // !NULL if a dup
                        this->res.dups++;
                    kid->parent = s;
                    kid->g = s->g + e.outcost;
                    kid->f = kid->g + kid->node->h;

                    assert (derror != 1);
                    double d = kid->node->derr / (1 - derror);
                    double h = kid->node->h + herror*d;
                    //kid->herr = std::max( 0.0, kid->f - s->f );
                    //kid->herr = s->herr + ( kid->herr - s->herr ) / kid->depth;
                    kid->fhat = kid->g + h;//kid->node->h + kid->herr * kid->node->d;

                    kid->op = e.op;
                    lssOpen.pushupdate(kid, kid->openind);
                }
                if (k->goal && (!goal || kid->g < goal->g))
                    goal = kid;

                if (!bestkid || kid->f < bestkid->f)
                    bestkid = kid;
            }

            if (bestkid) {
                s->best = bestkid;

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

        for( LssNode *s : open.data( ) )
            s->best = s;

        while (nclosed > 0 && !open.empty()) {
            LssNode *s = *open.pop();

            nclosed -= s->closed;

            for (auto e : s->node->preds) {
                Node *sprime = e.node;

                LssNode *sp = lssClosed.find(sprime->state);
                if (!sp || !sp->closed)
                    continue;

                if( sp->best == s )
                {
                    sp->best = s->best;
                }

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

    std::pair<Node*, double> move(D &d, Node *cur, LssNode *goal) {

        //Search along the path from s to s'
        LssNode *best = goal;
        LssNode *secondBest = goal;
        if (!best) {
            for (auto n : lssOpen.data()) {
                if (n->node == cur)
                    continue;
                if (best == NULL)
                    best = n;
                else if(LssNode::FSort::pred(n, best))
                {
                    secondBest = best;
                    best = n;
                }
                else if( secondBest == NULL || LssNode::FSort::pred(n, secondBest) )
                    secondBest = n;
            }
        }
        assert (best);
        assert (best->node != cur);
        std::vector<Oper> ops;

        LssNode *to = best;
        assert (to->parent);
        omitTimes.push_back(cputime());

        if( !goal )
        {
            LssNode *alpha, *beta;

            std::vector<LssNode *> path;

            for( LssNode *p = best; p->node != cur; p = p->parent )
                path.push_back( p );

            do
            {
                alpha = NULL;
                beta = NULL;
                to = path.back( );
                path.pop_back( );

                if( to == best )
                    break;

                for (auto n : lssOpen.data()) {

                    bool leaf = false;
                    for( LssNode *p = n; p->node != cur; p = p->parent )
                    {
                        if( p == to )
                        {
                            leaf = true;
                            break;
                        }
                    }

                    if( !leaf )
                        continue;

                    if (alpha == NULL)
                        alpha = n;
                    else if(LssNode::FSort::pred(n, alpha))
                    {
                        beta = alpha;
                        alpha = n;
                    }
                    else if( beta == NULL || LssNode::FSort::pred(n, beta) )
                        beta = n;
                }
            }
            //TODO doublecheck rhs if this doesn't work
            while(computeBenefit(d, cur, alpha, beta, to->g) <= 0.0);

            if( to!= best )
                ++nshort;

        }

        omitTimes.push_back(cputime());

        //to = maxUtility.first.first;

        /*if( to->node == cur )
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
            nshort++;*/

        for (LssNode *p = to; p->node != cur; p = p->parent) {
            assert (p->parent != to);   // no cycles
            ops.push_back(p->op);
        }

        /*ops.insert(ops.end(), P.begin(), P.end());

        P.clear();

        for (LssNode *p = to.second; p != to.first; p = p->parent) {
            assert (p->parent != to.second);   // no cycles
            P.push_back(p->op);
        }*/

        assert (ops.size() >= 1);

        this->res.ops.insert(this->res.ops.end(), ops.rbegin(), ops.rend());
        lengths.push_back(ops.size());
        return std::make_pair(to->node, to->g);
    }

    // Returns the utility of searching at a given state
    double computeBenefit(D &d, Node *n, LssNode *alpha, LssNode *beta, double g) {

        if( beta == NULL || alpha == beta )
            return 0.0;

        double x = lsslim->lookahead(g);

        int i = 0;
        double delay_alpha = 0.0;
        for (LssNode *p = alpha->parent; p->node != n && delay_alpha < x; p = p->parent) {
            ++i;
            delay_alpha += p->expCount - p->genCount;
        }
        delay_alpha /= i;
        i = 0;
        double delay_beta = 0.0;
        for (LssNode *p = beta->parent; p->node != n && delay_beta < x; p = p->parent) {
            ++i;
            delay_beta += p->expCount - p->genCount;
        }
        delay_beta /= i;

        double mu_alpha = alpha->fhat, mu_beta = beta->fhat;
        double bv_alpha = pow(alpha->fhat - alpha->f, 2),
                bv_beta  = pow(beta->fhat - beta->f, 2);
        double x_alpha = lsslim->lookahead(g) / delay_alpha,
                x_beta = lsslim->lookahead(g) / delay_beta;

        //TODO if not only considering identity actions, fix this
        double v_alpha = bv_alpha * (x_alpha / alpha->node->d), v_beta = bv_beta * (x_beta / beta->node->d);

        double start = ( mu_alpha - 2.0 * sqrt(v_alpha) < mu_beta - 2.0 * sqrt(v_beta)) ? mu_alpha - 2.0 * sqrt(v_alpha) : mu_beta - 2.0 * sqrt(v_beta);
        double end = ( mu_alpha + 2.0 * sqrt(v_alpha) > mu_beta + 2.0 * sqrt(v_beta)) ? mu_alpha + 2.0 * sqrt(v_alpha) : mu_beta + 2.0 * sqrt(v_beta);
        double a, b;
        double benefit = 0.0;

        for( int i = 0; i < 100; ++i )
        {
            a = start + ((end - start) * ((i + 0.5)/100.0));

            int j;
            double sum = 0.0;
            for( j = 0; j < 100; ++j )
            {
                b = start + ((end - start) * ((j + 0.5)/100.0));

                if( a < b )
                    break;

                sum += (a - b) * exp(-pow(b-mu_beta,2)/(2 * v_beta))/(sqrt(2 * PI * v_beta));
            }

            sum *= (end - start)/100.0;

            benefit += sum * exp(-pow(a-mu_alpha,2)/(2 * v_alpha))/(sqrt(2 * PI * v_alpha)) * ((end - start)/100.0);
        }

        return benefit;
    }

    // Expand returns the successor nodes of a state.
    std::vector<outedge> expand(D &d, Node *n) {
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
            Node *k = nodes.get(d, e.state);
            k->preds.emplace_back(n, e.cost);
            k->h = std::max(k->h, n->h - e.cost);
            n->succs.emplace_back(k, ops[i], e.revcost, e.cost);
        }
        n->expd = true;

        return n->succs;
    }

    BinHeap<typename LssNode::FSort, LssNode*> lssOpen;
    ClosedList<typename LssNode::Nodes, LssNode, D> lssClosed;
    Pool<LssNode> lssPool;
    unsigned int nclosed;

    double herror;
    double derror;

    LookaheadLimit *lsslim;
    Nodes nodes;

    bool onestep;

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
