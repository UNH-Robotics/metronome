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

            static bool comp(LssNode *a, LssNode *b) {
                return !FHatSort::pred(a, b);
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

    class TLA {
        public:
        TLA( LssNode *s ) : tla( s )
        {
            leaves = std::vector<LssNode *>( 1, s );
        }

        static bool comp( TLA *a, TLA *b )
        {
            return LssNode::FHatSort::comp( a->leaves.front( ), b->leaves.front( ) );
        }

        static void clear( TLA *tla )
        {
            tla->leaves.clear( );
        }

        static std::vector<LssNode *> combined( void )
        {
            std::vector<LssNode *> ret = std::vector<LssNode *>( );

            for( TLA *t : tlas )
            {
                for( LssNode *s : t->leaves )
                {
                    bool match = false;
                    for( LssNode *sPrime : ret )
                    {
                        if( s == sPrime )
                        {
                            match = true;
                            break;
                        }
                    }

                    if( !match )
                        ret.push_back( s );
                }
            }

            std::make_heap( ret.begin( ), ret.end( ), LssNode::FHatSort::comp );

            long count = 0;
            for( LssNode *s : ret )
            {
                s->openind = count++;
            }

            return ret;
        }

        static TLA * best_global( void )
        {
            TLA *best = tlas.front( );

            double alpha_reduction = 0.0;
            double best_reduction = 0.0;

            auto alpha_pair = tlas.front()->best_two( );
            LssNode *alpha1 = alpha_pair.first, *alpha2 = alpha_pair.second;
            double mean_alpha = alpha1->fhat, vertex_alpha = alpha2 ? alpha2->fhat : INFINITY;
            double variance_alpha = pow( ( alpha1->fhat - alpha1->f ), 2 );

            for( auto it = tlas.begin( ) + 1; it != tlas.end( ); ++it )
            {
                if( (*it)->front( ) == alpha1 )
                    (*it)->pop( );

                auto t_pair = (*it)->best_two( );
                LssNode *t1 = t_pair.first, *t2 = t_pair.second;
                double mean_t = t1->fhat, vertex_t = t2 ? t2->fhat : INFINITY;
                double variance_t = pow( ( t1->fhat - t1->f  ), 2 );

                double alpha_start = mean_alpha - 2 * sqrt( variance_alpha );
                double alpha_end = mean_alpha + 2 * sqrt( variance_alpha );

                double a = exp((-(pow(alpha_start - mean_alpha, 2)/(2 * variance_alpha)))/sqrt(2 * PI * variance_alpha));
                double sum = 0.0;
                for( int i = 1; i < 101; ++i )
                {
                    double x = alpha_start + ( ( alpha_end - alpha_start ) / 100 * i );
                    double b = exp((-(pow(x - mean_alpha, 2)/(2 * variance_alpha)))/sqrt(2 * PI * variance_alpha));

                    if( x < vertex_alpha )
                        sum += ( mean_t - x ) * ( ( alpha_end - alpha_start ) / 200 ) * ( b + a );
                    else
                        sum += ( mean_t - vertex_alpha ) * ( ( alpha_end - alpha_start ) / 200 ) * ( b + a );

                    a = b;
                }

                alpha_reduction += sum;
                if( alpha_reduction > best_reduction )
                {
                    best_reduction = alpha_reduction;
                    best = tlas.front( );
                }

                double t_start = mean_t - 2 * sqrt( variance_t );
                double t_end = mean_t + 2 * sqrt( variance_t );

                a = exp((-(pow(t_start - mean_t, 2)/(2 * variance_t)))/sqrt(2 * PI * variance_t));
                sum = 0.0;
                for( int i = 1; i < 101; ++i )
                {
                    double x = t_start + ( ( t_end - t_start ) / 100 * i );
                    double b = exp((-(pow(x - mean_t, 2)/(2 * variance_t)))/sqrt(2 * PI * variance_t));

                    if( x < vertex_t )
                        sum += ( x - mean_alpha ) * ( ( t_end - t_start ) / 200 ) * ( b + a );
                    else
                        sum += ( vertex_t - mean_alpha ) * ( ( t_end - t_start ) / 200 ) * ( b + a );

                    a = b;
                }

                if( sum > best_reduction )
                {
                    best_reduction = sum;
                    best = *it;
                }
            }

            return best;
        }

        LssNode * front( void )
        {
            LssNode *s = pop( );

            if( s )
                push( s );

            return s;
        }

        std::pair<LssNode *, LssNode *> best_two( void )
        {
            LssNode *first = pop( );
            LssNode *second = front( );
            push( first );

            return std::make_pair( first, second );
        }

        void push( LssNode *s )
        {
            for( auto it = leaves.begin( ); it != leaves.end( ); ++it )
            {
                if( *it == s )
                {
                    leaves.erase( it );
                    leaves.push_back( s );
                    std::make_heap( leaves.begin( ), leaves.end( ), LssNode::FHatSort::comp );
                    return;
                }
            }

            leaves.push_back( s );
            std::push_heap( leaves.begin( ), leaves.end( ), LssNode::FHatSort::comp );
        }

        LssNode * pop( void )
        {
            if( leaves.empty( ) )
                return NULL;

            LssNode *s;

            do
            {
                std::pop_heap( leaves.begin( ), leaves.end( ), LssNode::FHatSort::comp );
                s = leaves.back( );
                leaves.pop_back( );
            }
            while( !leaves.empty( ) && s && s->closed );

            return s;
        }

        private:
        std::vector<LssNode *> leaves;
        LssNode *tla;
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
        tlas.clear();
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
        for( TLA *t : tlas )
        {
            TLA::clear( t );
        }
        tlas.clear();
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
        lssClosed.add(a);

        LssNode *goal = NULL;

        double herrnext = 0;
        double derrnext = 0;


        unsigned int exp = 0;

        for ( auto e : expand(d, a->node) )
        {
            if( e.node == rootNode )
                continue;

            LssNode *kid = lssPool.construct();
            kid->parent = a;
            kid->node = e.node;
            kid->g = e.outcost;
            lssClosed.add(kid);
            kid->f = kid->g + kid->node->h;

            double d = kid->node->derr / (1 - derror);
            double h = kid->node->h + herror*d;
            kid->fhat = kid->g + h;

            kid->op = e.op;

            tlas.push_back( new TLA( kid ) );
            std::push_heap( tlas.begin(), tlas.end( ), TLA::comp );
        }

        while ( !lsslim->stop() && !this->limit() )
        {
            LssNode *s;

            TLA *t = TLA::best_global( );

            s = t->front( );

            nclosed += !s->closed;
            s->closed = true;
            s->expCount = ++exp;

            delaySum += s->expCount - s->genCount;
            ++delayCount;

            LssNode *bestkid = NULL;

            for( auto e : expand(d, s->node) )
            {
                Node *k = e.node;
                if (s->parent && k == s->parent->node)
                    continue;

                LssNode *kid = lssClosed.find(k->state);

                if (!kid) {
                    kid = lssPool.construct();
                    kid->node = k;
                    kid->g = geom2d::Infinity;
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
                    t->push(kid);
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
                t->push(s);
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

            open.append(TLA::combined( ));

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

    std::pair<Node*, double> move(D &d, Node *cur, LssNode *goal) {

        //Search along the path from s to s'
        LssNode *best = goal;
        if (!best) {
            best = tlas.front( )->front( );
        }
        assert (best);
        assert (best->node != cur);
        std::vector<Oper> ops;

        LssNode *to = best;
        assert (to->parent);
        omitTimes.push_back(cputime());

        //Reversing the tree pointers
        /*LssNode *parent = best;
        LssNode *next = best->parent;
        best->parent = NULL;
        for (LssNode *p = next; next != NULL; p = next)
        {
            next = p->parent;
            p->parent = parent;
            parent = p;
        }*/

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

    static std::vector<TLA *> tlas;

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

template <>
std::vector<LTL_CBE<GridNav>::TLA *> LTL_CBE<GridNav>::tlas = std::vector<LTL_CBE<GridNav>::TLA *>();
