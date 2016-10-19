#include "search.hpp"
//#include "idastar.hpp"
//#include "astar.hpp"
//#include "astar-dump.hpp"
//#include "wastar.hpp"
//#include "greedy.hpp"
//#include "bugsy.hpp"
//#include "bugsy-slim.hpp"
//#include "arastar.hpp"
//#include "gstar.hpp"
//#include "fhatlrtastar.hpp"
//#include "dtastar-dump.hpp"
//#include "dtastar.hpp"
//#include "multilrtastar.hpp"
//#include "lsslrtastar-rand.hpp"
//#include "frit-offline.hpp"
//#include "dalsslrtastar.hpp"
//#include "wlrtastar.hpp"
//#include "flrtastar.hpp"
//#include "rtastar.hpp"
#include "lsslrtastar2.hpp"
#include "fhatlrtastar2.hpp"
#include "fhatlrtastar3.hpp"
#include "ltl-cbe-ident.hpp"
#include "ltl-cbe-ident-fbl.hpp"
#include "ltl-cbe-ident-hybrid.hpp"
#include "ltl-cbe-ident-hybrid-fbl.hpp"
#include "ltl-cbe-noident.hpp"
#include "ltl-cbe-noident-fbl.hpp"
#include "ltl-cbe-noident-hybrid.hpp"
#include "ltl-cbe-noident-hybrid-fbl.hpp"
#include "morts.hpp"
//#include "distributions.hpp"
//#include "oracle-ident.hpp"

/*
#include "monstar-dump.hpp"
#include "lsslrtastar.hpp"
#include "flrtastar2.hpp"
#include "cautiouslrtastar.hpp"
#include "frankenlrtastar.hpp"
#include "dflrtastar.hpp"
#include "mrastar.hpp"
#include "greedylrtastar.hpp"
#include "uclrtastar.hpp"
*/

#include <cstddef>
#include <cstdio>

// Functions for conveniently defining a new main
// function for a domain.

void dfheader(FILE*);
void dffooter(FILE*);
void dfpair(FILE*, const char*, const char*, ...);
void dfprocstatus(FILE*);
void fatal(const char*, ...);

template<class D>
SearchAlgorithm<D>* getsearch(int argc, const char* argv[]);

template<class D>
Result<D> search(D& d, int argc, const char* argv[]) {
    return searchGet(getsearch, d, argc, argv);
}

template<class D>
Result<D> searchGet(SearchAlgorithm<D>* (* get)(int, const char* []), D& d, int argc, const char* argv[]) {
    SearchAlgorithm<D>* srch = get(argc, argv);
    if (!srch && argc > 1)
        fatal("Unknown search algorithm: %s", argv[1]);
    if (!srch)
        fatal("Must specify a search algorithm");

    bool printSolution = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-printsolution") == 0)
            printSolution = true;
    }

    typename D::State s0 = d.initialstate();
    dfpair(stdout, "initial heuristic", "%f", (double) d.h(s0));
    dfpair(stdout, "initial distance", "%f", (double) d.d(s0));
    dfpair(stdout, "algorithm", argv[1]);

    try {
        srch->search(d, s0);
    } catch (std::bad_alloc&) {
        dfpair(stdout, "out of memory", "%s", "true");
        srch->res.path.clear();
        srch->res.ops.clear();
        srch->finish();
    }
    if (srch->res.path.size() > 0) {
        dfpair(stdout, "final sol cost", "%f",
               (double) d.pathcost(srch->res.path, srch->res.ops, printSolution));
    } else {
        dfpair(stdout, "final sol cost", "%f", -1.0);
    }
    srch->output(stdout);

    Result<D> res = srch->res;
    delete srch;

    return res;
}

template<class D>
SearchAlgorithm<D>* getsearch(int argc, const char* argv[]) {
    if (argc < 2)
        fatal("No algorithm specified");

//	if (strcmp(argv[1], "idastar") == 0)
//		return new Idastar<D>(argc, argv);
//	else if (strcmp(argv[1], "astar") == 0)
//		return new Astar<D>(argc, argv);
//	else if (strcmp(argv[1], "astar-dump") == 0)
//		return new Astar_dump<D>(argc, argv);
//	else if (strcmp(argv[1], "wastar") == 0)
//		return new Wastar<D>(argc, argv);
//	else if (strcmp(argv[1], "greedy") == 0)
//		return new Greedy<D>(argc, argv);
//	else if (strcmp(argv[1], "speedy") == 0)
//		return new Greedy<D, true>(argc, argv);
//	else if (strcmp(argv[1], "bugsy") == 0)
//		return new Bugsy<D>(argc, argv);
//	else if (strcmp(argv[1], "bugsy-slim") == 0)
//		return new Bugsy_slim<D>(argc, argv);
//	else if (strcmp(argv[1], "arastar") == 0)
//		return new Arastar<D>(argc, argv);
//	else if (strcmp(argv[1], "arastarmon") == 0)
//		return new ArastarMon<D>(argc, argv);
//	else if (strcmp(argv[1], "dtastar-dump") == 0)
//		return new Dtastar_dump<D>(argc, argv);
//	else if (strcmp(argv[1], "lsslrtastar2") == 0)
//		return new Lsslrtastar2<D>(argc, argv);
//    else if (strcmp(argv[1], "gstar") == 0)
//        return new GStar<D>(argc, argv);
//	else if (strcmp(argv[1], "fhatlrtastar") == 0)
//		return new Fhatlrtastar<D>(argc, argv);
//    else if (strcmp(argv[1], "fhatlrtastar2") == 0)
//        return new Fhatlrtastar2<D>(argc, argv);
//    else if (strcmp(argv[1], "fhatlrtastar3") == 0)
//        return new Fhatlrtastar3<D>(argc, argv);
//	else if (strcmp(argv[1], "dtastar") == 0)
//		return new Dtastar<D>(argc, argv);
//	else if (strcmp(argv[1], "multilrtastar") == 0)
//		return new Multilrtastar<D>(argc, argv);
//	else if (strcmp(argv[1], "dalsslrtastar") == 0)
//		return new daLsslrtastar<D>(argc, argv);
//    else if (strcmp(argv[1], "wlrtastar") == 0)
//        return new Wlrtastar<D>(argc, argv);
//    else if (strcmp(argv[1], "flrtastar") == 0)
//        return new Flrtastar<D>(argc, argv);
//    else if (strcmp(argv[1], "rtastar") == 0)
//		return new Rtastar<D>(argc, argv);
    if (strcmp(argv[1], "ltl-cbe-ident") == 0)
        return new LTL_CBE_IDENT<D>(argc, argv);
    else if (strcmp(argv[1], "ltl-cbe-ident-fbl") == 0)
        return new LTL_CBE_IDENT_FBL<D>(argc, argv);
    else if (strcmp(argv[1], "ltl-cbe-ident-hybrid") == 0)
        return new LTL_CBE_IDENT_HYBRID<D>(argc, argv);
    else if (strcmp(argv[1], "ltl-cbe-ident-hybrid-fbl") == 0)
        return new LTL_CBE_IDENT_HYBRID_FBL<D>(argc, argv);
    else if (strcmp(argv[1], "ltl-cbe-ident-fbl") == 0)
        return new LTL_CBE_IDENT_FBL<D>(argc, argv);
    else if (strcmp(argv[1], "ltl-cbe-noident") == 0)
        return new LTL_CBE_NOIDENT<D>(argc, argv);
    else if (strcmp(argv[1], "ltl-cbe-noident-fbl") == 0)
        return new LTL_CBE_NOIDENT_FBL<D>(argc, argv);
    else if (strcmp(argv[1], "ltl-cbe-noident-hybrid") == 0)
        return new LTL_CBE_NOIDENT_HYBRID<D>(argc, argv);
    else if (strcmp(argv[1], "ltl-cbe-noident-hybrid-fbl") == 0)
        return new LTL_CBE_NOIDENT_HYBRID_FBL<D>(argc, argv);
    else if (strcmp(argv[1], "morts") == 0)
        return new MORTS<D>(argc, argv);
//    else if (strcmp(argv[1], "distributions") == 0)
//        return new Distributions<D>(argc, argv);
//    else if (strcmp(argv[1], "oracle-ident") == 0)
//        return new Oracle_Ident<D>(argc, argv);
//	else if (strcmp(argv[1], "frit-offline") == 0)
//		return new Frit_Offline<D>(argc, argv);

/*
	else if (strcmp(argv[1], "lsslrtastarrandom") == 0)
		return new Lsslrtastar2Random<D>(argc, argv);
	else if (strcmp(argv[1], "monstar-dump") == 0)
		return new Monstar_dump<D>(argc, argv);

	else if (strcmp(argv[1], "mrastar") == 0)
		return new Mrastar<D>(argc, argv);
	else if (strcmp(argv[1], "dflrtastar") == 0)
		return new Dflrtastar<D>(argc, argv);
	else if (strcmp(argv[1], "greedylrtastar") == 0)
		return new Greedylrtastar<D>(argc, argv);
	else if (strcmp(argv[1], "uclrtastar") == 0)
		return new Uclrtastar<D>(argc, argv);
	else if (strcmp(argv[1], "cautiouslrtastar") == 0)
		return new Cautiouslrtastar<D>(argc, argv);
	else if (strcmp(argv[1], "frankenlrtastar") == 0)
		return new Frankenlrtastar<D>(argc, argv);
	else if (strcmp(argv[1], "lsslrtastar") == 0)
		return new Lsslrtastar<D>(argc, argv);
	else if (strcmp(argv[1], "flrtastar2") == 0)
		return new Flrtastar2<D>(argc, argv);
*/

    fatal("Unknown algorithm: %s", argv[1]);
    return NULL;    // Unreachable
}
