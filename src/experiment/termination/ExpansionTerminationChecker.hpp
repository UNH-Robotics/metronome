#ifndef METRONOME_EXPANSIONTERMINATIONCHECKER_HPP
#define METRONOME_EXPANSIONTERMINATIONCHECKER_HPP

namespace metronome {

class ExpansionTerminationChecker {
public:
    ExpansionTerminationChecker() {}

    ExpansionTerminationChecker(const ExpansionTerminationChecker&) = delete;

    ExpansionTerminationChecker(ExpansionTerminationChecker&&) = delete;

    void resetTo(unsigned int expansionLimit) {
        this->expansionLimit = expansionLimit;
        expansionCount = 0;
    }

    ExpansionTerminationChecker& limit(double limit) {
        assert(limit > 0 && limit <= 1);
        this->limitRatio= limit;
        return *this;
    }

    bool reachedTermination() const {
        return expansionCount >= expansionLimit * limitRatio;
    }

    /**
     * Should be called by the planner at every expansion.
     */
    void notifyExpansion() {
        ++expansionCount;
    }

    unsigned int expansionsPerAction(unsigned long long actionDuration) const {
        return static_cast<unsigned int>(actionDuration);
    }
    
private:
    unsigned int expansionCount{0};
    unsigned int expansionLimit{0};
    double limitRatio{1};
};
}

#endif //METRONOME_EXPANSIONTERMINATIONCHECKER_HPP
