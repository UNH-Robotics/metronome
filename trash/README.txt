The search algorithms in question are all in search/search.  morts.hpp should be the main one with the variants all named with the format:
    ltl-cbe (Least time lost, cost based lookahead I believe is what I was working with)
    -ident or -noident depending on whether or not identity actions were considered (not sure why I didn't just use an argument for that)
    -hybrid to denote using the path based error model for computing belief distributions
    -fbl to denote a f-based lookahead rather than fhat based (I believe my experimentation showed no significant difference there)

