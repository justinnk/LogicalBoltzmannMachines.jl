
using LogicalBoltzmannMachines

"Reproduces the xor example from the paper."
function tweety_example()
    f = DNFFormula([
        # t ∧ b
        ConjunctiveClause([Lit('t'), Lit('b'), Lit('f')]),
        ConjunctiveClause([Lit('t'), Lit('b'), nLit('f')]),

        # ¬t
        ConjunctiveClause([nLit('t'), nLit('b'), Lit('f')]),
        ConjunctiveClause([nLit('t'), Lit('b'), nLit('f')]),
        ConjunctiveClause([nLit('t'), nLit('b'), nLit('f')]),

        # b ∧ f
        ConjunctiveClause([Lit('b'), Lit('f'), nLit('t')]),

        # ¬b
        ConjunctiveClause([nLit('b'), Lit('t'), Lit('f')]),
        ConjunctiveClause([nLit('b'), Lit('t'), nLit('f')]),
    ])
    @show isfull(f)
    lbm, lit2idx, idx2lit = sdnf2lbm(f)
    query = Dict(lit2idx['t'] => 1.0, lit2idx['f'] => 1.0)
    res = reason(lbm; samples=Int(1e2), query=query)
    print(query, "\n")
    for (key, val) in query
        print(idx2lit[key], " is ", val, "\n")
    end
    print("???\n")
    if isempty(res)
        print("No.")
    else
        print("Yes.\n")
        print("(In these cases:\n")
        for (i, _) in enumerate(collect(res))
            for (idx, val) in enumerate(collect(res)[i])
                print(idx2lit[idx], " is ", val, "\n")
            end
            print("OR\n")
        end
        print(")")
    end
end

tweety_example()