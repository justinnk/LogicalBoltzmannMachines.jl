
using LogicalBoltzmannMachines

function main()
    f = DNFFormula([
        ConjunctiveClause([nLit('x'), nLit('y'), nLit('z')]),
        ConjunctiveClause([nLit('x'),  Lit('y'),  Lit('z')]),
        ConjunctiveClause([ Lit('x'), nLit('y'),  Lit('z')]),
        ConjunctiveClause([ Lit('x'),  Lit('y'), nLit('z')]),
    ])
    @show is_full(f)
    @show dnf2wff(f)
    @show lbm_vals = sdnf2lbm(f)
    @show find_sat(lbm_vals...)
end

main()