
using LogicalBoltzmannMachines

"Reproduces the xor example from the paper."
function xor_example()
    f = DNFFormula([
        ConjunctiveClause([nLit("x"), nLit("y"), nLit("z")]),
        ConjunctiveClause([nLit("x"),  Lit("y"),  Lit("z")]),
        ConjunctiveClause([ Lit("x"), nLit("y"),  Lit("z")]),
        ConjunctiveClause([ Lit("x"),  Lit("y"), nLit("z")]),
    ])
    @show isfull(f)
    @show dnf2wff(f)
    @show lbm, _, _ = sdnf2lbm(f)
    # plot_fenergy_space(lbm, 5.0, 0.5)
    plot_energy_space!(lbm)
    res = reason!(lbm; samples=Int(1e2))
    @show res
end

@time xor_example()
