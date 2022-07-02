
using LogicalBoltzmannMachines

"Reproduces the xor example from the paper."
function tweety_example()
    # tweety -> bird
    # bird -> flies
    f = DNFFormula([
        # t ∧ b
        ConjunctiveClause([Lit("tweety"), Lit("bird(tweety)")]),
        # ¬
        ConjunctiveClause([nLit("tweety")]),
        # b ∧ f
        ConjunctiveClause([Lit("bird(tweety)"), Lit("flies(tweety)")]),
        # ¬ b
        ConjunctiveClause([nLit("bird(tweety)")])
    ])

    # convert f to an LBM
    lbm, lit2idx, idx2lit = sdnf2lbm(f; ϵ=0.5)

    # perform reasoning
    query = Dict{Int64, Float64}()
    query = Dict{Int64, Float64}(lit2idx["tweety"] => 1.0, lit2idx["flies(tweety)"] => 1.0)
    res = reason!(lbm; samples=Int(1e1), query=query, ϵ=0.5, factor=2)

    # visualisation and fancy output
    plot_energy_space!(lbm)
    # plot_fenergy_space!(lbm, 5.0, 0.5)
    @show lit2idx

    # visualize the underlying RBM
    open("tweety_example.dot", "w+") do io
       write(io, lbm2dot(lbm, idx2lit)) 
    end
    run(
        pipeline(
            `dot -Tpdf -Kdot tweety_example.dot`,
            stdout="tweety_example.pdf"
        )
    )

    # human-readable output
    if isempty(query)
        println("Find all satisfying valuations!")
    else
        for (key, val) in query
            print(idx2lit[key], " is ", val, "\n")
        end
        println("???")
    end
    if isempty(res)
        println("No.")
    else
        println("Yes, in these cases:")
        for (i, _) in enumerate(collect(res))
            for (idx, val) in enumerate(collect(res)[i])
                print(idx2lit[idx], " is ", val, "\n")
            end
            println("OR")
        end
    end

end

tweety_example()