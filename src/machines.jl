
using Plots

"A Logical Boltzmann Machine."
mutable struct LBM
    "Values of the input neurons."
    x::Vector{Float64}
    "Bias for each input neuron."
    a::Vector{Float64}
    "Values of the hidden neurons."
    h::Vector{Float64}
    "Bias for each hidden neuron."
    b::Vector{Float64}
    "Weights between h (rows) and x (columns)."
    W::Matrix{Float64}
    LBM(
        x::Vector{Float64},
        a::Vector{Float64},
        h::Vector{Float64},
        b::Vector{Float64},
        W::Matrix{Float64}
    ) = new(x, a, h, b, W)
end

"Creates a Restricted Boltzmann Machine with random visible and hidden values."
function LBM(W::Matrix{Float64}, a::Vector{Float64}, b::Vector{Float64})
    size(a, 1) == size(W, 2) || error(
        "Bias vector for the visible layer is incompatible with weight Matrix."
    )
    size(b, 1) == size(W, 1) || error(
        "Bias vector for the hidden layer is incompatible with weight Matrix."
    )
    LBM(rand([0.0, 1.0], size(a, 1)), a, rand([0.0, 1.0], size(b, 1)), b, W)
end

"Initialize the visible layer of rbm to a random vector."
rand_visible!(rbm::LBM) = rbm.x = rand([0.0, 1.0], size(rbm.a, 1))
"Initialize the visible layer to a random vector."
rand_hidden!(rbm::LBM) = rbm.h = rand([0.0, 1.0], size(rbm.b, 1))

"Calculate the energy of rbm the provided LBM."
energy(rbm::LBM) = -(rbm.W * rbm.x)' * rbm.h - rbm.a' * rbm.x - rbm.b' * rbm.h
"Calculate the free energy of the provided LBM."
fenergy(rbm::LBM, c::Float64) = 
    -sum(log.(1 .+ exp.(c *(rbm.W * rbm.x + rbm.b)))) - rbm.a' * rbm.x

"The sigmoid of x, scaled by T^-1."
sigmoid(x::Float64, T::Float64) = (1.0 / (1.0 + exp(-(1/T) * x)))
"The probability of hidden layer values given visible layer values, P(h|x)."
phx(W::Matrix{Float64}, x::Vector{Float64}, b::Vector{Float64}, T::Float64) =
    round.(sigmoid.(W * x + b, T))
"The probability of visible layer values given hidden layer values, P(x|h)."
pxh(W::Matrix{Float64}, h::Vector{Float64}, a::Vector{Float64}, T::Float64) =
    round.(sigmoid.(W' * h + a, T))
"Set visible layer values for the given RBM at certain indices."
function set_visible_values!(rbm::LBM, values::Dict{Int64, Float64})
    for (idx, f) in values
        rbm.x[idx] = f
    end
    rbm
end
"Perform Gibbs sampling on the given RBM, while fixating some visible values."
function gibbs!(rbm::LBM, fixed::Dict{Int64, Float64}=Dict{Int64, Float64}())
    # TODO: Implement parallel tempering, see https://christian-igel.github.io/paper/TRBMAI.pdf
    set_visible_values!(rbm, fixed)
    e = energy(rbm)
    old_e = e + 1.0
    T = 1.0
    while (abs(e - old_e) > eps())
        # TODO: make this prettier...
        idx = rand(1:length(rbm.h))
        rbm.h[idx] = phx(rbm.W, rbm.x, rbm.b, T)[idx]
        idx = rand(1:length(rbm.x))
        rbm.x[idx] = pxh(rbm.W, rbm.h, rbm.a, T)[idx]
        set_visible_values!(rbm, fixed)
        old_e = e
        e = energy(rbm)
    end 
end

"Find valuations for the lbm that are consistent with its knowledge."
function reason!(
    rbm::LBM
    ;
    samples::Int64=100,
    factor::Int64=1,
    ϵ::Float64=.5,
    query::Dict{Int64, Float64}=Dict{Int64, Float64}(),
    verbose::Bool=false,
)
    sats = Set{Vector}()
    for iter in range(1, samples)
        rand_visible!(rbm)
        rand_hidden!(rbm)
        if verbose
            println("iteration $iter")
        end
        gibbs!(rbm, query)
        # from Proof of Theorem 1
        # s(x) = nergy(rbm) <= -f * ϵ, where f is the number of "rules"
        satisfaction = -(1/(factor * ϵ)) * energy(rbm)
        if satisfaction >= 1 
            push!(sats, rbm.x)
        end
    end
    return collect(sats)
end

"Convert a boolean to a float value in {-1.0, 1.0}."
b2d(b::Bool) = if b 1.0 else -1.0 end

"""
Construct weights and biases of a Restricted Boltzmann Machine from the
provided DNF formula, which has to be a full DNF.
"""
function sdnf2lbm(f::DNFFormula; ϵ::Float64=.5)
    nliterals = length(f.literals)
    nclauses = length(f.clauses)
    W::Matrix{Float64} = zeros(nclauses, nliterals)
    b::Vector{Float64} = zeros(nclauses)
    a::Vector{Float64} = zeros(nliterals)
    lit2idx = Dict(lit => idx for (idx, lit) in enumerate(f.literals))
    for (j, clause) in enumerate(f.clauses)
        STj = 0.0
        for literal in clause.literals
            W[j, lit2idx[literal.literal]] = b2d(!literal.neg)
            if !literal.neg
                STj -= 1.0
            end
        end
        b[j] = STj + ϵ
    end
    return LBM(W, a, b), lit2idx, Dict(idx => lit for (lit, idx) in lit2idx)
end

"""
Turns an LBM into a dot language graph.
"""
function lbm2dot(lbm::LBM, idx2lit::Dict{Int64, String})
    dot = "graph LBM {\n"
    dot *= "splines=false\n"
    dot *= "node [shape=rect]\n"
    dot *= "{\nrank=source\n"
    dot *= "node [shape=oval]\n"
    for (idx, b) in enumerate(lbm.b)
        dot *= "b$idx [label=\"$b\"]\n"
    end
    dot *= "}\n"
    dot *= "{\nrank=2\n"
    for (idx, h) in enumerate(lbm.h)
        dot *= "h$idx [label=\"h$idx\"]\n"
    end
    dot *= "}\n"
    dot *= "{\nrank=sink\n"
    for (idx, x) in enumerate(lbm.x)
        dot *= "x$idx [label=\"$(idx2lit[idx])\"]\n"
    end
    dot *= "}\n"
    for (j, h) in enumerate(lbm.h)
        dot *= "b$j -- h$j\n"
        for (i, x) in enumerate(lbm.x)
            dot *= "x$i -- h$j [label=\"$(lbm.W[j,i])\"]\n"
        end
    end
    dot *= "\n}"
end

"""
Plots the energy function. Should only be used with small RBM.
"""
function plot_energy_space!(lbm::LBM)
    energies = []
    confs = []
    csv_confs = []
    hconfs = collect(Iterators.product([0.0:1.0 for _ in lbm.h]...))
    xconfs = collect(Iterators.product([0.0:1.0 for _ in lbm.x]...))
    tpl_to_str(tpl) = begin
        str = ""
        for e in tpl
            str *= "$(Int(e))"
        end
        return str
    end
    for (idx, hconf) in enumerate(hconfs)
        for xconf in xconfs
            lbm.x = collect(xconf)
            lbm.h = collect(hconf)
            push!(energies, energy(lbm))
            push!(confs, "h $hconf x $xconf")
            push!(csv_confs, "\"$(tpl_to_str((hconf...,xconf...)))\"")
        end
    end
    open("plot_data.csv", "w+") do io
        write(io, "x,y\n")
        for row in zip(csv_confs, energies)
            write(io, "$(row[1]), $(row[2])\n")
        end
    end
    plotly()
    p = plot(confs, energies; minorticks=true, xrotation=60, xticks=:all, size=(800, 600))
    display(p)
end

"""
Plots the free energy function. Should only be used with small RBM.
"""
function plot_fenergy_space!(lbm::LBM, c::Float64, ϵ::Float64)
    energies = Vector{Float64}([])
    confs = Vector([])
    xconfs = collect(Iterators.product([0.0:1.0 for _ in lbm.x]...))
    for xconf in xconfs
        lbm.x = collect(xconf)
        push!(energies, fenergy(lbm, c))
        push!(confs, "$xconf")
        end
    plotly()
    p = plot(confs, energies; minorticks=true, xrotation=60, xticks=:all, size=(800, 600))
    plot!(p, confs, fill(-log(1 + exp(c * ϵ)), size(confs, 1)))
    display(p)
end
