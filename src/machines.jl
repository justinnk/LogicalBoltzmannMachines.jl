
"A Restricted Boltzmann Machine."
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

"Initialize the visible layer to a random vector."
rand_visible!(rbm::LBM) = rbm.x = rand([0.0, 1.0], size(rbm.a, 1))
"Initialize the visible layer to a random vector."
rand_hidden!(rbm::LBM) = rbm.h = rand([0.0, 1.0], size(rbm.b, 1))

"Calculate the energy of the provided RBM."
energy(rbm::LBM) = -(rbm.W * rbm.x)' * rbm.h - rbm.a' * rbm.x - rbm.b' * rbm.h

"The sigmoid of x."
sigmoid(x::Float64, T::Float64) = (1.0 / (1.0 + exp(-(1/T) * x)))

"The probability of hidden layer values given visible layer values P(h|x)."
phx(W::Matrix{Float64}, x::Vector{Float64}, b::Vector{Float64}, T::Float64) =
    round.(sigmoid.(W * x + b, T))

"The probability of visible layer values given hidden layer values P(x|h)."
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
    set_visible_values!(rbm, fixed)
    e = energy(rbm)
    old_e = 0.0
    T = 1.0
    while (abs(e - old_e) > eps())
        rbm.h = phx(rbm.W, rbm.x, rbm.b, T)
        rbm.x = pxh(rbm.W, rbm.h, rbm.a, T)
        set_visible_values!(rbm, fixed)
        old_e = e
    end
end

"Find valuations for the lbm that are consistent with its knowledge."
function reason(
    rbm::LBM
    ;
    samples::Int64=100,
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
        satisfaction = -(1/ϵ) * energy(rbm)
        if satisfaction == 1.0
            push!(sats, rbm.x)
        end
    end
    return sats
end

"Convert a boolean to a float value in {-1.0, 1.0}."
b2d(b::Bool) = if b 1.0 else -1.0 end

"""
Construct weights and biases of a Restricted Boltzmann Machine from the
provided DNF formula, which has to be a full DNF.
"""
function sdnf2lbm(f::DNFFormula, ϵ::Float64=.5)
    isfull(f) || error("The provided formula has to be full.")
    nliterals = length(f.literals)
    nclauses = length(f.clauses)
    W::Matrix{Float64} = zeros(nclauses, nliterals)
    b::Vector{Float64} = zeros(nclauses)
    a::Vector{Float64} = zeros(nliterals)
    for (j, clause) in enumerate(f.clauses)
        STj = 0.0
        for (i, literal) in enumerate(clause.literals)
            W[j,i] = b2d(!literal.neg)
            if !literal.neg
                STj -= 1.0
            end
        end
        b[j] = STj + ϵ
    end
    return LBM(W, a, b)
end