
"A Restricted Boltzmann Machine."
struct RBM
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
    RBM(
        x::Vector{Float64},
        a::Vector{Float64},
        h::Vector{Float64},
        b::Vector{Float64},
        W::Matrix{Float64}
    ) = new(x, a, h, b, W)
end

"Creates a Restricted Boltzmann Machine with random visible and hidden values."
function RBM(a::Vector{Float64}, b::Vector{Float64}, W::Matrix{Float64})
    size(a, 1) == size(W, 2) || error(
        "Bias vector for the visible layer is incompatible with weight Matrix."
    )
    size(b, 1) == size(W, 1) || error(
        "Bias vector for the hidden layer is incompatible with weight Matrix."
    )
    RBM(rand([0.0, 1.0], size(a, 1)), a, rand([0.0, 1.0], size(b, 1)), b, W)
end

"Creates a new RBM with the provided visible layer values."
function update_visible(rbm::RBM, x::Vector{Float64})
    RBM(x, rbm.a, rbm.h, rbm.b, rbm.W)
end

"Creates a new RBM with the provided hidden layer values."
function update_hidden(rbm::RBM, h::Vector{Float64})
    RBM(rbm.x, rbm.a, h, rbm.b, rbm.W)
end

"Calculate the energy of the provided RBM."
energy(rbm::RBM) = -(rbm.W * rbm.x)' * rbm.h - rbm.a' * rbm.x - rbm.b' * rbm.h

"The sigmoid of x."
sigmoid(x::Float64, τ::Float64) = (1.0 / (1.0 + exp(-(1/τ) * x)))

"The probability of hidden layer values given visible layer values."
phx(W::Matrix{Float64}, x::Vector{Float64}, b::Vector{Float64}, τ::Float64) =
    round.(sigmoid.(W * x + b, τ))

"The probability of visible layer values given hidden layer values."
pxh(W::Matrix{Float64}, h::Vector{Float64}, a::Vector{Float64}, τ::Float64) =
    round.(sigmoid.(W' * h + a, τ))

"Set visible layer values for the given RBM at certain indices."
function set_visible_values!(rbm::RBM, values::Dict{Int64, Float64})
    for (idx, f) in values
        rbm.x[idx] = f
    end
    rbm
end

"Perform Gibbs sampling on the given RBM, while fixating some visible values."
function gibbs!(rbm::RBM, fixed::Dict{Int64, Float64}=Dict{Int64, Float64}())
    rbm = set_visible_values!(rbm, fixed)
    e = energy(rbm)
    old_e = 0.0
    τ = 1.0
    while (abs(e - old_e) > eps())
        rbm = update_hidden(rbm, phx(rbm.W, rbm.x, rbm.b, τ))
        rbm = update_visible(rbm, pxh(rbm.W, rbm.h, rbm.a, τ))
        rbm = set_visible_values!(rbm, fixed)
        old_e = e
    end
    return rbm
end

"Find valuations for the lbm that are consistent with its knowledge."
function reason(
    W::Matrix{Float64},
    a::Vector{Float64},
    b::Vector{Float64},
    ϵ::Float64=.5,
    query::Dict{Int64, Float64}=Dict{Int64, Float64}(),
    samples::Int64=100
)
    sats = Set{Vector}()
    for iter in range(1, samples)
        rbm = RBM(a, b, W)
        if iter % 10 == 0
            println("iteration $iter")
        end
        rbm = gibbs!(rbm, query)
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
    return W, a, b
end