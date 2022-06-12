"Logical Boltzmann Machines as in Tran 2021."
module LogicalBoltzmannMachines

#######
# WFF #
#######

include("formulae.jl")

# data structures

export
    AbstractFormula,
    EmptyFormula,
    SimpleFormula,
    And, Or,
    WFFormula,
    Valuation,
    DNF,
    ConjunctiveClause,
    DNFFormula

# functions

export 
    eval,
    is_full,
    dnf2wff,
    wff2dnf

#######
# lbm #
#######

include("machines.jl")

# data structures

export
    RBM

# functions

export
    energy,
    gibbs!,
    reason,
    sdnf2lbm

end