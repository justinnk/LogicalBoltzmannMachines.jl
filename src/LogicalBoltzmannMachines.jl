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
    Lit, nLit,
    wff_eval,
    isfull,
    dnf2wff,
    wff2dnf

#######
# lbm #
#######

include("machines.jl")

# data structures

export
    LBM

# functions

export
    energy,
    fenergy,
    gibbs!,
    reason!,
    sdnf2lbm,
    lbm2dot,
    plot_energy_space!,
    plot_fenergy_space!

end