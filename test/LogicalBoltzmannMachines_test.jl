using Test, LogicalBoltzmannMachines

# expected values
a = [0., 0., 0.]
b = [0.5, -1.5, -1.5, -1.5]
W = [-1. -1. -1.;
     1. 1. -1.;
     1. -1. 1.;
     -1. 1. 1.]

# corresponding formula
f_xor = DNFFormula([
    ConjunctiveClause([
        SimpleFormula('x', true),
        SimpleFormula('y', true),
        SimpleFormula('z', true)
    ]),
    ConjunctiveClause([
        SimpleFormula('x'),
        SimpleFormula('y'),
        SimpleFormula('z', true)
    ]),
    ConjunctiveClause([
        SimpleFormula('x'),
        SimpleFormula('y', true),
        SimpleFormula('z')
    ]),
    ConjunctiveClause([
        SimpleFormula('x', true),
        SimpleFormula('y'),
        SimpleFormula('z')
    ]),
])

@testset "LogicalBoltzmannMachines module tests" begin
    @testset "dnf2lbm tests" begin
        @test dnf2lbm(f_xor) = (W, a, b)
    end
end