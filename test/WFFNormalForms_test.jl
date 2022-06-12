using Test, LogicalBoltzmannMachines

f_full = DNFFormula([
    ConjunctiveClause([
        SimpleFormula('x', true),
        SimpleFormula('y', true),
        SimpleFormula('z', true)
    ]),
    ConjunctiveClause([
        SimpleFormula('x', true),
        SimpleFormula('y'),
        SimpleFormula('z')
    ])
])


f_impl = DNFFormula([
    ConjunctiveClause([SimpleFormula('x')]),
    ConjunctiveClause([SimpleFormula('y', true)])
])

@testset "DNF Module Tests" begin
    @testset "isfull Tests" begin
        @test isfull(f_full)
        @test !isfull(f_impl)
    end
end