using Printf: @printf

##############################
# Well-Formed Formulae (WFF) #
##############################

"An abstract propositional logic formula."
abstract type AbstractFormula end

"The empty formula."
struct EmptyFormula <: AbstractFormula
end

"An atomic proposition."
struct SimpleFormula <: AbstractFormula
    "The one-letter name of the literal"
    literal::Char
    "Whether the literal is in negated form"
    neg::Bool
    SimpleFormula(literal, neg) = new(literal, neg)
    SimpleFormula(literal) = new(literal, false)
end

# convenience constructors

"Creates a positive literal with name x."
Lit(x) = SimpleFormula(x)
"Creates a negative literal with name x."
nLit(x) = SimpleFormula(x, true)

"The conjunction of two propositions."
struct And <: AbstractFormula
    left::AbstractFormula
    right::AbstractFormula
end

"The disjunction of two propositions."
struct Or <: AbstractFormula
    left::AbstractFormula
    right::AbstractFormula
end

"A well-formed propositional logic formula."
struct WFFormula
    formula::AbstractFormula
    WFFormula(formula) = new(formula)
end

"An assignment of truth values for each name in a well-formed formula."
const Valuation = Dict{Char,Bool}

# evaluation

wff_eval(formula::SimpleFormula, v::Valuation) =
    v[formula.literal] ⊻ formula.neg
wff_eval(formula::And, v::Valuation) =
    wff_eval(formula.left, v) && wff_eval(formula.right, v)
wff_eval(formula::Or, v::Valuation) =
    wff_eval(formula.left, v) || wff_eval(formula.right, v)
wff_eval(formula::AbstractFormula, v::Valuation) = wff_eval(formula, v)
"Evaluate the truth value of the provided formula under the given valuation."
wff_eval(formula::WFFormula, v::Valuation) = wff_eval(formula.formula, v)

# pretty printing
        
function Base.show(io::IO, f::And)
    @printf(io, " ( ")
    show(io, f.left)
    @printf(io, " ∧ ")
    show(io, f.right)
    @printf(io, " ) ")
end
function Base.show(io::IO, f::Or)
    @printf(io, " ( ")
    show(io, f.left)
    @printf(io, " ∨ ")
    show(io, f.right)
    @printf(io, " ) ")
end
function Base.show(io::IO, f::SimpleFormula)
    neg = if f.neg "¬" else "" end 
    @printf(io, "%s%c", neg, f.literal)
end
Base.show(io::IO, f::EmptyFormula) = return
Base.show(io::IO, f::WFFormula) = show(io, f.formula)

###########################
# Disjunctive Normal Form #
###########################

"The Disjunctive Normal Form of a well-formed formula."
abstract type DNF <: AbstractFormula end

"A list of conjoined literals."
struct ConjunctiveClause <: DNF
    "The literals conjoined in this clause"
    literals::Array{SimpleFormula}
    ConjunctiveClause(literals::Array{SimpleFormula}) = new(literals)
end

"A disjunctive normal form of a well-formed formula."
struct DNFFormula <: AbstractFormula
    "The clauses disjoined in this dnf"
    clauses::Array{ConjunctiveClause}
    "The set of literals in this formula"
    literals::Set{Char}
    DNFFormula(clauses, literals) = new(clauses, literals)
end

function DNFFormula(clauses)
    literals::Set{Char} = Set{Char}()
    for c in clauses
        push!(literals, [l.literal for l in c.literals]...)
    end
    DNFFormula(clauses, literals)
end

"Returns whether the provided DNF is full (each clause contains all literals)"
function isfull(f::DNFFormula)
    for c in f.clauses
        for l in f.literals
            if l ∉ map(x -> x.literal, c.literals)
                return false
            end
        end
    end
    return true
end

"Converts the data structure for a DNF to the underlying well-formed formula."
function dnf2wff(f::DNFFormula)::WFFormula
    wff = EmptyFormula()
    for c in f.clauses
        conjf = c.literals[1]
        for cc in c.literals[2:end]
            conjf = And(conjf, cc)
        end
        if typeof(wff) == EmptyFormula
            wff = conjf
        else
            wff = Or(wff, conjf)
        end
    end
    WFFormula(wff)
end

"Converts a well-formed formula to disjunctive normal form."
function wff2dnf(f::WFFormula)::DNFFormula
    return
end

# pretty printing

function Base.show(io::IO, dnf::DNFFormula)
    @printf(io, "SKj: ")
    show(io, dnf.SK)
    @printf(io, "\nSTj: ")
    show(io, dnf.ST)
    @printf(io, "\n")
    show(io, dnf2wff(dnf))
end
