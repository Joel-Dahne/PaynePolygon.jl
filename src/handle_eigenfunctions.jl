"""
    save_eigenfunction(filename, u, λ; symmetry_class)

Stores the data for the eigenfunction and eigenvalue in a jlf-file
so that it can be loaded in later with [`load_eigenfunctions`](@ref).

It only handles the eigenfunctions defined by
`MethodOfParticularSolutions.example_domain_goal_v1`.

It stores the coefficients for the eigenfunction, the precision it
uses, the symmetry class and the eigenvalue. The values stored are
also returned.
"""
function save_eigenfunction(
    filename::AbstractString,
    u::AbstractEigenfunction,
    λ;
    symmetry_class::Integer,
)
    # Dump coefficients
    cs_dump = ArbTools.arb_dump.(u.domain.parent.(coefficients(u)))

    # Dump λ
    λ_dump = ArbTools.arb_dump(u.domain.parent(λ))

    save(
        filename,
        "cs_dump",
        cs_dump,
        "precision",
        precision(u.domain.parent),
        "symmetry_class",
        symmetry_class,
        "λ_dump",
        λ_dump,
    )

    return cs_dump, precision, symmetry_class, λ_dump
end

"""
    load_eigenfunction(filename, N = 27, d = 11, h = 6)

Load an eigenfunction which have been stored with
[`save_eigenfunctions`](@ref), see its documentation for more details.

It returns the domain, the eigenfunction and the eigenvalue.
"""
function load_eigenfunction(
    filename::AbstractString,
    N::Integer = 27,
    d::Integer = 11,
    h::Integer = 6;
    T = arb,
)
    values = load(filename)

    parent = RealField(values["precision"])

    domain, u = MethodOfParticularSolutions.example_domain_goal_v1(
        N,
        d,
        h,
        parent,
        symmetry_class = values["symmetry_class"];
        T,
    )

    set_eigenfunction!(u, ArbTools.arb_load_dump.(values["cs_dump"], Ref(parent)))

    λ = ArbTools.arb_load_dump(values["λ_dump"], parent)

    return domain, u, λ
end
