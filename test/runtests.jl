names = [
    "test_geometry.jl",
    "test_effective_interaction.jl",
    "test_system.jl",

    "test_meanfield_equations.jl",
    "test_mpc_equations.jl",
    "test_symmetric_equations.jl",

    "test_mpc.jl",
    "test_rotate.jl",
    "test_squeezing.jl",
    "test_interaction.jl",
    "test_reducedspin.jl"
]

detected_tests = filter(
    name->startswith(name, "test_") && endswith(name, ".jl"),
    readdir("."))

unused_tests = setdiff(detected_tests, names)
if length(unused_tests) != 0
    error("The following tests are not used:\n", join(unused_tests, "\n"))
end

unavailable_tests = setdiff(names, detected_tests)
if length(unavailable_tests) != 0
    error("The following tests could not be found:\n", join(unavailable_tests, "\n"))
end

for name=names
    if startswith(name, "test_") && endswith(name, ".jl")
        include(name)
    end
end
