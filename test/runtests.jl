names = readdir(".")

for name=names
    if startswith(name, "test_") && endswith(name, ".jl")
        println("Run $name")
        include(name)
    end
end
