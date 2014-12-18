module io

const escaperules = Vector{ASCIIString}[
    ["\"", "&quot;"],
    ["\'", "&apos;"],
    ["<",  "&lt;"],
    [">",  "&gt;"],
    ["&",  "&amp;"]
    ]

function escape_xml(s::AbstractString)
    for (character, escapesequence) in reverse(escaperules)
        s = replace(s, character, escapesequence)
    end
    return s
end

function unescape_xml(s::AbstractString)
    for (character, escapesequence) in escaperules
        s = replace(s,escapesequence, character)
    end
    return s
end

function git_commitid(path)
    return readall(`git --work-tree $path rev-parse HEAD`)[1:end-1]
end

function git_commitdate(path)
    return readall(`git --work-tree $path show -s --format=%ci`)[1:end-1]
end

function write_gitinfo(f, tag, gitpath)
    ID = escape_xml(git_commitid(gitpath))
    date = escape_xml(git_commitdate(gitpath))
    writeln(f, "<$(tag)>")
    writeln(f, "<commit-id>$(ID)</commit-id>")
    writeln(f, "<commit-date>$(date)</commit-date>")
    writeln(f, "</$(tag)>")
end

function write_libraryinfo_julia(f)
    writeln(f, "<julia>\n")
    writeln(f, escape_xml(versioninfo()))
    writeln(f, "</julia>")
end

function write_libraryinfo_quantumoptics(f)
    tag = "quantumoptics-library"
    searchname = "julia-quantumoptics"
    for path in LOAD_PATH
        if contains(path, searchname)
            write_gitinfo(f, tag, path)
        end
    end
end

function write_libraryinfo_collectivespins(f)
    tag = "collectivespins-library"
    searchname = "collectivespins-library"
    for path in LOAD_PATH
        if contains(path, searchname)
            write_gitinfo(f, tag, path)
        end
    end
end

function write_libraryinfo(f)
    writeln(f, "<libraries>")
    write_libraryinfo_julia(f)
    write_libraryinfo_quantumoptics(f)
    write_libraryinfo_collectivespins(f)
    writeln(f, "</libraries>")
end

function write_systemdescription(f, system)
    writeln(f, "<systemdescription>")
    writeln(f, escape_xml(string(system)))
    writeln(f, "</systemdescription>")
end

function write_state{T}(f, state::Vector{T}, parameter, statetype)
    writetime = escape_xml(now())
    statetype = escape_xml(statetype)
    datatype = escape_xml(string(T))
    parameter = escape_xml(parameter)
    writeln(f, "<state writetime=\"$(writetime)\" statetype=\"$(statetype)\" datatype=\"$(datatype)\", parameter=\"$(parameter)\"")
    for x in state
        writeln(f, string(x))
        writeln(f, ',')
    end
    writeln(f, "</state>")
end

function write_state_meanfield{T}(f, state::Vector{T}, parameter)
    write_state(f, state, parameter, "meanfield")
end

function write_state_mpc{T}(f, state::Vector{T}, parameter)
    write_state(f, state, parameter, "meanfield")
end

function write_state_densityoperator{T}(f, state::Matrix{T}, parameter)
    write_state(f, state, parameter, "densityoperator")
end

end # module