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
    write(f, "<$(tag)>\n")
    write(f, "<commit-id>$(ID)</commit-id>\n")
    write(f, "<commit-date>$(date)</commit-date>\n")
    write(f, "</$(tag)>\n")
end

function write_libraryinfo_julia(f)
    write(f, "<julia>\n")
    write(f, escape_xml(string(VERSION)))
    write(f, "\n")
    write(f, "</julia>\n")
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
    write(f, "<libraries>\n")
    write_libraryinfo_julia(f)
    write_libraryinfo_quantumoptics(f)
    write_libraryinfo_collectivespins(f)
    write(f, "</libraries>\n")
end

function write_systemdescription(f, system)
    write(f, "<systemdescription>\n")
    write(f, escape_xml(string(system)))
    write(f, "\n")
    write(f, "</systemdescription>\n")
end

function write_simulationparameters(f, parameters)
    write(f, "<simulationparameters>\n")
    write(f, escape_xml(string(parameters)))
    write(f, "\n")
    write(f, "</simulationparameters>\n")
end

function write_starttime(f)
    starttime = escape_xml(string(now()))
    write(f, "<starttime>\n")
    write(f, starttime)
    write(f, "\n")
    write(f, "</starttime>\n")
end

function write_head(f, system, parameters)
    write(f, "<head>\n")
    write_starttime(f)
    write_libraryinfo(f)
    write_systemdescription(f, system)
    write_simulationparameters(f, parameters)
    write(f, "</head>\n")
end

function write_state{T}(f, statetype, state::Vector{T}; parameter=nothing, time=true)
    writetime = time ? " writetime = \"$(escape_xml(string(now())))\"" : ""
    statetype = escape_xml(statetype)
    datatype = escape_xml(string(T))
    parameter = parameter!=nothing ? " parameter = \"$(escape_xml(parameter))\"" : ""
    write(f, "<state statetype=\"$(statetype)\" datatype=\"$(datatype)\"$writetime$parameter>")
    for (i, x) in enumerate(state)
        if i!=1
            write(f, ',')
        end
        write(f, string(x))
    end
    write(f, "</state>\n")
end

function write_state_meanfield{T}(f, state::Vector{T}, parameter)
    write_state(f, state, parameter, "meanfield")
end

function write_state_mpc{T}(f, state::Vector{T}, parameter)
    write_state(f, state, parameter, "mpc")
end

function write_state_densityoperator{T}(f, state::Matrix{T}, parameter)
    write_state(f, vec(state.data), parameter, "densityoperator")
end

function save_timeevolution(path, system, time, states, statetype, parameters)
    @assert length(time)==length(states)
    f = open(path, "w")
    write_head(f, system, parameters)
    for i=1:length(time)
        write_state(f, vec(states[i].data), "t=$(time[i])", statetype)
    end
    close(f)
end

end # module