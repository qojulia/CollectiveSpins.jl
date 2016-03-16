module meanfield

using ..interaction, ..system

using quantumoptics

#export symmetric_meanfield_timeevolution, meanfield_timeevolution, correlation_timeevolution


function symmetric_meanfield_timeevolution(T, omega_eff, gamma, gamma_eff)
    function f(t, y::Vector{Complex128}, dy::Vector{Complex128})
        sp = y[1]
        sm = y[2]
        sz = y[3]
        dsp = 1im*omega_eff*sp*sz - 0.5*gamma*sp - 0.5*gamma_eff*sp*sz
        dsm = -1im*omega_eff*sm*sz - 0.5*gamma*sm - 0.5*gamma_eff*sm*sz
        dsz = gamma*(1-sz) + 2*gamma_eff*sp*sm
        dy[1] = dsp
        dy[2] = dsm
        dy[3] = dsz
    end
    tout = Float64[]
    sp_out = Complex128[]
    sm_out = Complex128[]
    sz_out = Complex128[]
    function fout(t, y::Vector{Complex128})
        push!(tout, t)
        push!(sp_out, y[1])
        push!(sm_out, y[2])
        push!(sz_out, y[3])
    end

    y0 = Complex128[-0.5,-0.5,0.]
    #T = [0.:0.1:5.]
    quantumoptics.ode_dopri.ode(f, T, y0, fout=fout)
    return T, sp_out, sm_out, sz_out
end

function GammaMatrix(system)
    atoms = system.atoms
    N = length(atoms)
    Γ = zeros(Float64, N, N)
    for i=1:N, j=1:N
        polarization = atoms[i].polarization
        Γ[i,j] = Gamma(atoms[i].position, atoms[j].position, polarization, system.gamma)
    end
    return Γ
end

function OmegaMatrix(system)
    atoms = system.atoms
    N = length(atoms)
    Ω = zeros(Float64, N, N)
    for i=1:N, j=1:N
        if i==j
            continue
        end
        polarization = atoms[i].polarization
        Ω[i,j] = Omega(atoms[i].position, atoms[j].position, polarization, system.gamma)
    end
    return Ω
end

const I = identity(spinbasis)

function correlations(rho::Operator)
    N = length(rho.basis_l.bases)
    C = Dict{AbstractString, Matrix{Float64}}()
    m = zeros(Float64, N, N)
    C["xx"] = deepcopy(m); C["yx"] = deepcopy(m); C["zx"] = deepcopy(m);
    C["xy"] = deepcopy(m); C["yy"] = deepcopy(m); C["zy"] = deepcopy(m);
    C["xz"] = deepcopy(m); C["yz"] = deepcopy(m); C["zz"] = deepcopy(m);
    cor(s1, s2, k, l) = expect(quantumoptics.embed(rho.basis_l, [k,l], [s1,s2]), rho)
    expvalue(s,i) = expect(quantumoptics.embed(rho.basis_l, i, s), rho)
    cov(s1, s2, k, l) = cor(s1, s2, k, l) - expvalue(s1, k)*expvalue(s2, l)
    for k=1:N, l=1:N
        if k==l
            continue
        end
        C["xx"][k,l] = real(cov(sigmax, sigmax, k, l));
        C["xy"][k,l] = real(cov(sigmax, sigmay, k, l));
        C["xz"][k,l] = real(cov(sigmax, sigmaz, k, l));
        C["yx"][k,l] = real(cov(sigmay, sigmax, k, l));
        C["yy"][k,l] = real(cov(sigmay, sigmay, k, l));
        C["yz"][k,l] = real(cov(sigmay, sigmaz, k, l));
        C["zx"][k,l] = real(cov(sigmaz, sigmax, k, l));
        C["zy"][k,l] = real(cov(sigmaz, sigmay, k, l));
        C["zz"][k,l] = real(cov(sigmaz, sigmaz, k, l));
    end
    return C
end

function densityoperator{T<:Real}(σx::T, σy::T, σz::T)
    return 0.5*(I + σx*sigmax + σy*sigmay + σz*sigmaz)
end

function densityoperator{T<:Real}(σx::Vector{T}, σy::Vector{T}, σz::Vector{T})
    N = length(σx)
    result = densityoperator(σx[1], σy[1], σz[1])
    for k=2:N
        result = tensor(result, densityoperator(σx[k], σy[k], σz[k]))
    end
    return result
end

function densityoperator{T<:Real}(σx::Vector{T}, σy::Vector{T}, σz::Vector{T}, C::Dict{AbstractString, Matrix{T}})
    N = length(σx)
    ρ = densityoperator(σx, σy, σz)
    ρlist = DenseOperator[densityoperator(σx[k], σy[k], σz[k]) for k=1:N]
    function cor(s1,s2,k,l)
        if k==1
            result = s1
        elseif l==1
            result = s2
        else
            result = ρlist[1]
        end
        for i=2:N
            if k==i
                result = tensor(result, s1)
            elseif l==1i
                result = tensor(result, s2)
            else
                result = tensor(result, ρlist[i])
            end
        end
        return result
    end
    for k=1:N, l=k+1:N
        ρ += 0.25*(C["xx"][k,l]*cor(sigmax, sigmax, k, l) + C["yx"][k,l]*cor(sigmay, sigmax, k, l) + C["zx"][k,l]*cor(sigmaz, sigmax, k, l)
            + C["xy"][k,l]*cor(sigmax, sigmay, k, l) + C["yy"][k,l]*cor(sigmay, sigmay, k, l) + C["zy"][k,l]*cor(sigmaz, sigmay, k, l)
            + C["xz"][k,l]*cor(sigmax, sigmaz, k, l) + C["yz"][k,l]*cor(sigmay, sigmaz, k, l) + C["zz"][k,l]*cor(sigmaz, sigmaz, k, l))
    end
    return ρ
end

function meanfield_timeevolution(T, system, sz0, sp0)
    N = length(system.atoms)
    y0 = zeros(Complex128, 2*N)
    ω0 = Float64[atom.omega0 for atom in system.atoms]
    Ω = OmegaMatrix(system)
    Γ = GammaMatrix(system)
    function split_parameters(y)
        sz = y[1:N]
        sp = y[N+1:2*N]
        return sz, sp
    end
    function combine_parameters(dy, dsz, dsp)
        dy[1:N] = dsz
        dy[N+1:2*N] = dsp
    end
    combine_parameters(y0, sz0, sp0)

    function f(t, y::Vector{Complex128}, dy::Vector{Complex128})
        sz, sp = split_parameters(y)
        dsz, dsp = split_parameters(dy)
        sm = conj(sp)
        for k=1:N
            dsz[k] = -Γ[k,k]*(sz[k]-1)
            dsp[k] = -1im*ω0[k]*sp[k] - 0.5*Γ[k,k]*sp[k]
            for i=1:N
                if i==k
                    continue
                end
                dsz[k] += 2im*Ω[k,i]*(sp[k]*sm[i] - sp[i]*sm[k]) + Γ[k,i]*(sp[k]*sm[i] + sp[i]*sm[k])
                dsp[k] += 1im*Ω[k,i]*sz[k]*sp[i] -0.5*Γ[k,i]*sz[k]*sp[i]
            end
        end
        combine_parameters(dy, dsz, dsp)
    end

    tout = Float64[]
    sp_out = Vector{Complex128}[]
    sz_out = Vector{Complex128}[]
    function fout(t, y::Vector{Complex128})
        sz, sp = split_parameters(y)
        push!(tout, t)
        push!(sz_out, sz)
        push!(sp_out, sp)
    end

    quantumoptics.ode_dopri.ode(f, T, y0, fout=fout)
    return T, sz_out, sp_out
end

function correlation_timeevolution(T, system, sz0, sp0, Cpm0, Cpz0, Cpp0, Czz0)
    N = length(system.atoms)
    ω0 = Float64[atom.omega0 for atom in system.atoms]
    Ω = OmegaMatrix(system)
    Γ = GammaMatrix(system)
    function split_parameters(y)
        sz = y[0*N+1:1*N]
        sp = y[1*N+1:2*N]
        Cpm = reshape(y[2*N+0*N^2+1:2*N+1*N^2], N, N)
        Cpz = reshape(y[2*N+1*N^2+1:2*N+2*N^2], N, N)
        Cpp = reshape(y[2*N+2*N^2+1:2*N+3*N^2], N, N)
        Czz = reshape(y[2*N+3*N^2+1:2*N+4*N^2], N, N)
        return sz, sp, Cpm, Cpz, Cpp, Czz
    end
    function combine_parameters(y, sz, sp, Cpm, Cpz, Cpp, Czz)
        y[0*N+1:1*N] = sz
        y[1*N+1:2*N] = sp
        y[2*N+0*N^2+1:2*N+1*N^2] = reshape(Cpm, N^2)
        y[2*N+1*N^2+1:2*N+2*N^2] = reshape(Cpz, N^2)
        y[2*N+2*N^2+1:2*N+3*N^2] = reshape(Cpp, N^2)
        y[2*N+3*N^2+1:2*N+4*N^2] = reshape(Czz, N^2)
    end
    y0 = zeros(Complex128, 2*N + 4*N^2)
    combine_parameters(y0, sz0, sp0, Cpm0, Cpz0, Cpp0, Czz0)

    function f(t, y::Vector{Complex128}, dy::Vector{Complex128})
        sz, sp, Cpm, Cpz, Cpp, Czz = split_parameters(y)
        dsz, dsp, dCpm, dCpz, dCpp, dCzz = split_parameters(zeros(Complex128, 2*N + 4*N^2))
        Cmp = transpose(Cpm)
        Cmm = conj(Cpp)
        Cmz = conj(Cpz)
        Czp = transpose(Cpz)
        Czm = ctranspose(Cpz)
        sm = conj(sp)
        for k=1:N
            dsz[k] = Γ[k,k]*(1-sz[k])
            dsp[k] = -1im*ω0[k]*sp[k] - 0.5*Γ[k,k]*sp[k]
            for i=1:N
                if i==k
                    continue
                end
                dsz[k] += (1im*2*Ω[k,i]*(sp[k]*sm[i] + Cpm[k,i] - sp[i]*sm[k] - Cpm[i,k])
                            + Γ[k,i]*(sp[k]*sm[i] + Cpm[k,i] + sp[i]*sm[k] + Cpm[i,k])
                        )
                dsp[k] += (1im*Ω[k,i]*(sp[i]*sz[k] + Cpz[i,k])
                            - 0.5*Γ[k,i]*(sp[i]*sz[k] + Cpz[i,k])
                        )
            end
        end
        for k=1:N, l=1:N
            if k==l
                continue
            end
            #dCpm[k,l] = 1im*Ω[k,l]*(0.5*sz[k] - 0.5*sz[l] + sp[k]*(sm[k]*sz[l]+Cmz[k,l]) - sm[l]*(sz[k]*sp[l]+Czp[k,l]))-0.5*(Γ[k,k]+Γ[l,l])*Cpm[k,l] - 0.5*Γ[k,l]*(0.5*sz[k]+0.5*sz[l]-sz[k]*sz[l]-Czz[k,l])+0.5*Γ[k,l]*sp[k]*(sm[k]*sz[l]+Cmz[k,l]) + 0.5*Γ[k,l]*sm[l]*(sz[k]*sp[l]+Czp[k,l])
            dCpm[k,l] = (-1im*(ω0[k]-ω0[l])*Cpm[k,l]
                         -1im*Ω[k,l]*(0.5*sz[l] - 0.5*sz[k] - sp[k]*(sm[k]*sz[l]+Cmz[k,l]) + sm[l]*(sz[k]*sp[l]+Czp[k,l]))
                         -0.5*(Γ[k,k]+Γ[l,l])*Cpm[k,l]
                         -0.5*Γ[k,l]*(0.5*sz[l] + 0.5*sz[k] - (sz[k]*sz[l]+Czz[k,l]) - sp[k]*(sm[k]*sz[l]+Cmz[k,l]) - sm[l]*(sz[k]*sp[l]+Czp[k,l]))
                        )

            #dCpz[k,l] = -1im*ω0[k]*Cpz[k,l] - 1im*Ω[k,l]*(-sp[l] - 2*sp[k]*sp[k]*sm[l] + 2*sp[k]*sm[k]*sp[l] + sz[k]*sp[l]*sz[l]-2*sp[k]*Cpm[k,l] + 2*sp[k]*Cmp[k,l] + sz[l]*Czp[k,l])-0.5*(Γ[k,k]+2*Γ[l,l])*Cpz[k,l] - 0.5*Γ[k,l]*(-sp[l] + 2*sz[k]*sp[l] + 2*sp[k]*sm[k]*sp[l] + 2*sp[k]*sp[k]*sm[l] - sz[k]*sp[l]*sz[l]+ 2*Czp[k,l] + 2*sp[k]*Cmp[k,l] - sz[l]*Czp[k,l] + 2*sp[k]*Cpm[k,l])
            dCpz[k,l] = (-1im*ω0[k]*Cpz[k,l]
                         +1im*Ω[k,l]*(sp[l] + 2*sp[k]*(sp[k]*sm[l]+Cpm[k,l] - sm[k]*sp[l]-Cmp[k,l]) - sz[l]*(sz[k]*sp[l]+Czp[k,l]))
                         -0.5*(Γ[k,k]+2*Γ[l,l])*Cpz[k,l]
                         +0.5*Γ[k,l]*(sp[l] - 2*sp[k]*(sm[k]*sp[l]+Cmp[k,l] + sp[k]*sm[l]+Cpm[k,l])
                                            - 2*(sz[k]*sp[l]+Czp[k,l]) + sz[l]*(sz[k]*sp[l]+Czp[k,l]))
                        )

            #dCpp[k,l] = -1im*(ω0[k]+ω0[l])*Cpp[k,l] + (0.5*Γ[k,l]-1im*Ω[k,l])*(sp[k]*(sp[k]*sz[l]+Cpz[k,l]) + sp[l]*(sp[l]*sz[k]+Czp[k,l]))-0.5*(Γ[k,k]+Γ[l,l])*Cpp[k,l]
            dCpp[k,l] = (-1im*(ω0[k]+ω0[l])*Cpp[k,l]
                         -0.5*(Γ[k,k]+Γ[l,l])*Cpp[k,l]
                         +(0.5*Γ[k,l]-1im*Ω[k,l])*(sp[k]*(sp[k]*sz[l]+Cpz[k,l]) + sp[l]*(sz[k]*sp[l]+Czp[k,l]))
                        )

            #dCzz[k,l] = -2im*Ω[k,l]*sz[k]*(sm[k]*sp[l] - sp[k]*sm[l] + Cmp[k,l] - Cpm[k,l])-2im*Ω[k,l]*sz[l]*(sp[k]*sm[l] - sm[k]*sp[l] + Cpm[k,l] - Cmp[k,l]) -(Γ[k,k]+Γ[l,l])*Czz[k,l] + Γ[k,l]*(2*sm[k]*sp[l] + 2*sp[k]*sm[l] - sz[k]*sm[k]*sp[l] - sm[k]*sp[l]*sz[l]- sp[k]*sm[l]*sz[l] - sz[k]*sp[k]*sm[l]+ 2*Cmp[k,l] - sz[k]*Cmp[k,l] - sz[l]*Cmp[k,l]+ 2*Cpm[k,l] - sz[l]*Cpm[k,l] - sz[k]*Cpm[k,l])
            dCzz[k,l] = (-1im*2*Ω[k,l]*(sz[k]-sz[l])*(sm[k]*sp[l]+Cmp[k,l] - sp[k]*sm[l]-Cpm[k,l])
                         -(Γ[k,k]+Γ[l,l])*Czz[k,l]
                         +Γ[k,l]*(1-sz[k]+1-sz[l])*(sm[k]*sp[l]+Cmp[k,l] + sp[k]*sm[l]+Cpm[k,l])
                        )
            for i=1:N
                if i==l || i==k
                    continue
                end
                #dCpm[k,l] += -1im*Ω[l,i]*(sz[l]*Cpm[k,i]+sm[i]*Cpz[k,l]) +1im*Ω[k,i]*(sp[i]*Czm[k,l]+sz[k]*Cpm[i,l]) -0.5*Γ[l,i]*(sz[l]*Cpm[k,i] + sm[i]*Cpz[k,l]) -0.5*Γ[k,i]*(sp[i]*Czm[k,l]+sz[k]*Cpm[i,l])
                dCpm[k,l] += (1im*Ω[k,i]*(sz[k]*Cpm[i,l] + sp[i]*Czm[k,l])
                              -1im*Ω[l,i]*(sz[l]*Cpm[k,i] + sm[i]*Cpz[k,l])
                              -0.5*Γ[k,i]*(sz[k]*Cpm[i,l] + sp[i]*Czm[k,l])
                              -0.5*Γ[l,i]*(sz[l]*Cmp[i,k] + sm[i]*Cpz[k,l])
                              )

                #dCpz[k,l] += -1im*2*Ω[l,i]*(-sp[l]*Cpm[k,i]-sm[i]*Cpp[k,l]+sp[i]*Cpm[k,l]+sm[l]*Cpp[k,i]) + 1im*Ω[k,i]*(sp[i]*Czz[k,l]+sz[k]*Cpz[i,l]) + Γ[l,i]*(sp[i]*Cpm[k,l]+sm[l]*Cpp[i,k]+sp[l]*Cpm[k,i]+sm[i]*Cpp[k,l]) - 0.5*Γ[k,i]*(sp[i]*Czz[k,l]+sz[k]*Cpz[i,l])
                dCpz[k,l] += (1im*Ω[k,i]*(sz[k]*Cpz[i,l] + sp[i]*Czz[k,l])
                              +1im*2*Ω[l,i]*(sp[l]*Cmp[i,k] + sm[i]*Cpp[k,l] - sm[l]*Cpp[i,k] - sp[i]*Cpm[k,l])
                              -0.5*Γ[k,i]*(sz[k]*Cpz[i,l] + sp[i]*Czz[k,l])
                              +Γ[l,i]*(sp[l]*Cmp[i,k] + sm[i]*Cpp[k,l] + sm[l]*Cpp[i,k] + sp[i]*Cpm[k,l])
                              )

                #dCpp[k,l] += 1im*Ω[l,i]*(sp[i]*Cpz[k,l]+sz[l]*Cpp[k,i]) + 1im*Ω[k,i]*(sp[i]*Czp[k,l] + sz[k]*Cpp[i,l]) -0.5*Γ[l,i]*(sp[i]*Cpz[k,l]+sz[l]*Cpp[i,k]) - 0.5*Γ[k,i]*(sp[i]*Czp[k,l] + sz[k]*Cpp[i,l])
                dCpp[k,l] += (1im*Ω[k,i]*(sz[k]*Cpp[i,l] + sp[i]*Czp[k,l])
                              +1im*Ω[l,i]*(sz[l]*Cpp[i,k] + sp[i]*Cpz[k,l])
                              -0.5*Γ[k,i]*(sz[k]*Cpp[i,l] + sp[i]*Czp[k,l])
                              -0.5*Γ[l,i]*(sz[l]*Cpp[i,k] + sp[i]*Cpz[k,l])
                              )

                #dCzz[k,l] += -1im*2*Ω[l,i]*(-sp[l]*Czm[k,i]-sm[i]*Czp[k,l]+sp[i]*Czm[k,l]+sm[l]*Czp[k,i])-1im*2*Ω[k,i]*(sp[i]*Cmz[k,l]+sm[k]*Cpz[i,l]-sp[k]*Cmz[i,l]-sm[i]*Cpz[k,l]) + Γ[l,i]*(sp[i]*Czm[k,l]+sm[l]*Cpz[i,k]+sp[l]*Czm[k,i]+sm[i]*Czp[k,l])+Γ[k,i]*(sp[i]*Cmz[k,l]+sm[k]*Cpz[i,l]+sp[k]*Czm[l,i]+sm[i]*Cpz[k,l])
                dCzz[k,l] += (-1im*2*Ω[k,i]*(-sp[k]*Cmz[i,l] - sm[i]*Cpz[k,l] + sm[k]*Cpz[i,l] + sp[i]*Cmz[k,l])
                              -1im*2*Ω[l,i]*(-sp[l]*Cmz[i,k] - sm[i]*Czp[k,l] + sm[l]*Cpz[i,k] + sp[i]*Czm[k,l])
                              +Γ[k,i]*(sp[k]*Cmz[i,l] + sm[i]*Cpz[k,l] + sm[k]*Cpz[i,l] + sp[i]*Cmz[k,l])
                              +Γ[l,i]*(sp[l]*Cmz[i,k] + sm[i]*Czp[k,l] + sm[l]*Cpz[i,k] + sp[i]*Czm[k,l])
                              )
            end
        end
        combine_parameters(dy, dsz, dsp, dCpm, dCpz, dCpp, dCzz)
    end

    tout = Float64[]
    sp_out = Vector{Complex128}[]
    sz_out = Vector{Complex128}[]
    C_out = Dict{AbstractString, Matrix{Float64}}[]
    function fout(t, y::Vector{Complex128})
        sz, sp, Cpm, Cpz, Cpp, Czz = split_parameters(y)
        Cmp = transpose(Cpm)
        Cmm = conj(Cpp)
        Cmz = conj(Cpz)
        Czp = transpose(Cpz)
        Czm = ctranspose(Cpz)
        C = Dict{AbstractString, Matrix{Float64}}()
        C["xx"] = real(Cpp+Cpm+Cmp+Cmm)
        C["xy"] = -imag(Cpp-Cpm+Cmp-Cmm)
        C["xz"] = real(Cpz+Cmz)
        C["yx"] = -imag(Cpp-Cmp+Cpm-Cmm)
        C["yy"] = -real(Cpp-Cpm-Cmp+Cmm)
        C["yz"] = -imag(Cpz-Cmz)
        C["zx"] = real(Czp+Czm)
        C["zy"] = -imag(Czp-Czm)
        C["zz"] = real(Czz)
        push!(tout, t)
        push!(sz_out, sz)
        push!(sp_out, sp)
        push!(C_out, C)
    end

    quantumoptics.ode_dopri.ode(f, T, y0, fout=fout)
    return tout, sz_out, sp_out, C_out
end

type State
    sx::Vector{Float64}
    sy::Vector{Float64}
    sz::Vector{Float64}
end

function blochstate(phi, theta, N::Int=1)
    sx = ones(Float64, N)*cos(phi)*sin(theta)
    sy = ones(Float64, N)*sin(phi)*sin(theta)
    sz = ones(Float64, N)*cos(theta)
    return State(sx, sy, sz)
end

function meanfield_timeevolution2(T, system, sx0, sy0, sz0)
    N = length(system.atoms)
    Ω = OmegaMatrix(system)
    Γ = GammaMatrix(system)
    γ = system.gamma
    s0 = zeros(Float64, 3, N)
    s0[1,:] = sx0
    s0[2,:] = sy0
    s0[3,:] = sz0
    function f(t, s::Vector{Float64}, ds::Vector{Float64})
        s = reshape(s, 3, N)
        ds = reshape(ds, 3, N)
        for k=1:N
            ds[1,k] = -0.5*γ*s[1,k]
            ds[2,k] = -0.5*γ*s[2,k]
            ds[3,k] = γ*(1-s[3,k])
            for j=1:N
                if j==k
                    continue
                end
                ds[1,k] += Ω[k,j]*s[2,j]*s[3,k] - 0.5*Γ[k,j]*s[1,j]*s[3,k]
                ds[2,k] += -Ω[k,j]*s[1,j]*s[3,k] - 0.5*Γ[k,j]*s[2,j]*s[3,k]
                ds[3,k] += Ω[k,j]*(s[1,j]*s[2,k] - s[2,j]*s[1,k]) + 0.5*Γ[k,j]*(s[1,j]*s[1,k] + s[2,j]*s[2,k])
            end
        end
    end

    t_out = Float64[]
    sx_out = Vector{Float64}[]
    sy_out = Vector{Float64}[]
    sz_out = Vector{Float64}[]
    function fout(t, y::Vector{Float64})
        y = reshape(y, 3, N)
        push!(t_out, t)
        push!(sx_out, vec(y[1,:]))
        push!(sy_out, vec(y[2,:]))
        push!(sz_out, vec(y[3,:]))
    end

    quantumoptics.ode_dopri.ode(f, T, reshape(s0, 3*N), fout=fout)
    return t_out, sx_out, sy_out, sz_out
end


function correlation_timeevolution2(T, system, sx0, sy0, sz0, C0)
    N = length(system.atoms)
    Ω = OmegaMatrix(system)
    Γ = GammaMatrix(system)
    γ = system.gamma
    function split_parameters(y)
        sx = y[0*N+1:1*N]
        sy = y[1*N+1:2*N]
        sz = y[2*N+1:3*N]
        a = 3*N
        C = Dict{AbstractString, Matrix{Float64}}()
        C["xx"] = reshape(y[a+0*N^2+1:a+1*N^2], N, N)
        C["yy"] = reshape(y[a+1*N^2+1:a+2*N^2], N, N)
        C["zz"] = reshape(y[a+2*N^2+1:a+3*N^2], N, N)
        C["xy"] = reshape(y[a+3*N^2+1:a+4*N^2], N, N)
        C["xz"] = reshape(y[a+4*N^2+1:a+5*N^2], N, N)
        C["yz"] = reshape(y[a+5*N^2+1:a+6*N^2], N, N)
        return sx, sy, sz, C
    end
    function combine_parameters(y, sx, sy, sz, C)
        y[0*N+1:1*N] = sx
        y[1*N+1:2*N] = sy
        y[2*N+1:3*N] = sz
        a = 3*N
        y[a+0*N^2+1:a+1*N^2] = reshape(C["xx"], N^2)
        y[a+1*N^2+1:a+2*N^2] = reshape(C["yy"], N^2)
        y[a+2*N^2+1:a+3*N^2] = reshape(C["zz"], N^2)
        y[a+3*N^2+1:a+4*N^2] = reshape(C["xy"], N^2)
        y[a+4*N^2+1:a+5*N^2] = reshape(C["xz"], N^2)
        y[a+5*N^2+1:a+6*N^2] = reshape(C["yz"], N^2)
    end
    y0 = zeros(Float64, 3*N + 6*N^2)
    combine_parameters(y0, sx0, sy0, sz0, C0)
    function cor3(s, C, name, indices)
        name = Char[c for c=name]
        perm = sortperm(name)
        name = name[perm]
        indices = indices[perm]
        x = -2*s[name[1]][indices[1]]*s[name[2]][indices[2]]*s[name[3]][indices[3]]
        x += s[name[1]][indices[1]]*C[string(name[2], name[3])][indices[2], indices[3]]
        x += s[name[2]][indices[2]]*C[string(name[1], name[3])][indices[1], indices[3]]
        x += s[name[3]][indices[3]]*C[string(name[1], name[2])][indices[1], indices[2]]
        return x
    end

    function f(t, y::Vector{Float64}, dy::Vector{Float64})
        sx, sy, sz, C = split_parameters(y)
        s = Dict('x'=>sx, 'y'=>sy, 'z'=>sz)
        dsx, dsy, dsz, dC = split_parameters(zeros(Float64, 3*N + 6*N^2))
        for k=1:N
            dsx[k] = -0.5*γ*sx[k]
            dsy[k] = -0.5*γ*sy[k]
            dsz[k] = γ*(1-sz[k])
            for j=1:N
                if j==k
                    continue
                end
                dsx[k] += Ω[k,j]*C["yz"][j,k] - 0.5*Γ[k,j]*C["xz"][j,k]
                dsy[k] += -Ω[k,j]*C["xz"][j,k] - 0.5*Γ[k,j]*C["yz"][j,k]
                dsz[k] += Ω[k,j]*(C["xy"][j,k]-C["xy"][k,j]) + 0.5*Γ[k,j]*(C["xx"][j,k] + C["yy"][j,k])
            end
        end
        for k=1:N, l=1:N
            if k==l
                continue
            end
            dC["xx"][k,l] = -γ*C["xx"][k,l] + Γ[k,l]*(C["zz"][k,l]-0.5*sz[k]-0.5*sz[l])
            dC["yy"][k,l] = -γ*C["yy"][k,l] + Γ[k,l]*(C["zz"][k,l]-0.5*sz[k]-0.5*sz[l])
            dC["zz"][k,l] = -2*γ*C["zz"][k,l] + γ*(sz[k]+sz[l]) + Γ[k,l]*(C["yy"][k,l]+C["xx"][k,l])
            dC["xy"][k,l] = Ω[k,l]*(sz[k]-sz[l]) - γ*C["xy"][k,l]
            dC["xz"][k,l] = Ω[k,l]*sy[l] - 1.5*γ*C["xz"][k,l] + γ*sx[k] - Γ[k,l]*(C["xz"][l,k]-0.5*sx[l])
            dC["yz"][k,l] = -Ω[k,l]*sx[l] - 1.5*γ*C["yz"][k,l] + γ*sy[k] - Γ[k,l]*(C["yz"][l,k]-0.5*sy[l])
            for j=1:N
                if j==l || j==k
                    continue
                end
                c(name) = cor3(s, C, name, [k,l,j])
                dC["xx"][k,l] += Ω[k,j]*c("zxy") + Ω[l,j]*c("xzy") - 0.5*Γ[k,j]*c("zxx") - 0.5*Γ[l,j]*c("xzx")
                dC["yy"][k,l] += -Ω[k,j]*c("zyx") - Ω[l,j]*c("yzx") - 0.5*Γ[k,j]*c("zyy") - 0.5*Γ[l,j]*c("yzy")
                dC["zz"][k,l] += (Ω[k,j]*(c("yzx")-c("xzy")) + Ω[l,j]*(c("zyx")-c("zxy"))
                                + 0.5*Γ[k,j]*(c("xzx")+c("yzy")) + 0.5*Γ[l,j]*(c("zxx")+c("zyy")))
                dC["xy"][k,l] += Ω[k,j]*c("zyy") - Ω[l,j]*c("xzx") - 0.5*Γ[k,j]*c("zyx") - 0.5*Γ[l,j]*c("xzy")
                dC["xz"][k,l] += (Ω[k,j]*c("zzy") + Ω[l,j]*(c("xyx")-c("xxy"))
                                - 0.5*Γ[k,j]*c("zzx") + 0.5*Γ[l,j]*(c("xxx")+c("xyy")))
                dC["yz"][k,l] += (-Ω[k,j]*c("zzx") + Ω[l,j]*(c("yyx")-c("yxy"))
                                - 0.5*Γ[k,j]*c("zzy") + 0.5*Γ[l,j]*(c("yxx")+c("yyy")))
            end
        end
        combine_parameters(dy, dsx, dsy, dsz, dC)
    end

    t_out = Float64[]
    sx_out = Vector{Float64}[]
    sy_out = Vector{Float64}[]
    sz_out = Vector{Float64}[]
    C_out = Dict{AbstractString, Matrix{Float64}}[]
    function fout(t, y::Vector{Float64})
        sx, sy, sz, C = split_parameters(y)
        push!(t_out, t)
        push!(sx_out, sx)
        push!(sy_out, sy)
        push!(sz_out, sz)
        #push!(C_out, C)
    end

    quantumoptics.ode_dopri.ode(f, T, y0, fout=fout)
    return t_out, sx_out, sy_out, sz_out, C_out
end

function correlation_timeevolution3(T, system, sx0, sy0, sz0, C0)
    N = length(system.atoms)
    Ω = OmegaMatrix(system)
    Γ = GammaMatrix(system)
    γ = system.gamma
    function split_parameters(y)
        sx = y[0*N+1:1*N]
        sy = y[1*N+1:2*N]
        sz = y[2*N+1:3*N]
        a = 3*N
        Cxx = reshape(y[a+0*N^2+1:a+1*N^2], N, N)
        Cyy = reshape(y[a+1*N^2+1:a+2*N^2], N, N)
        Czz = reshape(y[a+2*N^2+1:a+3*N^2], N, N)
        Cxy = reshape(y[a+3*N^2+1:a+4*N^2], N, N)
        Cxz = reshape(y[a+4*N^2+1:a+5*N^2], N, N)
        Cyz = reshape(y[a+5*N^2+1:a+6*N^2], N, N)
        return sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz
    end
    function combine_parameters(y, sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz)
        y[0*N+1:1*N] = sx
        y[1*N+1:2*N] = sy
        y[2*N+1:3*N] = sz
        a = 3*N
        y[a+0*N^2+1:a+1*N^2] = reshape(Cxx, N^2)
        y[a+1*N^2+1:a+2*N^2] = reshape(Cyy, N^2)
        y[a+2*N^2+1:a+3*N^2] = reshape(Czz, N^2)
        y[a+3*N^2+1:a+4*N^2] = reshape(Cxy, N^2)
        y[a+4*N^2+1:a+5*N^2] = reshape(Cxz, N^2)
        y[a+5*N^2+1:a+6*N^2] = reshape(Cyz, N^2)
    end
    y0 = zeros(Float64, 3*N + 6*N^2)
    combine_parameters(y0, sx0, sy0, sz0, C0["xx"], C0["yy"], C0["zz"], C0["xy"], C0["xz"], C0["yz"])
    function correlation(s1::Float64, s2::Float64, s3::Float64, C12::Float64, C13::Float64, C23::Float64)
        return -2.*s1*s2*s3 + s1*C23 + s2*C13 + s3*C12
    end

    function f(t, y::Vector{Float64}, dy::Vector{Float64})
        sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = split_parameters(y)
        dsx, dsy, dsz, dCxx, dCyy, dCzz, dCxy, dCxz, dCyz = split_parameters(zeros(Float64, 3*N + 6*N^2))
        Cyx = transpose(Cxy)
        Czx = transpose(Cxz)
        Czy = transpose(Cyz)
        for k=1:N
            dsx[k] = -0.5*γ*sx[k]
            dsy[k] = -0.5*γ*sy[k]
            dsz[k] = γ*(1-sz[k])
            for j=1:N
                if j==k
                    continue
                end
                dsx[k] += Ω[k,j]*Cyz[j,k] - 0.5*Γ[k,j]*Cxz[j,k]
                dsy[k] += -Ω[k,j]*Cxz[j,k] - 0.5*Γ[k,j]*Cyz[j,k]
                dsz[k] += Ω[k,j]*(Cxy[j,k]-Cxy[k,j]) + 0.5*Γ[k,j]*(Cxx[j,k] + Cyy[j,k])
            end
        end
        for k=1:N, l=1:N
            if k==l
                continue
            end
            dCxx[k,l] = -γ*Cxx[k,l] + Γ[k,l]*(Czz[k,l]-0.5*sz[k]-0.5*sz[l])
            dCyy[k,l] = -γ*Cyy[k,l] + Γ[k,l]*(Czz[k,l]-0.5*sz[k]-0.5*sz[l])
            dCzz[k,l] = -2*γ*Czz[k,l] + γ*(sz[k]+sz[l]) + Γ[k,l]*(Cyy[k,l]+Cxx[k,l])
            dCxy[k,l] = Ω[k,l]*(sz[k]-sz[l]) - γ*Cxy[k,l]
            dCxz[k,l] = Ω[k,l]*sy[l] - 1.5*γ*Cxz[k,l] + γ*sx[k] - Γ[k,l]*(Cxz[l,k]-0.5*sx[l])
            dCyz[k,l] = -Ω[k,l]*sx[l] - 1.5*γ*Cyz[k,l] + γ*sy[k] - Γ[k,l]*(Cyz[l,k]-0.5*sy[l])
            for j=1:N
                if j==l || j==k
                    continue
                end
                Cxxx = correlation(sx[k], sx[l], sx[j], Cxx[k,l], Cxx[k,j], Cxx[l,j])
                Cxxy = correlation(sx[k], sx[l], sy[j], Cxx[k,l], Cxy[k,j], Cxy[l,j])
                Cxyx = correlation(sx[k], sy[l], sx[j], Cxy[k,l], Cxx[k,j], Cyx[l,j])
                Cxyy = correlation(sx[k], sy[l], sy[j], Cxy[k,l], Cxy[k,j], Cyy[l,j])
                Cxzx = correlation(sx[k], sz[l], sx[j], Cxz[k,l], Cxx[k,j], Czx[l,j])
                Cxzy = correlation(sx[k], sz[l], sy[j], Cxz[k,l], Cxy[k,j], Czy[l,j])
                Cyxx = correlation(sy[k], sx[l], sx[j], Cyx[k,l], Cyx[k,j], Cxx[l,j])
                Cyxy = correlation(sy[k], sx[l], sy[j], Cyx[k,l], Cyy[k,j], Cxy[l,j])
                Cyyx = correlation(sy[k], sy[l], sx[j], Cyy[k,l], Cyx[k,j], Cyx[l,j])
                Cyyy = correlation(sy[k], sy[l], sy[j], Cyy[k,l], Cyy[k,j], Cyy[l,j])
                Cyzx = correlation(sy[k], sz[l], sx[j], Cyz[k,l], Cyx[k,j], Czx[l,j])
                Cyzy = correlation(sy[k], sz[l], sy[j], Cyz[k,l], Cyy[k,j], Czy[l,j])
                Czxx = correlation(sz[k], sx[l], sx[j], Czx[k,l], Czx[k,j], Cxx[l,j])
                Czxy = correlation(sz[k], sx[l], sy[j], Czx[k,l], Czy[k,j], Cxy[l,j])
                Czyx = correlation(sz[k], sy[l], sx[j], Czy[k,l], Czx[k,j], Cyx[l,j])
                Czyy = correlation(sz[k], sy[l], sy[j], Czy[k,l], Czy[k,j], Cyy[l,j])
                Czzx = correlation(sz[k], sz[l], sx[j], Czz[k,l], Czx[k,j], Czx[l,j])
                Czzy = correlation(sz[k], sz[l], sy[j], Czz[k,l], Czy[k,j], Czy[l,j])

                dCxx[k,l] += Ω[k,j]*Czxy + Ω[l,j]*Cxzy - 0.5*Γ[k,j]*Czxx - 0.5*Γ[l,j]*Cxzx
                dCyy[k,l] += -Ω[k,j]*Czyx - Ω[l,j]*Cyzx - 0.5*Γ[k,j]*Czyy - 0.5*Γ[l,j]*Cyzy
                dCzz[k,l] += (Ω[k,j]*(Cyzx-Cxzy) + Ω[l,j]*(Czyx-Czxy)
                                + 0.5*Γ[k,j]*(Cxzx+Cyzy) + 0.5*Γ[l,j]*(Czxx+Czyy))
                dCxy[k,l] += Ω[k,j]*Czyy - Ω[l,j]*Cxzx - 0.5*Γ[k,j]*Czyx - 0.5*Γ[l,j]*Cxzy
                dCxz[k,l] += (Ω[k,j]*Czzy + Ω[l,j]*(Cxyx-Cxxy)
                                - 0.5*Γ[k,j]*Czzx + 0.5*Γ[l,j]*(Cxxx+Cxyy))
                dCyz[k,l] += (-Ω[k,j]*Czzx + Ω[l,j]*(Cyyx-Cyxy)
                                - 0.5*Γ[k,j]*Czzy + 0.5*Γ[l,j]*(Cyxx+Cyyy))
            end
        end
        combine_parameters(dy, dsx, dsy, dsz, dCxx, dCyy, dCzz, dCxy, dCxz, dCyz)
    end

    t_out = Float64[]
    sx_out = Vector{Float64}[]
    sy_out = Vector{Float64}[]
    sz_out = Vector{Float64}[]
    C_out = Dict{AbstractString, Matrix{Float64}}[]
    function fout(t, y::Vector{Float64})
        sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = split_parameters(y)
        push!(t_out, t)
        push!(sx_out, sx)
        push!(sy_out, sy)
        push!(sz_out, sz)
        #push!(C_out, C)
    end

    quantumoptics.ode_dopri.ode(f, T, y0, fout=fout)
    return t_out, sx_out, sy_out, sz_out, C_out
end

end # module
