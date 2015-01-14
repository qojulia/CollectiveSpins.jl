module mpc

using ArrayViews
module Optim
    try
        require("Optim")
        global optimize = Main.Optim.optimize
    catch
        println("Optim package not available. (Needed for calculation of squeezing parameter)")
    end
end
using quantumoptics
using ..interaction, ..system, ..meanfield, ..quantum


function blochstate(phi, theta, N::Int=1)
    state = zeros(Float64, (3*N)*(2*N+1))
    sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = splitstate(state)
    for i=1:N
        sx[i] = cos(phi)*sin(theta)
        sy[i] = sin(phi)*sin(theta)
        sz[i] = cos(theta)
    end
    for i=1:N, j=1:N
        if i==j
            continue
        end
        Cxx[i,j] = sx[i]*sx[j]
        Cyy[i,j] = sy[i]*sy[j]
        Czz[i,j] = sz[i]*sz[j]
        Cxy[i,j] = sx[i]*sy[j]
        Cxz[i,j] = sx[i]*sz[j]
        Cyz[i,j] = sy[i]*sz[j]
    end
    return vec(state)
end

function integersqrt(N::Int)
    n = sqrt(N)
    if abs(int(n)-n)>10*eps(n)
        error("N is not a square of an integer.")
    end
    return int(n)
end

function dim(state::Vector{Float64})
    x, rem = divrem(length(state), 3)
    @assert rem==0
    N, rem = divrem(-1 + integersqrt(1+8*x), 4)
    @assert rem==0
    return N
end

function splitstate(state::Vector{Float64})
    N = dim(state)
    state = reshape(state, 3*N, 2*N+1)
    sx = view(state, 0*N+1:1*N, 2*N+1)
    sy = view(state, 1*N+1:2*N, 2*N+1)
    sz = view(state, 2*N+1:3*N, 2*N+1)
    Cxx = view(state, 0*N+1:1*N, 0*N+1:1*N)
    Cyy = view(state, 1*N+1:2*N, 0*N+1:1*N)
    Czz = view(state, 2*N+1:3*N, 0*N+1:1*N)
    Cxy = view(state, 0*N+1:1*N, 1*N+1:2*N)
    Cxz = view(state, 1*N+1:2*N, 1*N+1:2*N)
    Cyz = view(state, 2*N+1:3*N, 1*N+1:2*N)
    return sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz
end

function covarianceoperator(productstate::Vector{Operator}, operators::Vector{Operator}, indices::Vector{Int})
    x = Operator[(i in indices ? operators[findfirst(indices, i)] : productstate[i]) for i=1:length(productstate)]
    return tensor(x...)
end

function correlation2covariance(corstate::Vector{Float64})
    N = dim(corstate)
    covstate = zeros(Float64, size(corstate)...)
    sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = splitstate(corstate)
    covsx, covsy, covsz, Covxx, Covyy, Covzz, Covxy, Covxz, Covyz = splitstate(covstate)
    for k=1:N
        covsx[k] = sx[k]
        covsy[k] = sy[k]
        covsz[k] = sz[k]
    end
    for k=1:N, l=1:N
        if k==l
            continue
        end
        Covxx[k,l] = Cxx[k,l] - sx[k]*sx[l]
        Covyy[k,l] = Cyy[k,l] - sy[k]*sy[l]
        Covzz[k,l] = Czz[k,l] - sz[k]*sz[l]
        Covxy[k,l] = Cxy[k,l] - sx[k]*sy[l]
        Covxz[k,l] = Cxz[k,l] - sx[k]*sz[l]
        Covyz[k,l] = Cyz[k,l] - sy[k]*sz[l]
    end
    return covstate
end

function covariance2correlation(covstate::Vector{Float64})
    N = dim(covstate)
    corstate = zeros(Float64, size(covstate)...)
    sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = splitstate(corstate)
    covsx, covsy, covsz, Covxx, Covyy, Covzz, Covxy, Covxz, Covyz = splitstate(covstate)
    for k=1:N
        sx[k] = covsx[k]
        sy[k] = covsy[k]
        sz[k] = covsz[k]
    end
    for k=1:N, l=1:N
        if k==l
            continue
        end
        Cxx[k,l] = Covxx[k,l] + sx[k]*sx[l]
        Cyy[k,l] = Covyy[k,l] + sy[k]*sy[l]
        Czz[k,l] = Covzz[k,l] + sz[k]*sz[l]
        Cxy[k,l] = Covxy[k,l] + sx[k]*sy[l]
        Cxz[k,l] = Covxz[k,l] + sx[k]*sz[l]
        Cyz[k,l] = Covyz[k,l] + sy[k]*sz[l]
    end
    return corstate
end

function densityoperator(state::Vector{Float64})
    N = dim(state)
    covstate = correlation2covariance(state)
    sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = splitstate(covstate)
    productstate = Operator[meanfield.densityoperator(sx[k], sy[k], sz[k]) for k=1:N]
    C(op1,op2,index1,index2) = covarianceoperator(productstate, [op1,op2], [index1,index2])
    ρ = reduce(tensor, productstate)
    for k=1:N, l=k+1:N
        ρ += 0.25*(
              Cxx[k,l]*C(sigmax,sigmax,k,l) + Cxy[l,k]*C(sigmay,sigmax,k,l) + Cxz[l,k]*C(sigmaz,sigmax,k,l)
            + Cxy[k,l]*C(sigmax,sigmay,k,l) + Cyy[k,l]*C(sigmay,sigmay,k,l) + Cyz[l,k]*C(sigmaz,sigmay,k,l)
            + Cxz[k,l]*C(sigmax,sigmaz,k,l) + Cyz[k,l]*C(sigmay,sigmaz,k,l) + Czz[k,l]*C(sigmaz,sigmaz,k,l))
    end
    return ρ
end

mpcstate(N::Int) = zeros(Float64, (3*N)*(2*N+1))

function mpcstate(rho::Operator)
    N = quantum.dim(rho)
    basis = rho.basis_l
    state = mpcstate(N)
    sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = splitstate(state)
    f(ind, op) = real(expect(embed(basis, ind, op), rho))
    for k=1:N
        sx[k] = f(k, sigmax)
        sy[k] = f(k, sigmay)
        sz[k] = f(k, sigmaz)
        for l=1:N
            if k==l
                continue
            end
            Cxx[k,l] = f([k,l], [sigmax, sigmax])
            Cyy[k,l] = f([k,l], [sigmay, sigmay])
            Czz[k,l] = f([k,l], [sigmaz, sigmaz])
            Cxy[k,l] = f([k,l], [sigmax, sigmay])
            Cxz[k,l] = f([k,l], [sigmax, sigmaz])
            Cyz[k,l] = f([k,l], [sigmay, sigmaz])
        end
    end
    return state
end

sx(state::Vector{Float64}) = begin N=dim(state); view(reshape(state, 3*N, 2*N+1), 0*N+1:1*N, 2*N+1) end
sy(state::Vector{Float64}) = begin N=dim(state); view(reshape(state, 3*N, 2*N+1), 1*N+1:2*N, 2*N+1) end
sz(state::Vector{Float64}) = begin N=dim(state); view(reshape(state, 3*N, 2*N+1), 2*N+1:3*N, 2*N+1) end

function correlation(s1::Float64, s2::Float64, s3::Float64, C12::Float64, C13::Float64, C23::Float64)
    return -2.*s1*s2*s3 + s1*C23 + s2*C13 + s3*C12
end

function timeevolution(T, S::system.SpinCollection, state0::Vector{Float64}; fout=nothing)
    N = length(S.spins)
    Ω = interaction.OmegaMatrix(S)
    Γ = interaction.GammaMatrix(S)
    γ = S.gamma

    function f(t, y::Vector{Float64}, dy::Vector{Float64})
        sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = splitstate(y)
        dsx, dsy, dsz, dCxx, dCyy, dCzz, dCxy, dCxz, dCyz = splitstate(dy)
        @inbounds for k=1:N
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
        @inbounds for k=1:N, l=1:N
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
                Cxyx = correlation(sx[k], sy[l], sx[j], Cxy[k,l], Cxx[k,j], Cxy[j,l])
                Cxyy = correlation(sx[k], sy[l], sy[j], Cxy[k,l], Cxy[k,j], Cyy[l,j])
                Cxzx = correlation(sx[k], sz[l], sx[j], Cxz[k,l], Cxx[k,j], Cxz[j,l])
                Cxzy = correlation(sx[k], sz[l], sy[j], Cxz[k,l], Cxy[k,j], Cyz[j,l])
                Cyxx = correlation(sy[k], sx[l], sx[j], Cxy[l,k], Cxy[j,k], Cxx[l,j])
                Cyxy = correlation(sy[k], sx[l], sy[j], Cxy[l,k], Cyy[k,j], Cxy[l,j])
                Cyyx = correlation(sy[k], sy[l], sx[j], Cyy[k,l], Cxy[j,k], Cxy[j,l])
                Cyyy = correlation(sy[k], sy[l], sy[j], Cyy[k,l], Cyy[k,j], Cyy[l,j])
                Cyzx = correlation(sy[k], sz[l], sx[j], Cyz[k,l], Cxy[j,k], Cxz[j,l])
                Cyzy = correlation(sy[k], sz[l], sy[j], Cyz[k,l], Cyy[k,j], Cyz[j,l])
                Czxx = correlation(sz[k], sx[l], sx[j], Cxz[l,k], Cxz[j,k], Cxx[l,j])
                Czxy = correlation(sz[k], sx[l], sy[j], Cxz[l,k], Cyz[j,k], Cxy[l,j])
                Czyx = correlation(sz[k], sy[l], sx[j], Cyz[l,k], Cxz[j,k], Cxy[j,l])
                Czyy = correlation(sz[k], sy[l], sy[j], Cyz[l,k], Cyz[j,k], Cyy[l,j])
                Czzx = correlation(sz[k], sz[l], sx[j], Czz[k,l], Cxz[j,k], Cxz[j,l])
                Czzy = correlation(sz[k], sz[l], sy[j], Czz[k,l], Cyz[j,k], Cyz[j,l])

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
    end

    if fout==nothing
        t_out = Float64[]
        state_out = Vector{Float64}[]
        function fout_(t, state::Vector{Float64})
            push!(t_out, t)
            push!(state_out, deepcopy(state))
        end

        quantumoptics.ode_dopri.ode(f, T, state0, fout=fout_)
        return t_out, state_out
    else
        return quantumoptics.ode_dopri.ode(f, T, state0, fout=fout)
    end
end

# function rotate(rotationaxis::Vector{Float64}, angles::Vector{Float64}, state::Vector{Float64})
#     N = dim(state)
#     @assert length(rotationaxis)==3
#     @assert length(angles)==N
#     w = rotationaxis/norm(rotationaxis)
#     covstate = correlation2covariance(state)
#     sx, sy, sz, Covxx, Covyy, Covzz, Covxy, Covxz, Covyz = splitstate(covstate)
#     for i=1:N
#         v = [sx[i], sy[i], sz[i]]
#         θ = angles[i]
#         sx[i], sy[i], sz[i] = cos(θ)*v + sin(θ)*(w × v) + (1-cos(θ))*(w ⋅ v)*w
#     end
#     return covariance2correlation(covstate)
# end

function rotate(rotationaxis::Vector{Float64}, angles::Vector{Float64}, state::Vector{Float64})
    N = dim(state)
    @assert length(rotationaxis)==3
    @assert length(angles)==N
    w = rotationaxis/norm(rotationaxis)
    S_total = splitstate(state)
    rotstate = deepcopy(state)
    S_total_rot = splitstate(rotstate)
    pstate = zeros(Float64, (3*2)*(2*2+1))
    S_2 = splitstate(pstate)
    for k=1:N,l=1:N
        if k==l
            continue
        end
        for i=1:3
            S_2[i][1] = S_total[i][k]
            S_2[i][2] = S_total[i][l]
        end
        for i=4:9
            S_2[i][1,2] = S_total[i][k,l]
            S_2[i][2,1] = S_total[i][l,k]
        end
        rho_p = densityoperator(pstate)
        rho_p_rot = quantum.rotate(rotationaxis, Float64[angles[k], angles[l]], rho_p)
        pstate_rot = mpcstate(rho_p_rot)
        S_2rot = splitstate(pstate_rot)
        for i=1:3
            S_total_rot[i][k] = S_2rot[i][1]
            S_total_rot[i][l] = S_2rot[i][2]
        end
        for i=4:9
            S_total_rot[i][k,l] = S_2rot[i][1,2]
            S_total_rot[i][l,k] = S_2rot[i][2,1]
        end
    end
    return rotstate
end

rotate(rotationaxis::Vector{Float64}, angle::Float64, state::Vector{Float64}) = rotate(rotationaxis, ones(Float64, dim(state))*angle, state)

function squeeze_sx(χT::Float64, state0::Vector{Float64})
    T = [0,1.]
    N = dim(state0)
    χeff = 4*χT/N^2
    function f(t, y::Vector{Float64}, dy::Vector{Float64})
        sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = splitstate(y)
        dsx, dsy, dsz, dCxx, dCyy, dCzz, dCxy, dCxz, dCyz = splitstate(dy)
        for k=1:N
            dsx[k] = 0.
            dsy[k] = 0.
            dsz[k] = 0.
            for j=1:N
                if j==k
                    continue
                end
                dsy[k] += -χeff*Cxz[j,k]
                dsz[k] += χeff*Cxy[j,k]
            end
        end
        for k=1:N, l=1:N
            if k==l
                continue
            end
            dCxx[k,l] = 0.
            dCyy[k,l] = 0.
            dCzz[k,l] = 0.
            dCxy[k,l] = -χeff*sz[l]
            dCxz[k,l] = χeff*sy[l]
            dCyz[k,l] = 0.
            for j=1:N
                if j==l || j==k
                    continue
                end
                Cxyx = correlation(sx[k], sy[l], sx[j], Cxy[k,l], Cxx[k,j], Cxy[j,l])
                Cxzx = correlation(sx[k], sz[l], sx[j], Cxz[k,l], Cxx[k,j], Cxz[j,l])
                Cyyx = correlation(sy[k], sy[l], sx[j], Cyy[k,l], Cxy[j,k], Cxy[j,l])
                Cyzx = correlation(sy[k], sz[l], sx[j], Cyz[k,l], Cxy[j,k], Cxz[j,l])
                Czyx = correlation(sz[k], sy[l], sx[j], Cyz[l,k], Cxz[j,k], Cxy[j,l])
                Czzx = correlation(sz[k], sz[l], sx[j], Czz[k,l], Cxz[j,k], Cxz[j,l])

                dCyy[k,l] += -χeff*(Czyx+Cyzx)
                dCzz[k,l] += χeff*(Cyzx+Czyx)
                dCxy[k,l] += -χeff*Cxzx
                dCxz[k,l] += χeff*Cxyx
                dCyz[k,l] += χeff*(Czzx-Cyyx)
            end
        end
    end

    state_out = Vector{Float64}[]
    function fout_(t, state::Vector{Float64})
        push!(state_out, deepcopy(state))
    end

    quantumoptics.ode_dopri.ode(f, T, state0, fout=fout_)
    return state_out[end]
end

function orthogonal_vectors(n::Vector{Float64})
    n = n/norm(n)
    v = (n[1]<n[2] ? [1.,0.,0.] : [0.,1.,0.])
    e1 = v - dot(n,v)*n
    e1 = e1/norm(e1)
    e2 = cross(n, e1)
    e2 = e2/norm(e2)
    return e1, e2
end

function squeezingparameter(state::Vector{Float64})
    N = dim(state)
    sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz = map(sum, splitstate(state))
    n = 1./N*[sx, sy, sz]
    e1, e2 = orthogonal_vectors(n)
    function f(phi)
        nphi = cos(phi)*e1 + sin(phi)*e2
        nx = nphi[1]
        ny = nphi[2]
        nz = nphi[3]
        Sphi2 = 1./N^2*(N+nx*nx*Cxx + ny*ny*Cyy + nz*nz*Czz +
                2*nx*ny*Cxy + 2*nx*nz*Cxz + 2*ny*nz*Cyz)
        return Sphi2
    end
    varSmin = Optim.optimize(f, 0., 2.pi).f_minimum
    return sqrt(N*varSmin)/norm(n)
end

end # module
