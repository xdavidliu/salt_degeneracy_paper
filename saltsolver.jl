type Laser{T<:AbstractFloat}
    ɛ::Vector{Complex{T}}
    F::Vector{T}
    ωa::T
    γ⟂::T

    Laser(ɛ, F, ωa, γ⟂) = begin
        length(ɛ) == length(F) || error("incorrect length(F)")
        new(ɛ, F, ωa, γ⟂)
    end
end
Laser(ɛ, F, ωa, γ⟂) = Laser{typeof(ωa)}(ɛ, F, ωa, γ⟂)

function normalize!(E::AbstractVector)
    imax = indmax(abs(E))
    imax != 0 && scale!(E, 1/E[imax])
    imax
end

type Mode{T<:AbstractFloat, V<:Number, U<:Integer}
    E::Vector{Complex{T}}
    ω::V
    c²::T
    imax::U

    Mode(E, ω, c²) = begin
        imax = normalize!(E)
        new(E, ω, c², imax)
    end
    # force imax
    Mode(E, ω, c², imax) = begin
        scale!(E, 1/E[imax])
        new(E, ω, c², imax)
    end
end
Mode(E, ω, c²) = Mode{typeof(c²), typeof(ω), Int}(E, ω, c²)
Mode(E, ω, c², imax) = Mode{typeof(c²), typeof(ω), typeof(imax)}(E, ω, c², imax)

function salt_operator!(md::Mode, las::Laser, D::Real, J::AbstractMatrix)
    γ = las.γ⟂ / (md.ω - las.ωa + im*las.γ⟂)
    N = length(md.E)
    for i=1:N
        H = 1/(1 + md.c²*abs(md.E[i])^2)
        z = md.ω^2 * (las.ɛ[i] + D*γ*H*las.F[i])
        J[i,i]     +=  real(z)
        J[i+N,i+N] +=  real(z)
        J[i,i+N]   += -imag(z)
        J[i+N,i]   +=  imag(z)
    end
end

function copy_component!(v::AbstractVector, E::AbstractVector)
    N = length(E)
    length(v) == 2N || error("incorrect length(v)")
    for i=1:N
        v[i]   = real(E[i])
        v[i+N] = imag(E[i])
    end
end

function residual!(laplacian!::Function, md::Mode, las::Laser, D::Real, J::AbstractMatrix, f::AbstractVector)
    fill!(J, 0.0)
    fill!(f, 0.0)
    N = length(md.E)
    laplacian!(sub(J, 1:N, 1:N))
    laplacian!(sub(J, N+1:2N, N+1:2N))
    salt_operator!(md, las, D, J)

    # use last column of J as scratch space
    v = sub(J, 1:2N, 2N+2)
    copy_component!(v, md.E)
    A_mul_B!(sub(f,1:2N), sub(J,1:2N,1:2N), v)
    f[end-1] = real(md.E[md.imax]) - 1
    f[end]   = imag(md.E[md.imax])
    norm(f)
end

function derivative!(md::Mode, las::Laser, D::Real, J::AbstractMatrix)
    N = length(md.E)
    J[end-1, md.imax] = 1
    J[end, md.imax+N] = 1

    ∂f∂ω  = sub(J, 1:2N, 2N+1)
    ∂f∂c² = sub(J, 1:2N, 2N+2)
    γ = las.γ⟂ / (md.ω - las.ωa + im*las.γ⟂)
    ω²D∂γ∂ω = md.ω^2*D*(-γ^2/las.γ⟂)
    for i=1:N
        # columns must be set using =
        H = 1/(1 + md.c²*abs(md.E[i])^2)
        ∂M∂ω = 2md.ω*(las.ɛ[i]+D*γ*las.F[i]*H) + ω²D∂γ∂ω*las.F[i]*H
        ∂f∂ω[i]   = real(∂M∂ω * md.E[i])
        ∂f∂ω[i+N] = imag(∂M∂ω * md.E[i])

        ω²γDFH²E = md.ω^2*γ*D*las.F[i]*H^2*md.E[i]
        ∂M∂c²E = -abs(md.E[i])^2*ω²γDFH²E
        ∂f∂c²[i]   = real(∂M∂c²E)
        ∂f∂c²[i+N] = imag(∂M∂c²E)

        # blocks must be incremented using +=
        ∂M∂E = -2md.c²*ω²γDFH²E
        J[i,i]     += real(∂M∂E) * real(md.E[i])
        J[i,i+N]   += real(∂M∂E) * imag(md.E[i])
        J[i+N,i]   += imag(∂M∂E) * real(md.E[i])
        J[i+N,i+N] += imag(∂M∂E) * imag(md.E[i])
    end
end

function update!(md::Mode, dv::AbstractVector)
    md.ω  += dv[end-1]
    md.c² += dv[end]
    N = length(md.E)
    for i=1:N
        md.E[i] += dv[i] + im*dv[i+N]
    end
end

type MaxIterationsReached <: Exception end

function solve!(laplacian!::Function, md::Mode, las::Laser, D::Real;
                maxits=9, tol=1e-8, isprint=false)
    D != 0.0 || error("D=0 unsolvable")
    N = length(md.E)
    length(las.ɛ) == N || error("incorrect length(ɛ)")
    J = zeros(2N+2, 2N+2)
    f = zeros(2N+2)

    for its=0:maxits
        res = residual!(laplacian!, md, las, D, J, f)
        isprint && println("|f| = ", res)
        res < tol && return (res, its)
        its < maxits || throw(MaxIterationsReached())
        derivative!(md, las, D, J)
        dv = scale!(J\f, -1)
        update!(md, dv)
    end
end

function periodic!(M::AbstractMatrix, h::Real)
    vals = [1, -2, 1]/h^2
    N = size(M, 1)
    rows(i::Int) = [1+mod(i-2, N), i, 1+mod(i, N)]
    for i=1:N
        M[rows(i), i] = vals
    end
end

function ring_eigenpair(N::Integer, ℓ::Integer, h::Real)
    (0 < abs(ℓ) <= (N-1)÷2) || error("incorrect ℓ")
    ω = √(2-2cos(2π*ℓ/N))/h
    E = exp(-2π*im*ℓ/N*(1:N))
    E, ω
end

# solution created with ring_eigenpair transforms
# with eigenvalue exp(2πiℓ/N) under 1/N rotations
function rotate(E::AbstractVector, frac::Rational)
    shift = frac*length(E)
    den(shift) == 1 || error("incorrect rotation")
    circshift(E, num(shift))
end

function get_ell(E::AbstractVector, nsym::Integer)
    imax = indmax(abs(E))
    rE = rotate(E, 1//nsym)
    λ = rE[imax] / E[imax]
    ℓ = angle(λ) / (2π/nsym)
    ℓ, norm(rE-λ*E) / norm(E)
end

function project_chiral(E::AbstractVector, nsym::Integer, ℓ::Integer)
    Echiral = zeros(length(E))
    for m=1:nsym
        z = exp(-2π*im*ℓ*m/nsym)
        Echiral += z/nsym * rotate(E, m//nsym)
        # 1/nsym factor means project(Eℓ) = Eℓ
    end
    Echiral
end

function coefficients{T<:AbstractVector}(E::AbstractVector, basis::Tuple{T, T})
    overlap(F) = dot(F, E) / norm(F)^2
    map(overlap, basis)
end

function cn_profile(F::AbstractVector, n::Integer)
    len = length(F)
    prof = zeros(eltype(F), len*n)
    for i=1:n
        start = (i-1)*len+1
        prof[start:start+len-1] = F
    end
    prof
end

function ring_flip(F::AbstractVector)
    vcat(F[1:1], flipdim(F[2:end], 1))
end

# must use this because ring_flip then normalize!
# will change c² if imax changes
function mode_flip(md::Mode)
    N = length(md.E)
    imax = (N-md.imax+1)%N+1
    Mode(ring_flip(md.E), md.ω, md.c², imax)
end

# simple root-finding using (extra/inter)polation
function root{T<:Number}(func::Function, xs::Tuple{T, T}; tol=1e-10, maxits=10, isprint=false)
    ys = zeros(T, 2)
    for (i,x) in enumerate(xs)
        ys[i] = func(x)
        abs(ys[i]) < tol && return x
    end
    ifar = indmax(abs(ys))
    iclose = 3-ifar
    xclose, xfar = xs[[iclose, ifar]]
    yclose, yfar = ys[[iclose, ifar]]

    for its=1:maxits
        abs(yclose-yfar) < 1e-9 && error("small denominator in root function")
        xnext = xclose - (xclose-xfar)/(yclose-yfar)*yclose
        ynext = func(xnext)
        if isprint
            println("xnext = ", xnext, ", ynext = ", ynext)
        end
        abs(ynext) < tol && return xnext
        xfar, yfar = xclose, yclose
        xclose, yclose = xnext, ynext
    end
    throw(MaxIterationsReached())
end

function threshold!{T<:Real}(laplacian!::Function, md::Mode, las::Laser, Ds::Tuple{T, T};
                    init_incr=0.05, tol=1e-10, maxits=10, newtontol=1e-8, newton_maxits=9, isprint=false)
    func(D) = begin
        solve!(laplacian!, md, las, D, tol=newtontol, maxits=newton_maxits)
        md.c²
    end
    root(func, Ds, tol=tol, maxits=maxits, isprint=isprint)
end

function overlap_integrals{U<:AbstractVector}(Es::Tuple{U, U}, ωt::Real, Dt::Real, las::Laser, nsym::Integer, ℓ::Integer, Lcav::Real)
    h = Lcav / length(Es[1])
    Gɛ = h*sum(las.ɛ .* Es[1] .* Es[2])
    GD = h*Dt*sum(las.F .* Es[1] .* Es[2])
    I = zeros(Complex{Float64}, 2)
    J = copy(I)
    K = copy(I)
    other(i) = 3-i
    for a in (1, 2)
        I[a] = h*Dt*sum(las.F .* abs(Es[a]).^2 .* Es[1] .* Es[2])
        J[a] = h*Dt*sum(las.F .* Es[a].^2 .* conj(Es[a]) .* Es[other(a)])
        K[a] = h*Dt*sum(las.F .* Es[a].^2 .* conj(Es[other(a)]) .* Es[a])
    end
    γ = las.γ⟂ / (ωt - las.ωa + im*las.γ⟂)
    H = 2Gɛ/(ωt*γ) + GD*(2/ωt - γ/las.γ⟂)
    H, I, J, K, GD
end

function standing_roots_test(H, I, J, K, GD, iscnv, nsym, ℓ; isplot=true, maxits=10)
    other(i) = 3-i
    T(θ, i) = begin
        z = exp(im*θ)
        z2 = [z^2, z^(-2)]
        den = (I[1]+J[1]+z^(-2)*K[1])*(I[2]+J[2]+z^2*K[2])-I[2]*I[1]
        num = J[other(i)] + z2[i]*K[other(i)]
        num / den
    end

    ω1(θ, i) = begin
        Ti = T(θ, i)
        -imag(GD*Ti) / imag(H*Ti)
    end

    # imag part is machine epsilon by construction
    a²(θ, i) = real((ω1(θ, i)*H + GD) * T(θ, i))
    a²pair(θ) = (a²(θ, 1), a²(θ, 2))

    θroots = Float64[]
    ω1roots = Float64[]
    a²standing = nothing
    if nsym==4abs(ℓ)
        print("n = 4|ℓ|; expecting zeros at mπ/2 ")
        println(iscnv ? "exactly because Cnv" : "shifted because Cn only")
        θs = linspace(-π, π, 100)
        ω1s1 = map(θ->ω1(θ, 1), θs)
        ω1s2 = map(θ->ω1(θ, 2), θs)
        if isplot
            plot(θs, ω1s1, θs, ω1s2)
            xlabel("phase angle")
            ylabel("omega_1")
            title("two expressions for omega1")
        end
        δωs = ω1s2 - ω1s1
        δω(θ) = ω1(θ,2) - ω1(θ,1)
        for i=1:length(δωs)-1
            # if two points straddle, there's a root between
            if δωs[i]*δωs[i+1] < 0
                θ = root(δω, (θs[i], θs[i+1]), maxits=maxits)
                push!(θroots, θ)
                push!(ω1roots, ω1(θ,1))
            end
        end
        println("a² differs across roots because ±i and ±1 solutions")
    elseif iscnv
        println("Cnv but not n=4|ℓ|; zeros everywhere")
        println("|K₊|+|K₋| = ", sum(abs(K)))
        println("if K=0, then T independent of z")
        println("and ω1₊ = ω1₋ satisfied for all z because")
        println("Cnv -> I₊ = I₋)")
        println("outputting random phase angle")
        push!(θroots, rand())
        push!(ω1roots, ω1(θroots[1], 1))
        # stands for all θ allowed
    else
        println("not Cnv and not n=4|ℓ|; zeros nowhere")
        println("|K₊|+|K₋| = ", sum(abs(K)))
        println("ω1₊ = ", ω1(0.0, 1))
        println("ω1₋ = ", ω1(0.0, 2))
        println("^^ if these two not equal at one θ,")
        println("then also non-zero for all θ")
    end
    a²standing = map(a²pair, θroots)
    println("standing lasing mode E = |a₊|E₊ exp(iθ)|a₋|E₋")
    for i=1:length(θroots)
        println("θ", i, " = ", θroots[i]/(π/2), " × π/2")
        println("ω1[", i, "] = ", ω1roots[i])
        println("  |a₊|²[",i,"] = ", a²standing[i][1])
        println("  |a₋|²[",i,"] = ", a²standing[i][2])
    end
    if iscnv
        println("|a₊|²=|a₋|² because mirror symmetry")
    else
        println("|a₊|²!=|a₋|² because not mirror symmetry")
    end
    θroots, ω1roots, a²standing
end

function C_matrix(ɛ::AbstractVector)
    N = length(ɛ)
    C = zeros(5N, 5N)
    complex_block_add!(C, 1, 1, ɛ, scal=-1)
    complex_block_add!(C, 1, 2, ones(N), scal=-1)
    C
end

function B_matrix(ɛ::AbstractVector, ωsalt)
    N = length(ɛ)
    B = zeros(5N, 5N)
    complex_block_add!(B, 1, 1, ɛ, scal=2ωsalt*im)
    complex_block_add!(B, 1, 2, ones(N), scal=2ωsalt*im)
    complex_block_add!(B, 2, 2, ones(N), scal=-im)

    for i=1:N
        B[4N+i, 4N+i] = 1.0
    end
    B
end

function A_matrix(laplacian!::Function, Esalt::AbstractVector, ωsalt::Real, las::Laser, Dlasing::Real, γpar::Real)
    N = length(Esalt)
    A = zeros(5N, 5N)
    laplacian!(sub(A, 1:N, 1:N))
    laplacian!(sub(A, N+1:2N, N+1:2N))

    complex_block_add!(A, 1, 1, las.ɛ, scal = ωsalt^2)
    complex_block_add!(A, 1, 2, ones(N), scal = ωsalt^2)
    complex_block_add!(A, 2, 2, ones(N), scal = las.ωa-ωsalt-im*las.γ⟂)
    for i=1:N
        A[i+4N, i+4N] = γpar
    end

    let x = zeros(N)
        for i=1:N
            x[i] = las.γ⟂*Dlasing*las.F[i] / (1+abs(Esalt[i])^2)
        end
        complex_block_add!(A, 2, 1, x)
    end

    # Esalt satisfies H = 1 / (1 + |Esalt|^2)
    γ = las.γ⟂ / (ωsalt - las.ωa + im*las.γ⟂)
    for i=1:N
        E = Esalt[i] * √γpar / abs(γ)
        P = Dlasing * γ * las.F[i] / (1+abs(Esalt[i])^2) * E
        A[2N+i, 4N+i] = las.γ⟂ * real(E)
        A[3N+i, 4N+i] = las.γ⟂ * imag(E)
        A[4N+i, 2N+i] = imag(E)
        A[4N+i, 3N+i] = -real(E)
        A[4N+i, i]    = -imag(P)
        A[4N+i, N+i]  = real(P)
    end
    Emb = Esalt * √γpar / abs(γ)
    A
end

function complex_block_add!{T<:Real, U<:Number}(A::AbstractMatrix{T}, iblk::Integer, jblk::Integer, x::AbstractVector{U}; scal::Number=1.0)
    N = length(x)
    rows = 2N*(iblk-1)+1:2N*iblk
    cols = 2N*(jblk-1)+1:2N*jblk
    Ablk = sub(A, rows, cols)
    for i=1:N
        z = scal*x[i]
        Ablk[i,i]      +=  real(z)
        Ablk[i+N, i+N] +=  real(z)
        Ablk[i, i+N]   += -imag(z)
        Ablk[i+N, i]   +=  imag(z)
    end
end

function quadeig(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    Z = zeros(A)
    I = eye(A)
    P = [A B; Z -I]
    Q = [Z -C; -I Z]
    Λ, X = eig(P, Q)
    notnan = !isnan(Λ)
    Λ[notnan], X[1:size(A,1), notnan]
end

function smallest_stability_pairs(laplacian!, Esalt, ωsalt, Dlasing, γpar, ɛ)
    A = A_matrix(laplacian!, Esalt, ωsalt, las, Dlasing, γpar)
    B = B_matrix(las.ɛ, ωsalt)
    C = C_matrix(las.ɛ);

    Λ, X = quadeig(A, B, C)
    smallest_ind = sortperm(abs(Λ))[1:4]
    Λ[smallest_ind], X[:, smallest_ind]
end

function smallest_stability_eigs(laplacian!, Esalt, ωsalt, Dlasing, γpar, ɛ)
    Λ, X = smallest_stability_pairs(laplacian!, Esalt, ωsalt, Dlasing, γpar, ɛ)
    Λ
end
