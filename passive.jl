PassiveMode(E::AbstractVector, ω::Number) = Mode(E, ω+0.0im, 0.0)

function passive_threshold!{T<:Real}(laplacian!::Function, md::Mode, las::Laser, Ds::Tuple{T, T};
                    init_incr=0.05, tol=1e-10, maxits=10, newtontol=1e-8, newton_maxits=9, isprint=false)
    func(D) = begin
        passive_solve!(laplacian!, md, las, D, tol=newtontol, maxits=newton_maxits)
        imag(md.ω)
    end
    root(func, Ds, tol=tol, maxits=maxits, isprint=isprint)
end    

function passive_update!(md::Mode, dv::AbstractVector)
    md.ω  += dv[end]
    N = length(md.E)
    for i=1:N
        md.E[i] += dv[i]
    end
end

function passive_derivative!(md::Mode, las::Laser, D::Real, J::AbstractMatrix)
    N = length(md.E)
    J[end, md.imax] = 1
    
    ∂f∂ω  = sub(J, 1:N, N+1)
    γ = las.γ⟂ / (md.ω - las.ωa + im*las.γ⟂)
    ω²D∂γ∂ω = md.ω^2*D*(-γ^2/las.γ⟂)
    for i=1:N
        # columns must be set using =
        H = 1
        ∂M∂ω = 2md.ω*(las.ɛ[i]+D*γ*las.F[i]*H) + ω²D∂γ∂ω*las.F[i]*H
        ∂f∂ω[i]   = ∂M∂ω * md.E[i]
    end
end

function passive_salt_operator!(md::Mode, las::Laser, D::Real, J::AbstractMatrix)
    γ = las.γ⟂ / (md.ω - las.ωa + im*las.γ⟂)
    N = length(md.E)    
    for i=1:N
        H = 1
        # limitation: no hole-burning of other lasing mode
        z = md.ω^2 * (las.ɛ[i] + D*γ*H*las.F[i])
        J[i,i] +=  z
    end    
end

function passive_residual!(laplacian!::Function, md::Mode, las::Laser, D::Real, J::AbstractMatrix, f::AbstractVector)
    fill!(J, 0.0)
    fill!(f, 0.0)    
    N = length(md.E)
    laplacian!(sub(J, 1:N, 1:N))
    passive_salt_operator!(md, las, D, J)
    
    A_mul_B!(sub(f, 1:N), sub(J,1:N,1:N), md.E)
    f[end]   = md.E[md.imax] - 1
    norm(f)
end

function passive_solve!(laplacian!::Function, md::Mode, las::Laser, D::Real;
                maxits=9, tol=1e-8, isprint=false)
    N = length(md.E)
    isa(md.ω, Complex) || error("incorrect type of mode")
    length(las.ɛ) == N || error("incorrect length(ɛ)")
    J = zeros(Complex{Float64}, N+1, N+1)
    f = zeros(Complex{Float64}, N+1)
    
    for its=0:maxits
        res = passive_residual!(laplacian!, md, las, D, J, f)
        isprint && println("|f| = ", res)
        res < tol && return (res, its)
        its < maxits || throw(MaxIterationsReached())
        passive_derivative!(md, las, D, J)
        dv = scale!(J\f, -1)
        passive_update!(md, dv)
    end
end
