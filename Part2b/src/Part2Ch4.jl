#This file includes all the necessary functions from Chapter 4
#It excludes Runge-Kutta functions and Drivers because they usually change

#Ch 4.1.1 pg 97
#algorithm 41 
#The Fourier Collocation Time Derivative for the Advection-Diffusion Equation

function FourierCollocationTimeDerivative(Phi::Array{Float64}, D::Matrix{Float64}, nu::Float64)

    F = D*Phi
    @. F = nu*F - Phi
    F = D*F #Dphi
    return F


end




#Ch 4.2.1 pg 103 
#alg 44
#Advection Diffusion Time Derivative for Fourier Galerkin
function AdvectionDiffusionTimeDerivative!(Phi::Array{ComplexF64}, DPhi::Array{ComplexF64}, k::Array{Float64}, nu::Float64)

    @. DPhi = -(1.0im * k + nu*k^2) * Phi
    return nothing

end


#alg 46
#Ch 4.2.1 pg 104
#Direct Synthesis of Fourier Galerkin Solution
function EvaluateFourierGalerkinSolution!(Phi::Array{ComplexF64},Phik::Array{ComplexF64},k::Array{Float64}, x::Array{Float64})

    No = length(x)
    phi = zeros(ComplexF64, size(Phik))
    for j in 1:No
        @. phi = Phik*exp(1.0im * k * x[j])
        Phi[j] = sum(phi)
    end

    return nothing

    

end


#alg 49
function FastConvolutionSum!(Q::Matrix{ComplexF64}, V::Matrix{ComplexF64}, W::Matrix{ComplexF64}, N::Int64, planf, planb)

    M = 2*N
    #initialize vectors
    Vkpad = zeros(ComplexF64,M)
    Wkpad = zeros(ComplexF64,M)
    Vx = zeros(ComplexF64,M)
    Wx = zeros(ComplexF64,M)
    Qx = zeros(ComplexF64,M)
    Qkpad = zeros(ComplexF64,M)

    for i in 1:N
        Vx[N + i] = V[i]
        Wx[N + i] = W[i]
    end 

    FFTW.ifftshift!(Vkpad,Vx)
    FFTW.ifftshift!(Wkpad, Wx)

    Vx = planb*Vkpad
    Wx = planb*Wkpad

    @. Qx = Vx*Wx
    Qkpad = planf*Qx
    FFTW.fftshift!(Qx,Qkpad)

    for i in 1:N
        Q[i] = Qx[N + i]
    end 

    return nothing
    

end



#Alg 40b
#Ch 3.5.5 pg 86
# Compuation of Derivative of Chebyshev Transform
function mthFastChebyshevDerivative(m::Int64,f::Array{Float64}, GaussLob::Bool, planf, planb)

    N = length(f)-1
    fk = zeros(Float64, N+1)
    dfk = zeros(Float64, N+1)
    DF = zeros(Float64, N+1)

    if GaussLob
        fk = FastChebyshevTransform2(f, planf)
    elseif !Gausslob
        fk = FastChebyshevTransform1(f, planf)
    end

    dfk = mthChebyshevDerivativeCoefficients(m,fk);
    
    if GaussLob
        DF = InvFastChebyshevTransform2(dfk, planb)
        #DF = InvFastChebyshevTransform2(dfk, planb)
    elseif !Gausslob
        DF = InvFastChebyshevTransform1(dfk, planb)
    end

    return DF    

end


#Alg 40c
#find derivative in Chebyshev spectrum
function mthFastChebyshevDerivativeC(m::Int64,f::Array{Float64}, GaussLob::Bool, planf, planb)

    N = length(f)-1
    #fk = zeros(Float64, N+1)
    dfk = zeros(Float64, N+1)
    #DF = zeros(Float64, N+1)

    dfk = mthChebyshevDerivativeCoefficients(m,f);
    

    return dfk   

end


function LegendreTransformNaive(Phi::Array{Float64}, x::Array{Float64}, w::Array{Float64})

    N = length(Phi)
    Phik = zeros(Float64, N)
    #Lk = zeros(Float64,N)
    #Lkk = zeros(Float64,N)
    ll = 0.0
    den = 0.0

    for k in 1:N
        #@. Lk = 0.0
        #@. Lkk = 0.0
        #calculate legendre functions at x[k]
        den = 0.0
        for j in 1:N
            #Lk[j] = LegendrePolynomial(k-1, x[j])
            ll = LegendrePolynomial(k-1, x[j])
            Phik[k] = Phik[k] + w[j]*Phi[j]*ll
            den = den + ll^2 *w[j]
        end

        Phik[k] = Phik[k]/den
        # @. Lkk = w*Phi*Lk
        # Phik[k] = sum(Lkk)
        # if k<N
        #     Phik[k] = Phik[k]*((k-1)+0.5)
        #     #Phik[k] = Phik[k]*((k-1)+0.5)
        # else
        #     Phik[k] = Phik[k]*(N*0.5)
        # end

    end

    return Phik

end



function InverseLegendreTransformNaive(Phik::Array{Float64}, x::Array{Float64}, w::Array{Float64})

    N = length(Phik)
    Phi = zeros(Float64, N)
    #Lk = zeros(Float64,N)
    #Lkk = zeros(Float64,N)
    ll = 0.0
    kk = 0.0
    den = 0.0

    for j in 1:N
        #@. Lk = 0.0
        #@. Lkk = 0.0
        #calculate legendre functions at x[k]
        den = 0.0
        for k in 1:N
            #Lk[j] = LegendrePolynomial(k-1, x[j])
            ll = LegendrePolynomial(k-1, x[j])
            #kk = 2.0*(k-1)+1
            Phi[j] = Phi[j] + Phik[k]*ll
            #Phi[j] = Phi[j] + w[j]*Phik[k]*ll
            #den = den + ll^2 *w[j]
        end

        #@. Phi = Phi*N
        #Phik[k] = Phik[k]/den
        # @. Lkk = w*Phi*Lk
        # Phik[k] = sum(Lkk)
        # if k<N
        #     Phik[k] = Phik[k]*((k-1)+0.5)
        #     #Phik[k] = Phik[k]*((k-1)+0.5)
        # else
        #     Phik[k] = Phik[k]*(N*0.5)
        # end

    end

    return Phi

end


function FastLegendreTransform(Phi::Array{Float64}, planf, planchebtoleg)

    N = length(Phi)
    phikcheb = zeros(Float64, N)
    phikleg = zeros(Float64, N)

    phikcheb = FastChebyshevTransform2(Phi, planf);
    phikleg = planchebtoleg*phikcheb;

    return phikleg

end



function FastInverseLegendreTransform(Phik::Array{Float64}, planb, planlegtocheb)

    N = length(Phik)
    phikcheb = zeros(Float64, N)
    phi = zeros(Float64, N)

    phikcheb = planlegtocheb*Phik;
    phi = InvFastChebyshevTransform2(phikcheb, planb)#FastChebyshevTransform2(phikcheb, planb);
    #@. phi = phi*N;

    return phi

end


#using only FastTransforms.jl
function FastTransformsLegendre(Phi::Array{Float64}, plancheb, planchebtoleg)

    N = length(Phi)
    phikcheb = zeros(Float64, N)
    phikleg = zeros(Float64, N)

    phikcheb = plancheb*Phi;
    phikleg = planchebtoleg*phikcheb;

    return phikleg

end


function FastInverseTransformsLegendre(Phik::Array{Float64}, planinvcheb, planlegtocheb)

    N = length(Phik)
    phikcheb = zeros(Float64, N)
    phi = zeros(Float64, N)

    phikcheb = planlegtocheb*Phik;
    phi = planinvcheb*phikcheb;

    return phi

end


#alg 52
#ch 4.5.1 pg 127
#Legendre basis modified to vanish at endpoints
function ModifiedLegendreBasis(k::Int64, x::Float64)

    phik = (LegendrePolynomial(k, x) - LegendrePolynomial(k+2, x))/sqrt(4.0*k + 6.0)

    return phik   

end


#alg 53
#ch 4.5.1 pg 127
#Synthesis of Legendre Galerkin Solution at different x from original x
function EvaluateLegendreGalerkinSolution(N::Int64, x::Array{Float64}, Phik::Array{Float64})

    M = length(x)
    Phiout = zeros(Float64, M)
    for i in 1:M
        for k in 0:N-2
            Phiout[i] =Phiout[i] + Phik[k+1] * ModifiedLegendreBasis(k, x[i])
        end
    end

    return Phiout

end


#eqns 4.96 on pg 125
alpha(n::Int64) = 1.0/sqrt(4.0*n + 6.0)
gamma(n::Int64) = -2.0/(2.0*n + 1.0)
mu(n::Int64) = -2.0/(2.0*n + 5.0)
beta(n::Int64) = -(gamma(n) + mu(n))


#alg 54
#ch 4.5.1 pg 128
#Legendre Galerkin Tridiagonal Matrix
function initTMatrix(N::Int64, p::Int64)

    l = zeros(Float64, N) 
    d = zeros(Float64, N+1)
    u = zeros(Float64, N)

    for j in 0:N
        jj = j+1
        d[jj] = beta(2*j+p)*alpha(2*j+p)^2
    end

    for j in 1:N
        l[j] = gamma(2*j+p) * alpha(2*j+p) * alpha(2*(j-1)+p)
        jj = j-1
        u[j] = mu(2*jj+p) * alpha(2*jj+p) * alpha(2*(jj+1)+p)
    end
    
    return l, d, u
end


#alg 54b
#Legendre Galerkin Bidiagonal Matrix
function initTMatrixUpperBidiagonal(N::Int64)

    #l = zeros(Float64, N) 
    # d = zeros(Float64, N)
    # u = zeros(Float64, N-1)
    bidim = zeros(Float64, (N-1, N+1))

    for k in 0:N
        kk = k+1
        for m in 0:(N-2)
            mm = m+1
            if m==k
                bidim[mm,kk] = -alpha(m)*gamma(m)
            elseif m+2==k
                bidim[mm,kk] = mu(m) * alpha(m)
            end

        end
    end

    bidim = sparse(bidim)

    # for j in 0:N-1
    #     jj = j+1
    #     d[jj] = -alpha(2*j+p)*gamma(2*j+p)
    # end

    # for j in 1:N-1
    #     # l[j] = gamma(2*j+p) * alpha(2*j+p) * alpha(2*(j-1)+p)
    #     jj = j-1
    #     u[j] = mu(2*jj+p) * alpha(2*jj+p) 
    # end
    
    return bidim#d, u
end


#alg 55
#ch 4.5.1 pg 128
#Modified Legendre Coefficients from Legendre Coefficients
function ModifiedCoefsFromLegendreCoefs(N::Int64, phik::Array{Float64}, Te::Tridiagonal{Float64, Vector{Float64}}, To::Tridiagonal{Float64, Vector{Float64}})

    #even indices
    Me = floor(Int64, (N-2)/2)
    Mo = floor(Int64, (N-2+1)/2) -1
    #L, D, U = initTMatrix(M, 0)

    #computing RHS into modified basis from Legendre coefficients
    rhse = zeros(Float64, Me+1)
    rhso = zeros(Float64, Mo+1)
    Phik = zeros(Float64, N-1)
    
    for j in 0:Me
        #jj = j+1
        rhse[j+1] = mu(2*j) * alpha(2*j) * phik[(2*j+2)+1] - alpha(2*j) * gamma(2*j) * phik[(2*j)+1]

    end

    #T= Tridiagonal(L, D, U)
    B = Te\rhse;

    for j in 0:Me
        jj = j+1
        Phik[(2*j)+1] = B[jj]
    end

    #odd indices
    #M = floor(Int64, (N-2+1)/2) -1
    #L, D, U = initTMatrix(M, 1)
    
    for j in 0:Mo
        jj = j+1
        rhso[jj] = mu(2*j+1) * alpha(2*j+1) * phik[(2*j+3)+1] - alpha(2*j+1) * gamma(2*j+1) * phik[(2*j+1)+1]

    end

    #T = Tridiagonal(L, D, U)
    B = To\rhso;

    for j in 0:Mo
        jj = j+1
        Phik[(2*j+1)+1] = B[jj]
    end

    return Phik
    
end

#alg 55b
function ModifiedCoefsFromLegendreCoefsMatrix(N::Int64, phik::Array{Float64}, Te::Tridiagonal{Float64, Vector{Float64}}, To::Tridiagonal{Float64, Vector{Float64}})

    #modified PhiK
    Phik = zeros(Float64, N-1)
    #rhs bidiagonal sparse matrix
    bidiagm = initTMatrixUpperBidiagonal(N);
    rhsk = bidiagm*phik
    #even indices
    rhse = @view rhsk[1:2:end]
    # phike = @view phik[1:2:end] #note julia indices start from 1
    # Me = length(phike)
    Me = length(rhse)
    # #odd indices
    rhso = @view rhsk[2:2:end]
    # phiko = @view phik[2:2:end]#note julia indices start from 1
    # Mo = length(phiko)
    Mo = length(rhso)
    #L, D, U = initTMatrix(M, 0)

    #computing RHS into modified basis from Legendre coefficients
    #using bidiagonal matrix
    #db, ub = initTMatrixUpperBidiagonal(Me, 0);
     
    # be = Bidiagonal(db, ub, :U);
    # rhse = be*phike;
    

    #T= Tridiagonal(L, D, U)
    B = Te\rhse;

    for j in 0:Me-1
        jj = j+1
        Phik[(2*j)+1] = B[jj]
    end

    #using bidiagonal matrix
    # db, ub = initTMatrixUpperBidiagonal(Mo, 1);
    # bo = Bidiagonal(db, ub, :U);
    # rhso = bo*phiko;

    #T = Tridiagonal(L, D, U)
    B = To\rhso;

    for j in 0:Mo-1
        jj = j+1
        Phik[(2*j+1)+1] = B[jj]
    end

    return Phik
    
end


#Alg 57
#ch 4.6 pg 133
#Matrix of Legendre Galerkin Approximation
function CGDerivativeMatrix(N::Int64, x::Array{Float64}, w::Array{Float64})

    G = zeros(Float64, (N+1, N+1))
    D = PolynomialDerivativeMatrix(x)
    s = 0.0

    for j in 0:N
        jj = j+1
        for n in 0:N
            nn = n+1
            s = 0.0
            for k in 0:N
                kk = k+1
                s = s+ D[kk,nn]*D[kk,jj]*w[kk]
            end
            G[jj,nn] = -s/w[jj]
        end

    end
    return G
end



#alg 59
#ch 4.7.1 pg 139
#Julia function to create necessary variables
function NodalDiscontinuousGalerkinConstruct(N::Int64, x::Array{Float64}, w::Array{Float64})
    #Inputs
    #N - Number of nodes from 0 to N
    #x - positions of nodes (0 to N) in [-1,1]
    #w - Gauss-Legendre Weights of nodes x
    
    #barycentric weights
    wb = BarycentricWeights(x)
    #interpolation coefficients for l_j(+-1)
    ln1 = LagrangeInterpolatingPolynomials(-1.0, x, wb,1.0e-12)
    lp1 = LagrangeInterpolatingPolynomials(1.0, x, wb,1.0e-12)

    #derivative matrix
    D = PolynomialDerivativeMatrix(x)
    Dhat = zeros(Float64, (N+1,N+1))
    Dhattotal = zeros(Float64, (N+1,N+1))

    for j in 0:N
        jj=j+1
        for i in 0:N
            ii = i+1 
            Dhat[ii,jj] = -D[jj,ii]*w[jj]/w[ii]
            #total Dhat with transpose
            Dhattotal[jj,ii] = -(lp1[ii]*lp1[jj]  - D[ii,jj]*w[ii])/w[jj]
            #total Dhat without transpose
            # Dhattotal[ii,jj] = -(lp1[ii]*lp1[jj] - D[ii,jj]*w[jj]/w[ii])
        end
    end

    return wb, ln1, lp1, D, Dhat, Dhattotal

end

#Alg 60 
#ch 4.7.1 pg 139
#First spatial derivative via Galerkin Approximation
function ComputeDGDerivative(Phi::Array{Float64}, w::Array{Float64}, phiL::Float64, phiR::Float64, Dhat::Array{Float64}, Lp::Array{Float64}, Ln::Array{Float64})

    NN = length(Phi)-1
    dPhi = Dhat*Phi
    
    for j in 0:NN
        jj=j+1
        dPhi[jj] = dPhi[jj] + (phiR * Lp[jj] - phiL * Ln[jj])/w[jj]
        # dPhi[jj] = dPhi[jj] + (phiR * Lp[jj] - phiL * Ln[jj])/w[jj]
    end

    return dPhi
end


#alg 61a
#ch 4.7.1 pg 140
#Interpolation to Boundary for Gauss nodes
function InterpolateToBoundary(Phi::Array{Float64}, Lj::Array{Float64})

    interpolatedValue = 0.0

    for j in 0:N
        jj=j+1
        interpolatedValue = interpolatedValue + Phi[jj]*Lj[jj]
    end

    return interpolatedValue
end
#alg 61
#ch 4.7.1 pg 140
#Time derivative via Discontinuous Galerkin Approximation
function DGTimeDerivative(Phi::Array{Float64},w::Array{Float64}, Dhat::Array{Float64},gt::Float64, Lp::Array{Float64}, Ln::Array{Float64},t::Float64, c::Float64)

    Phil = 0.0
    Phir = 0.0

    if c>0.0
        Phil = gt
        #Phir = InterpolateToBoundary(Phi, Lp)
        Phir = dot(Phi,Lp)
    else
        Phir = gt
        Phil = dot(Phi,Ln)#InterpolateToBoundary(Phi, Ln)
    end
    

    dPhi = ComputeDGDerivative(Phi, w, Phil, Phir, Dhat, Lp, Ln)
    @. dPhi =  -c*dPhi

    return dPhi
    
end







































