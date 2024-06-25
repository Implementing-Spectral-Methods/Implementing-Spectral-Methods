module Part1

using FFTW
using Base.Threads
using LinearAlgebra
using SparseArrays
using Test
using TimerOutputs
using HDF5
using Printf
using FastGaussQuadrature
using ApproxFun
using FastTransforms
using FastChebInterp
using Polynomials
using SpecialPolynomials
using BenchmarkTools
using Latexify
using Plots
using LaTeXStrings
using IncompleteLU
using Preconditioners 
using IterativeSolvers
using LinearMaps
using RecursiveFactorization
using LinearOperators
using Krylov
using IJulia


export FourierInterpolantFromModes, FourierInterpolantFromNodes!, LegendreDerivativeCoefficients!, ChebyshevDerivativeCoefficients!, FourierDerivativeByFFT, FourierDerivativeMatrix, LegendrePolynomial, LegendreSum, ChebyshevPolynomialTrig, ChebyshevPolynomialIter, SPLegendre, LegendrePolynomialandDerivative, LegendreGaussNodesAndWeights, qAndLEvaluation, LegendreGaussLobattoNodesAndWeights, ChebyshevGaussNodesAndWeights, ChebyshevGaussLobattoNodesAndWeights, FastChebyshevTransform1, InvFastChebyshevTransform1, FastChebyshevTransform2, InvFastChebyshevTransform2, mthChebyshevDerivativeCoefficients, BarycentricWeights, LagrangeInterpolation, PolynomialInterpolationMatrix, LagrangeInterpolatingPolynomials, CoarseToFineInterpolation2D, LagrangeInterpolantDerivative, PolynomialDerivativeMatrix, mthOrderPolynomialDerivativeMatrix, EOMatrixDerivative, FastChebyshevDerivative
#Ch 1.6.1, Page 18
#algorithm 2
#Direct Evaluation of the Fourier Interpolant from Its Modes
function FourierInterpolantFromModes(FK::Array{ComplexF64}, k::Array{ComplexF64}, x::Float64)
    N = length(k)
    s = 0.5 * ( FK[1] * exp(-1.0im * N * x/2)  + FK[end] * exp(1.0im * N * x/2) )
    for i in 2:N-1
        s = s + FK[i]*exp(1.0im * k[i] * x)
        
    end
    return s

end



#Ch 1.6.1, Page 18
#algorithm 3
#Direct Evaluation of the Fourier Interpolant from Its Nodes
function FourierInterpolantFromNodes!(Fj::Array{ComplexF64}, xj::Array{ComplexF64}, x::Float64, diffprec::Float64)
    N = length(xj)
    t = 0.0
    s = 0.0
    for j in 1:N
        if abs(x-xj[j]) < diffprec
            return Fj[j]
            
        else
            t = (x - xj)/2.0
            s = s + Fj[j]*sin(N*t)*cot(t)/2
        end
        
    end
    return s

end


#Ch 1.10.1, pg 31
#algorithm 4
#Evaluate the Legendre Coefficients of the Derivative of a Polynomial
function LegendreDerivativeCoefficients!(dFk1::Array{ComplexF64},Fk::Array{ComplexF64})
    N = length(Fk1)-1
    dFk1[end] = 0
    dFk1[end-1] = (2*N - 1) * Fk[end]

    for k in N-2:-1:0
        kk = k+1
        dFk1[kk] = (2*k +1) * (Fk[kk+1] + dFk1[kk+2]/(2*k+5))
        
    end

    return nothing

end


#Ch 1.10.1, pg 31
#algorithm 5
#Evaluate the LChebyshev Coefficients of the Derivative of a Polynomial
function ChebyshevDerivativeCoefficients!(dFk1::Array{ComplexF64},Fk::Array{ComplexF64})
    N = length(Fk1)-1
    dFk1[end] = 0
    dFk1[end-1] = (2*N) * Fk[end]

    for k in N-2:-1:1
        kk = k+1
        dFk1[kk] = 2.0*(k +1) *Fk[kk+1] + dFk1[kk+2]
        
    end

    DFk1[1] = Fk[1] + DFk1[2]/2.0

    return nothing

end



#Ch 2.3, pg 54
#algorithm 17
#Fast Evaluation of Fourier Polynomial Derivative
function FourierDerivativeByFFT(f::Array{ComplexF64}, k::Array{Float64}, dn::Int64, planf, plani)

    #fourier transform using FFTW plans
    ftran = planf*f
    #create d/dx derivative with wavenumber
    ikderiv = @. 1.0im*kx^dn
    #multiply wavenumbers by the FT(f)
    ftran .= ikderiv.*ftran
    #below is the part from the book that is unnecessary with Julia
    # if mod(dn,2.0) != 0
    #     ftran[end] = 0.0 + 0.0im
    # end
    #inverse fourier transform
    g = plani * ftran
    # @. g = g/(length(g)^2)
    
    return g

end


#ch 2.4, pg 55
#algorithm 18
#Computation of Fourier Derivative Matrix using negative sum trick
function FourierDerivativeMatrix(N::Int64)

    D = zeros(Float64,(N,N))
    ii= 0
    jj = 0
    for i in 0:(N-1)
        ii = i+1
        D[ii,ii] = 0.0
        for j in 0:(N-1)
            jj = j+1
            if j!=i
                D[ii,jj] = 0.5 * (-1.0)^(i+j) * cot((i-j)*pi/N)
                D[ii,ii] = D[ii,ii] - D[ii,jj]
            end

        end

    end
    return D

end



#ch 3.1 pg 60
#Algorithm 20
#evaluate Legendre Polynomials of Degree k using three term recursion
function LegendrePolynomial(k::Int64, x::Float64)
    if k==0
        return 1.0

    elseif k==1
        return x
    end
    Lk = 0.0
    Lkm2 = 1 # L_{k-2}
    Lkm1 = x #L_{k-1}
    for j in 2:k
        Lk = (2.0*j - 1)/j *x *Lkm1 - (j-1)/j * Lkm2
        Lkm2 = Lkm1
        Lkm1 = Lk
    end

    return Lk
    
end


#sum of Legendre series polynimoals at point x
#
function LegendreSum(k::Int64, x::Float64)

    if k==0
    return 1.0

    elseif k==1
        return x
    end
    Lk = 0.0
    
    Lkm2 = 1 # L_{k-2}
    Lkm1 = x #L_{k-1}
    Lksum = Lkm2 + Lkm1
    for j in 2:k
        Lk = (2.0*j - 1)/j *x *Lkm1 - (j-1)/j * Lkm2
        Lkm2 = Lkm1
        Lkm1 = Lk
        Lksum = Lksum + Lk
    end

    return Lksum
end


# ch 3.1, page 60
#Algorithm 21
#evaluate Chebyshev Polynomials of Degree k using three term recursion and trigonometric forms
function ChebyshevPolynomialTrig(k::Int64, x::Float64)
    if k==0
        return 1.0

    elseif k==1
        return x
    end

    return cos(k*acos(x))
    
end


function ChebyshevPolynomialIter(k::Int64, x::Float64, Ks::Int64)
    if k==0
        return 1.0

    elseif k==1
        return x
    end
    Tk = 0.0
    Tkm2 = 1 # T_{k-2}
    Tkm1 = x #T_{k-1}
    for j in 2:k
        Tk = 2.0*x *Tkm1 - Tkm2
        Tkm2 = Tkm1
        Tkm1 = Tk
    end

    return Tk
    
end

function SPLegendre(j::Int64,x::Float64)
    oo = ones(Float64, j+1);
    Ls = SpecialPolynomials.Legendre(oo)
    return Ls(x)
    

end


#Ch 3.2.1 pg 63
#Algorithm 22
#The Legendre Polynomial of degree k and its Derivative using Three Term Recursion
function LegendrePolynomialandDerivative(N::Int64, x::Float64)
    Ln =0.0
    Lnm1 = 0.0
    Lnm2 = 0.0
    dLn =0.0
    dLnm1 = 0.0
    dLnm2 = 0.0
    if N==0
        Ln = 1
        dLn = 0
    elseif N==1
        Ln = x
        dLn = 1
    else
        Lnm2 = 1
        Lnm1 = x
        dLnm2 = 0
        dLnm1 = 1

        for j in 2:N
            Ln = (2.0*j - 1)/j *x *Lnm1 - (j-1)/j * Lnm2
            dLn = dLnm2 + (2.0*j-1.0)*Lnm1
            Lnm2 = Lnm1
            Lnm1 = Ln
            dLnm2 = dLnm1
            dLnm1 = dLn
            #Lnsum = Lksum + Lk

        end
        
    end

    return Ln, dLn

end


#Ch 3.2.1 pg 64
#algorithm 23
function LegendreGaussNodesAndWeights(N::Int64, nit::Int64, TOL::Float64)

    x = zeros(Float64,N+1)
    w = zeros(Float64,N+1)

    if N==0
        x[1] = 0.0
        w[1] = 2.0
    elseif N ==1
        x[1] = -1.0*sqrt(1/3)
        w[1] = 1.0
        x[2] = -x0
        w[2] = w0
    else
        for j in 0:(floor(Int64,((N+1)/2))-1)
            x[j+1] = - cos((2*j+1)/(2*N+2)*pi)
            for k in 0:nit
                Lnp1, dLnp1 = LegendrePolynomialandDerivative(N+1, x[j+1])
                delta = -Lnp1/dLnp1
                x[j+1] = x[j+1] + delta
                if abs(delta) < TOL*abs(x[j+1])
                    break
                end
            end
            Lnp1, dLnp1 = LegendrePolynomialandDerivative(N+1, x[j+1])
            x[end-j] =-x[j+1]
            w[j+1] = 2.0/(1.0 - x[j+1]^2) / dLnp1^2
            w[end-j] = w[j+1]

        end
        
    end
    if mod(N,2) == 0
        Lnp1, dLnp1 = LegendrePolynomialandDerivative(N+1, 0.0)
        x[ceil(Int64,N/2)+1] = 0.0
        w[ceil(Int64,N/2)+1] = 2.0 / dLnp1^2
        
    end

    return x, w


end


# Ch 3.2.2 pg 65
#algorithm 24
#Combined algorithm to compute Ln(x), q(x) = L_{N+1}(x) - L_{N-1}, and q'(x)
function qAndLEvaluation(N::Int64, x::Float64)

    k=2
    Lnm2 = 1.0
    Lnm1 = x
    Ln =0.0
    dLn = 0.0
    dLnm2 = 0.0
    dLnm1 = 1.0


    for j in 2:N
        Ln = (2.0*j - 1)/j *x *Lnm1 - (j-1)/j * Lnm2
        dLn = dLnm2 + (2.0*j-1.0)*Lnm1
        Lnm2 = Lnm1
        Lnm1 = Ln
        dLnm2 = dLnm1
        dLnm1 = dLn
        #Lnsum = Lksum + Lk

    end

    q = dLn
    #second derivative from https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas
    dq = (2*x*dLn - N*(N+1)*Ln)/(1-x^2) #ddLn
    
    # k = N+1
    # Lnp1 = (2*k-1)/k * x * Ln - (k-1)/k * Lnm1
    # dLnp1 = dLnm1 + (2*k-1)*Ln

    # q = Lnp1 - Lnm1
    # dq = dLnp1 - dLnm1
    
        
    

    return q, dq, Ln



end


# Ch 3.2.2 pg 66
#algorithm 25
#calculates the gauss lobatto nodes and weights
function LegendreGaussLobattoNodesAndWeights(N::Int64, nit::Int64, TOL::Float64)

    x = zeros(Float64,N+1)
    w = zeros(Float64,N+1)
    xold = 0.0

    if N ==1
        x[1] = -1.0
        w[1] = 1.0
        x[2] = 1.0
        w[2] = w[1]
    else
        x[1] = -1.0
        w[1] = 2.0/(N*(N+1))
        x[end] = 1.0
        w[end] = w[1]
        #xold = 0.0
        
        for j in 1:(floor(Int64,((N+1)/2))-1)
            x[j+1] = - cos((j+0.25)/N*pi - 3.0/(8.0*N*pi) /(j + 0.25))
            #println(j)
            #x[j+1] = -cos(pi*j/N)
            for k in 0:nit
                q, dq, Ln = qAndLEvaluation(N, x[j+1])
                delta = -q/dq
                xold = x[j+1]
                x[j+1] = x[j+1] + delta
                if abs(x[j+1] - xold) < TOL #abs(delta) < TOL*abs(x[j])
                    break
                end
                if k==nit-1
                    println("Reached Iterator Max")
                end
            end
            q, dq, Ln = qAndLEvaluation(N, x[j+1])
            x[end-j] =-x[j+1]
            w[j+1] = 2.0/(N*(N+1)) / Ln^2
            w[end-j] = w[j+1]

        end
        
    end
    if mod(N,2) == 0
        q, dq, Ln = qAndLEvaluation(N, 0.0)
        x[ceil(Int64,N/2)+1] = 0.0
        w[ceil(Int64,N/2)+1] = 2.0/(N*(N+1)) / Ln^2
        
    end

    return x, w



end


# Ch 3.2.3 pg 67
#algorithm 26
#Calculates Chebyshev Guass Nodes and Weights
function ChebyshevGaussNodesAndWeights(N::Int64)
    x = zeros(Float64,N+1)
    w = zeros(Float64,N+1)

    for j in 0:N
        x[j+1] = -cos((2*j+1)*pi/(2*N+2))
        w[j+1] = pi/(N+1)

    end

    return x, w

end


# Ch 3.2.3 pg 68
#algorithm 27
#Calculates Chebyshev Guass Nodes and Weights
function  ChebyshevGaussLobattoNodesAndWeights(N::Int64)

    x = zeros(Float64,N+1)
    w = zeros(Float64,N+1)

    for j in 0:N
        x[j+1] = -cos(pi*j/N)
        w[j+1] = pi/(N)

    end
    w[1] = w[1]/2.0
    w[end] = w[end]/2.0
    
    return x, w

end


#gauss points
function FastChebyshevTransform1(g::Array{Float64}, plan)
    n = length(g)
    ftran = plan*g #pr*g;
    @. ftran = sqrt(2.0/n)*ftran;
    ftran[1] = ftran[1]/sqrt(2.0)
    return ftran

end



function InvFastChebyshevTransform1(g::Array{Float64}, plan)
    n = length(g)
    ftran = plan*g #pr*g;
    #@. ftran = sqrt(2.0/n)*ftran;
    ftran[1] = ftran[1]*sqrt(2.0)
    return ftran

end


#gauss lobatto points
function FastChebyshevTransform2(g::Array{Float64}, plan)
    #println("FCT2")
    n = length(g)
    nscale = 1.0/(n-1)#1/sqrt(2*(n-1)) /(n-1)  #
    ft = deepcopy(g)
    # ft[1] = g[1]/2.0
    # ft[end] = g[end]/2.0
    ftran  = plan*ft 
    @. ftran = nscale*ftran
    # @. ftran = 1.0/(n-1)* ftran;
    ftran[1] = ftran[1]/2.0#*sqrt(0.5)
    ftran[end] = ftran[end]/2.0#*sqrt(0.5)
    return ftran
end


function InvFastChebyshevTransform2(g::Array{Float64}, plan)
    #println("IFCT2")
    n = length(g)
    nscale = 1.0#/(2*(n-1))#1.0/sqrt(2*(n-1))
    ft = deepcopy(g)
    @. ft[2:end-1] = ft[2:end-1]/2.0
    #ft[1] = g[1]*2.0
    #ft[end] = g[end]*2.0
    ftran  = plan*ft 
    @. ftran = nscale*ftran
    # @. ftran = 1.0/(n-1)* ftran;
    #ftran[1] = ftran[1]*2.0
    #ftran[end] = ftran[end]*2.0
    return ftran

end


#algorithm 5b
#this algorithm takes algorithm 5 and changes to give the mth derivative of the Chebyshev coefficients
function mthChebyshevDerivativeCoefficients(m::Int64, Fk::Array{Float64})
    N = length(Fk)-1
    #dFk = zeros(Float64, N+1)
    #@. dFk = Fk
    dFk = deepcopy(Fk)
    dFkm = zeros(Float64,N+1)

    for l in 1:m
        @. dFkm = 0.0
        dFkm[end] = 0.0
        dFkm[end-1] = 2.0*N * dFk[end]
        # dFkm[end-l+1] = 0
        # dFkm[end-l] = (2*(N-l+1)) * dFk[end-l+1]
    
        #for k in N-2:-1:1
        for k in N-2:-1:1
            kk = k+1
            dFkm[kk] = 2.0*(k +1) *dFk[kk+1] + dFkm[kk+2]
            
        end
    
        dFkm[1] = dFk[2] + dFkm[3]/2.0
        @. dFk = dFkm
        
    end

    return dFkm

end


#ch 3.4 pg 75
#Algorithm 30
#Weights for Lagrange Interpolation
function BarycentricWeights(x::Array{Float64})

    n = length(x)-1
    w = ones(length(x))

    for j in 1:n
        jj=j+1
        for k in 0:j-1
            kk=k+1
            w[kk] = w[kk]*(x[kk] - x[jj])
            w[jj] = w[jj]*(x[jj] - x[kk])
            
        end
    end

    @. w = 1.0 / w

    return w


end


#ch 3.4 pg 75
#Algorithm 31 
#Lagrange Interpolant from Barycentric Form
function LagrangeInterpolation(x::Float64, xx::Array{Float64}, f::Array{Float64}, w::Array{Float64},xtol::Float64)

    numerator = zeros(Float64, length(f))
    denominator = zeros(Float64, length(f))
    t = 0.0

    for j in 1:N
        if abs(x-xx[j]) < xtol
            return f[j]
        end
        t = w[j]/(x-xx[j])
        numerator = numerator + t * f[j]
        denominator = denominator + t

    end

    return numerator/denominator

end

#ch 3.4 pg 76
#Algorithm 32
#matrix for interpolation between two sets of points
function PolynomialInterpolationMatrix(x::Array{Float64}, xx::Array{Float64}, w::Array{Float64},xtol::Float64)

    T = zeros((length(xx), length(x)))
    rowhasmatch = false
    M = length(xx)-1
    N = length(x)-1
    s = 0.0
    t = 0.0
    for k in 0:M
        kk = k+1
        rowhasmatch = false
        for j in 0:N
            jj = j+1
            T[kk,jj] = 0.0
            if abs(x[jj] - xx[kk]) < xtol
                rowhasmatch = true
                T[kk,jj] = 1.0
            end
        end
        if !rowhasmatch
            s = 0.0
            for j in 0:N
                jj=j+1
                t = w[jj]/(xx[kk] - x[jj])
                T[kk,jj] = t
                s = s+t

            end
            for j in 0:N
                jj=j+1
                T[kk,jj] = T[kk,jj]/s
            end
        end
    end

    return T

end


#ch 3.4 pg 77
#Algorithm 34
# 
function LagrangeInterpolatingPolynomials(x::Float64, xx::Array{Float64}, w::Array{Float64},xtol::Float64)
    N = length(xx)-1
    lag = zeros(Float64, N+1)
    xMatchesNode = false
    s = 0.0
    t=0.0
    for j in 0:N
        jj=j+1
        lag[jj] = 0.0
        if abs(x - xx[jj]) < xtol
            lag[jj] = 1.0
            xMatchesNode = true
        end
        
    end
    if xMatchesNode
        return lag
    end
    s = 0.0
    for j in 0:N
        jj=j+1
        t = w[jj]/(x-xx[jj])
        lag[jj] = t
        s = s + t
    end

    @. lag = lag/s

    return lag

end


#Ch 3.5 pg 79
#Algorithm 35
#interpolation from Coarse to Fine grid in 2D
function CoarseToFineInterpolation2D(x::Array{Float64},  y::Array{Float64}, f::Matrix{Float64}, xx::Array{Float64},  yy::Array{Float64},xtol::Float64)

    w = BarycentricWeights(x)
    #Tnewold = PolynomialInterpolationMatrix(x, xx, w,xtol)
    T = PolynomialInterpolationMatrix(x, xx, w,xtol)
    Nold = length(x)
    Mold = length(y)
    Nnew = length(xx)
    Mnew = length(yy)

    #for j in 1:Mold
    #end
    #just need matrix multiplication by vector
    #Fnewold = [Tnewold * f[:,j] for j in 1:Mold]
    F = [T * f[:,j] for j in 1:Mold]

    w = BarycentricWeights(xx)
    T = PolynomialInterpolationMatrix(y, yy, w,xtol)
    F = [T * F[j,:] for j in 1:Nnew]

    return F
    

end


#Ch 3.5.1 pg 80
#Algorithm 36
#Direct Computation of the Polynomial derivative in Barycentric Form
function LagrangeInterpolantDerivative(x::Float64, xx::Array{Float64}, f::Array{Float64}, w::Array{Float64},xtol::Float64)

    numerator = 0.0
    denominator = 0.0
    atNode = false
    i = zero(Int64)
    p = 0.0
    t = 0.0
    
    N = length(x)
    for j in 1:N
        if abs(x - xx[j]) < xtol
            atNode = true
            p = f[j]
            denominator = -w[j]
            i = j
            
        end
    end

    if atNode
        for j in 1:N
            if j != i
                numerator = numerator +w[j] *(p-f[j])/(x - xx[j])
            end

        end
    else
        denominator = 0.0
        p = LagrangeInterpolation(x, xx, f, w,xtol)
        for j in 1:N
            t = w[j]/(x - xx[j])
            numerator = numerator + t*(p-f[j])/(x - xx[j])
            denominator  = denominator +t
        end

    end

    return numerator/denominator
    

end


#Ch 3.5.2 pg 82
#Algorithm 37
#First derivative approximation matrix
function PolynomialDerivativeMatrix(x::Array{Float64})

    N = length(x)-1
    w = BarycentricWeights(x)
    D = zeros(Float64, (N+1,N+1))

    for i in 0:N
        ii = i+1
        D[ii,ii] = 0.0
        for j in 0:N
            jj=j+1
            if j != i
                D[ii,jj] = w[jj]/w[ii] /(x[ii] - x[jj])
                D[ii,ii] = D[ii,ii] - D[ii,jj] #+ 1/(x[jj] - x[ii]) #
            end
        end

    end
    return D

end


#Algorithm 38
function mthOrderPolynomialDerivativeMatrix(m::Int64, x::Array{Float64})
    N = length(x)-1
    w = BarycentricWeights(x)
    #D = zeros(Float64, (N,N))

    D = PolynomialDerivativeMatrix(x)

    if m == 1
        return D
    end

    Dm = zeros(Float64, size(D))
    
    for k in 2:m
        for i in 0:N
            ii=i+1
            Dm[ii,ii] = 0.0
            for j in 0:N
                jj = j+1
                if j != i
                    Dm[ii,jj] = k /(x[ii] - x[jj]) *  (w[jj]/w[ii] * D[ii,ii] - D[ii,jj] )
                    Dm[ii,ii] = Dm[ii,ii] - Dm[ii,jj] # this might be Dm[ii,ii] - Dm[ii,jj]
                end
            end
        end
        @. D = Dm
    end

    return Dm
end


#alg 39
function EOMatrixDerivative(f::Array{Float64}, D::Matrix{Float64})

    N = length(f)
    Df = zeros(Float64, N)


    #M = floor(Int64,(N+1)/2)
    #N will be to total number of points including 0th point
    M = floor(Int64,N/2)
    MM = M+1 #because Julia starts from 1
    ej = zeros(Float64, M)
    oj = zeros(Float64, M)
    De = zeros(Float64, M)
    Do = zeros(Float64, M)

    jj = zero(Int64)
    ii = zero(Int64)

    for j in 0:M
        jj = j+1
        ej[jj] =  0.5 * (f[jj] + f[N -j])
        oj[jj] =  0.5 * (f[jj] - f[N -j])
    end
    for i in 0:M-1
        ii = i+1
        De[ii] = 0.0
        Do[ii] =  0.0
        for j in 0:M-1
            jj = j+1
            De[ii] = De[ii] + (D[ii,jj] + D[ii,N-j])*ej[jj]
            Do[ii] = Do[ii] + (D[ii,jj] - D[ii,N-j])*oj[jj]
        end
    end

    if mod(N,2) != 0
        for i in 0:M-1
            ii = i+1
            De[ii] = De[ii] + D[ii,M] *ej[M]
        end
        De[M+1] = 0.0
        for j in 0:M-1
            jj = j+1
            Do[M] = Do[M] + (D[M,jj] - D[M,N-j]) *oj[jj]
        end
    end

    @views @. Df[1:M] = De + Do
    @views @. Df[M+1:end] = -De[end:-1:1] + Do[end:-1:1]

    if mod(N,2) != 0
        Df[M] = De[M] + Do[M]
    end

    return Df

end


#Alg 40
function FastChebyshevDerivative(f::Array{Float64}, GaussLob::Bool, planf, planb)

    N = length(f)-1
    fk = zeros(Float64, N+1)
    dfk = zeros(Float64, N+1)
    DF = zeros(Float64, N+1)

    if GaussLob
        fk = FastChebyshevTransform2(f, planf)
    elseif !Gausslob
        fk = FastChebyshevTransform1(f, planf)
    end

    ChebyshevDerivativeCoefficients!(dfk,fk);
    if GaussLob
        DF = InvFastChebyshevTransform2(dfk, planb)
    elseif !Gausslob
        DF = InvastChebyshevTransform1(dfk, planb)
    end

    return DF    

end












































end # module Part1
