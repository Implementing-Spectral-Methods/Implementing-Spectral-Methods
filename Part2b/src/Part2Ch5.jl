#Algorithm 37b
function PolynomialDerivativeMatrixElement(x::Array{Float64}, w::Array{Float64},II::Int64,JJ::Int64)
    #need barycentric weights
    N = length(x)-1
    Dij = 0.0
    #w = BarycentricWeights(x)
    # D = zeros(Float64, (N+1,N+1))
    if II!=JJ 
        Dij = w[JJ+1]/w[II+1] /(x[II+1] - x[JJ+1])

    elseif II==JJ
        # for i in 0:N
        #     ii = i+1
        #     Dij = 0.0
        for j in 0:N
            jj=j+1
            if j != II
                # D[ii,jj] = w[jj]/w[ii] /(x[ii] - x[jj])
                # D[ii,ii] = D[ii,ii] - D[ii,jj] #+ 1/(x[jj] - x[ii]) #
                Dij = Dij -  w[jj]/w[II+1] /(x[II+1] - x[jj])
            end
        end
    
        # end
    end
    return Dij

end


#algorithm 38b
function mthPolynomialDerivativeMatrixElement(m::Int64,x::Array{Float64},w::Array{Float64}, II::Int64,JJ::Int64)
    #need barycentric weights
    N = length(x)-1
    #Dij = 0.0
    #w = BarycentricWeights(x)
    # D = zeros(Float64, (N+1,N+1))
    mDij = [ PolynomialDerivativeMatrixElement(x, w,II,jjj) for jjj in 0:N]
    #Dii = [ PolynomialDerivativeMatrixElement(x, w,jjj,jjj) for jjj in 0:N]
    Dij = zeros(Float64,N+1)
    #mDii = zeros(Float64,N+1)

    if m==1
        return Dij
    else
        for mm in 2:m
            @. Dij = mDij #copy at start of loop so it will end when derivative is calculated rather than copy
            for j in 0:N
                jj=j+1
                if jj!=II+1 #exclude when i=j
                    mDij[jj] = mm/(x[II+1]-x[jj])*(w[jj]/w[II+1]*Dij[II+1] -Dij[jj])
                end
            end
            #calculate diagonal
            mDij[II+1] = -1.0 * sum(mDij[1:end .!=(II+1)])
            #
        end
        
    end
    
    
    return mDij[JJ+1]

end


#alg 64
#ch 5.2.1.1 pg 155
# Function to construct necessary variables for 2D 
function NodalPotentialConstruct(x::Array{Float64},y::Array{Float64}, ww::Array{Float64})

    #uses Algorithm 38 (mthOrderPolynomialDerivativeMatrix)
    #inputs
    #xy - 2D lobatto points of grid, first index  is x-coordinates, second is y-coord
    #ww - 2D weights of lobatto points

    # N = lengt(x[:,1])-1
    # M = lengt(y[1,:])-1
    N = length(x)-1
    M = length(y)-1
    
    # Dx = mthOrderPolynomialDerivativeMatrix(1, xy[:,1])
    # DDx = mthOrderPolynomialDerivativeMatrix(2, xy[:,1])
    # Dy = mthOrderPolynomialDerivativeMatrix(1, xy[1,:])
    # DDy = mthOrderPolynomialDerivativeMatrix(1, xy[1,:])
    Dx = mthOrderPolynomialDerivativeMatrix(1, x)
    DDx = mthOrderPolynomialDerivativeMatrix(2, x)
    Dy = mthOrderPolynomialDerivativeMatrix(1, y)
    DDy = mthOrderPolynomialDerivativeMatrix(1, y)

    return Dx, DDx, Dy, DDy

end


#Alg 66
#ch 5.2.1.1 pg 156
# Anonymous Function for Collocation Laplace Operator
#will be initialized after inputs x,y,wbx, wby are created
#Note that this operators works on a 1D Vector so any 2D vector must be flattened to 1D
CollocationLapOp = (x::Array{Float64}, y::Array{Float64}, wbx::Array{Float64}, wby::Array{Float64}) -> LinearMap((length(x))*(length(y)); ismutating=true) do dPhi,Phi

    #uses barycentric weights
    N = length(x)-1
    M = length(y)-1

    #set all to zeros
    @. dPhi = 0.0
    #println(length(dPhi))

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii=i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            for k in 0:N #sum over x
                kk = k+1
                nx = kk + (jj-1)*(N+1)
                dPhi[nn] = dPhi[nn] + mthPolynomialDerivativeMatrixElement(2,x,wbx, i,k)*Phi[nx]
            end
            for k in 0:M #sum over y
                kk = k+1
                ny = ii + (kk-1)*(N+1)
                dPhi[nn] = dPhi[nn] + mthPolynomialDerivativeMatrixElement(2,y,wby, j,k)*Phi[ny]
                #dPhi[nn] = dPhi[nn] + PolynomialDerivativeMatrixElement(y, j,k)*PolynomialDerivativeMatrixElement(y, j,k)*Phi[ny]
            end

        end
    end


    return dPhi


end


#alg 69
#ch 5.2.1.3 pg 160
#RHS construction of equation
#This function creates a 1D Vector of the RHS from a 2D input
function CollocationRHSComputation!(N::Int64, M::Int64, Phi::Array{Float64},RHS::Array{Float64},S::Array{Float64}, Dx::Array{Float64},Dy::Array{Float64})

    #inputs
    #N,M - number of nodes minus one for each dimension (x,y)
    #Phi - current potential (2D)
    #RHS - Flattened vector of size L = (N-1)*(M-1)
    #S - 2D steady state source function
    #Dx,Dy - polynomial derivative matrix in x,y direction
    
    L = (N-1)*(M-1)
    #RHS = zeros(Float64, L+1)
    jj = zero(Int64)
    ii = zero(Int64)
    nn = zero(Int64)

    for j in 1:(M-1)
        jj=j+1
        for i in 1:(N-1)
            ii = i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            RHS[nn] = S[ii,jj] - Dx[ii, 0+1]*Phi[0+1,jj] - Dx[ii,N+1]*Phi[N+1,jj] - Dy[jj,0+1]*Phi[ii,0+1] - Dy[jj,M+1]*Phi[ii,M+1]

        end
    end

    return nothing

end


#alg 70 
#ch 5.2.1.3 pg 161
#Direct Matrix Construction of Laplace Operator
function LaplaceCollocationMatrix!(N::Int64, M::Int64, A::Array{Float64}, Dx::Array{Float64},Dy::Array{Float64})
    #inputs
    #N,M - number of nodes minus one for each dimension (x,y)
    #A - LHS matrix (2D) for A*Phi = RHS with size LxL
    #Dx,Dy - polynomial derivative matrix in x,y direction

    L = (N-1)*(M-1)
    @. A = 0.0
    nn = zero(Int64)
    mm = zero(Int64)

    for j in 1:M-1
        jj = j+1
        for i in 1:N-1
            ii = i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            for k in 0:N#1:N-1
                kk = k+1
                mm = kk + (jj-1)*(N+1) #columnwise/fortran
                # mm = j + (i-1)*(M-1) #rowwise/c
                A[nn,mm] = Dx[ii,kk]
            end
            for k in 0:M#1:M-1
                kk=k+1
                mm = ii + (kk-1)*(N+1) #columnwise/fortran
                # mm = j + (i-1)*(M-1) #rowwise/c
                A[nn,mm] = A[nn,mm] + Dy[jj,kk]
            end

        end
    end
    #boundary conditions
    @. A[1,:] = 0.0
    @. A[end,:] = 0.0
    @. A[:,1] = 0.0
    @. A[:,end] = 0.0

    return nothing

    
end


#eqn 5.37 & 5.38
#pg 163
Aij(dxi::Float64, dxip1::Float64,dyj::Float64, dyjp1::Float64) = -2*(1/(dxi*dxip1) + 1/(dyj*dyjp1))
Bij(dxi::Float64, dxip1::Float64) = 2*1/(dxi*(dxi+dxip1))
Cij(dyj::Float64, dyjp1::Float64) = 2*1/(dyj*(dyj+dyjp1))
Eij(dxi::Float64, dxip1::Float64) = 2*1/(dxip1*(dxi+dxip1))
Fij(dyj::Float64, dyjp1::Float64) = 2*1/(dyjp1*(dyj+dyjp1))

#Alg 72 & Alg 73
#Ch 5.2.1.5 pg 166
#Anonymous Function for Finite Difference Preconditioner
FDPreconditionerMul = (x::Array{Float64}, y::Array{Float64}) -> LinearMap((length(x))*(length(y)); ismutating=true) do Z,U

    N = length(x) - 1
    M = length(y) - 1
    @. Z = 0.0

    for j in 1:M-1
        jj = j+1
        for i in 1:N-1
            ii = i+1
            nn = ii + (jj-1)*(N+1)
            nnim1 = ii-1 + (jj-1)*(N+1)
            nnjm1 = ii + (jj-1-1)*(N+1)
            nnip1 = ii+1 + (jj-1)*(N+1)
            nnjp1 = ii + (jj-1+1)*(N+1)
            
            Z[nn] = Aij( (x[ii]-x[ii-1]) , (x[ii+1]-x[ii]), (y[jj]-y[jj-1]), (y[jj+1]-y[jj]) ) * U[nn] + Bij( (x[ii]-x[ii-1]) , (x[ii+1]-x[ii]) ) * U[nnim1] +  Cij( (y[jj]-y[jj-1]), (y[jj+1]-y[jj]) ) * U[nnjm1] + Eij( (x[ii]-x[ii-1]) , (x[ii+1]-x[ii])  ) * U[nnip1] + Fij( (y[jj]-y[jj-1]), (y[jj+1]-y[jj])  ) * U[nnjp1]

        end
    end

    return Z
    

end


#eqn 5.74
#ch 5.2.2 pg 177
#This function produces one elment at the correct indices
function CNGLaplacianElement(x::Array{Float64},w::Array{Float64},wb::Array{Float64}, II::Int64,NN::Int64)

    N = length(x)-1

    G = 0.0

    for k in 0:N
        kk=k+1
        G = G + PolynomialDerivativeMatrixElement(x, wb,k,NN)*PolynomialDerivativeMatrixElement(x, wb,k,II)*w[kk]

    end

    return G   
        
end


#eqn 5.74
#ch 5.2.2 pg 177
#This function produces the entire matrix
function CNGDerivativeMatrix(N::Int64, x::Array{Float64}, w::Array{Float64})

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
            G[jj,nn] = s
        end

    end
    return G
end


#alg 69a
function NodalGalerkinRHSComputation!(N::Int64, M::Int64, Phi::Array{Float64},x::Array{Float64}, y::Array{Float64},RHS::Array{Float64},S::Array{Float64}, wx::Array{Float64}, wy::Array{Float64},wbx::Array{Float64},wby::Array{Float64})


    #inputs
    #N,M - number of nodes minus one for each dimension (x,y)
    #Phi - current potential (2D)
    #RHS - Flattened vector of size L = (N-1)*(M-1)
    #S - 2D steady state source function
    #Dx,Dy - polynomial derivative matrix in x,y direction
    
    
    L = (N-1)*(M-1)
    #RHS = zeros(Float64, L+1)
    jj = zero(Int64)
    ii = zero(Int64)
    nn = zero(Int64)
    @. RHS = 0.0

    for j in 1:(M-1)
        jj=j+1
        for i in 1:(N-1)
            ii = i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            RHS[nn] = -wx[ii]*wy[jj]*S[ii,jj] - wy[jj]*CNGLaplacianElement(x,wx,wbx, i,0)*Phi[0+1,jj] - wy[jj]*CNGLaplacianElement(x,wx,wbx, i,N)*Phi[N+1,jj] - wx[ii]*CNGLaplacianElement(y,wy,wby, j,0)*Phi[ii,0+1] - wx[ii]*CNGLaplacianElement(y,wy,wby, j,M)*Phi[ii,M+1]
            #RHS[nn] = S[ii,jj] - Dx[ii, 0+1]*Phi[0+1,jj] - Dx[ii,N+1]*Phi[N+1,jj] - Dy[jj,0+1]*Phi[ii,0+1] - Dy[jj,M+1]*Phi[ii,M+1]

        end
    end
 
    return nothing

end


#Alg 66
#eqn 5.77
#and eqn 5.79

CNGLapOperatorLinMaps = (x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64},wbx::Array{Float64}, wby::Array{Float64}) -> LinearMap((length(x))*(length(y)); ismutating=true) do dPhi,Phi

    #uses barycentric weights
    N = length(x)-1
    M = length(y)-1

    #set all to zeros
    @. dPhi = 0.0
    #println(length(dPhi))

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii=i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            #dPhi[nn] = 0.0 #dPhi[ii,jj] = 0.0
            for k in 1:(N-1) #sum over x
                kk = k+1
                nx = kk + (jj-1)*(N+1) #for Phi[k,j]
                dPhi[nn] = dPhi[nn] - wy[jj]*CNGLaplacianElement(x,wx,wbx, i,k)*Phi[nx]
            end
            for k in 1:(M-1) #sum over y
                kk = k+1
                ny = ii + (kk-1)*(N+1)
                dPhi[nn] = dPhi[nn] - wx[ii]*CNGLaplacianElement(y,wy,wby, j,k)*Phi[ny]
                #dPhi[nn] = dPhi[nn] + PolynomialDerivativeMatrixElement(y, j,k)*PolynomialDerivativeMatrixElement(y, j,k)*Phi[ny]
            end

        end
    end


    return dPhi


end


#Alg 66
#eqn 5.77
#and eqn 5.79

CNGLapOperatorLinOp = (x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64},wbx::Array{Float64}, wby::Array{Float64}) -> (dPhi,Phi,alpha,beta) -> begin

    #uses barycentric weights
    N = length(x)-1
    M = length(y)-1

    #set all to zeros
    @. dPhi = 0.0
    #println(length(dPhi))

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii=i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            #dPhi[nn] = 0.0 #dPhi[ii,jj] = 0.0
            for k in 1:(N-1) #sum over x
                kk = k+1
                nx = kk + (jj-1)*(N+1) #for Phi[k,j]
                dPhi[nn] = dPhi[nn] - wy[jj]*CNGLaplacianElement(x,wx,wbx, i,k)*Phi[nx]
            end
            for k in 1:(M-1) #sum over y
                kk = k+1
                ny = ii + (kk-1)*(N+1)
                dPhi[nn] = dPhi[nn] - wx[ii]*CNGLaplacianElement(y,wy,wby, j,k)*Phi[ny]
                #dPhi[nn] = dPhi[nn] + PolynomialDerivativeMatrixElement(y, j,k)*PolynomialDerivativeMatrixElement(y, j,k)*Phi[ny]
            end

        end
    end


end


#alg 78
#ch 5.2.2.4 pg 183
function Psi_Xi(k::Int64, l::Int64, eta::Int64)

    dPsi = zero(Float64)
    dPsi = -(-1.0)^k * (1.0-l-(-1.0)^l * eta)

    return dPsi

end

function Psi_Eta(k::Int64, l::Int64, xi::Int64)

    dPsi = zero(Float64)
    dPsi = -(-1.0)^l * (1.0-k-(-1.0)^k * xi)

    return dPsi

end

function gnmkl(dx::Float64, dy::Float64, n::Int64, m::Int64,k::Int64, l::Int64, r::Int64, s::Int64)

    g = Psi_Xi(k, l, s)*Psi_Xi(n, m, s)/dx^2 + Psi_Eta(k, l, r)*Psi_Eta(n, m, r)/dy^2
    return g

end

function LocalStiffnessMatrix(dx::Float64, dy::Float64)

    S = zeros(Float64, (4,4))
    t = 0.0

    for m in 0:1
        for n in 0:1
            p = n + 2*m + 1
            for l in 0:1
                for k in 0:1
                    q = k+2*l+1
                    # t = 0.0
                    # for s in 0:1
                    #     for r in 0:1

                    #         t = t + Psi_Xi(k, l, s)*Psi_Xi(n, m, s)/dx^2 + Psi_Eta(k, l, r)*Psi_Eta(n, m, r)/dy^2

                    #     end
                    # end
                    t = gnmkl(dx, dy,n,m, k, l, 0,0) + gnmkl(dx, dy,n,m, k, l, 1,0) + gnmkl(dx, dy, n,m,k, l, 0,1)+ gnmkl(dx, dy, n,m,k, l, 1,1)
                    S[p,q] = dx*dy*t/4.0
                end
            end
        end
    end

    return S

end

function StencilCoefficients(II::Int64, JJ::Int64, x::Array{Float64}, y::Array{Float64})


    Sp = Array{Any}(nothing, 4)
    C = zeros(Float64, (3,3))
    
    for m in 0:1
        for n in 0:1
            #p = n + 2*m + 1
            Sp[n+2*m+1] = LocalStiffnessMatrix(x[II-n+1+1] - x[II-n+1], y[JJ-m+1+1] - y[JJ-m+1])

        end
    end
    for m in 0:1
        for n in 0:1
            p = n + 2*m+1#p = n + 2*m+1
            # for l in -m:(-m+1)
            #     for k in -n:(-n+1)
            for k in -n:(-n+1)
                for l in -m:(-m+1)
                    q= (k+n) + 2*(l+m) + 1
                    #println(m,":::",n,"|||",k,":::",l,"|||",p,":::",q,"|||",k+2, ":::",l+2)
                    C[k+2,l+2] = C[k+2,l+2] + Sp[p][p,q]
                end
            end
        end
    end

    return C

end


#function FiniteElementPreconditioner!(X,B)
FEPreLinMaps = ( x::Array{Float64}, y::Array{Float64}) -> LinearMap((length(x))*(length(y)); issymmetric=true, ismutating=true) do U,B
    
    N = length(x) - 1
    M = length(y) - 1
    Cijkl = zeros(Float64,(3,3))

    @. B = 0.0
    
    for j in 1:M-1 #do not do y boundaries
        jj = j+1
        for i in 1:N-1 #do not do x boundaries
            ii = i+1
            n = ii + (jj-1)*(N+1)
            Cijkl .= StencilCoefficients(i, j, x, y)
            for l in -1:1
                ll=l+2
                for k in -1:1
                    #eqn 5.96
                    kk = k+2
                    nn = ii+k + (jj+l-1)*(N+1)
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        
                        #A[n,nnn] = A[n,nnn] +  Cijkl[kk,ll]
                        B[n] = B[n] + Cijkl[kk,ll] * U[nn]
                    # end

                end
            end
            # B[n] = B[n] + Cijkl[1,1] * U[ii-1 + (jj-1-1)*(N+1)] + Cijkl[2,1] + U[ii-0 + (jj-1-1)*(N+1)]+ Cijkl[3,1] * U[ii+1 + (jj-1-1) * (N+1)] + Cijkl[1,2] * U[ii-1 + (jj-0-1)*(N+1)] + Cijkl[2,2] * U[ii-0 + (jj-0-1)*(N+1)] + Cijkl[3,2] * U[ii+1 + (jj-0-1)*(N+1)] + Cijkl[1,3] * U[ii-1 + (jj+1-1)*(N+1)] + Cijkl[2,3] * U[ii-0 + (jj+1-1)*(N+1)] + Cijkl[3,3] * U[ii+1 + (jj+1-1)*(N+1)]

        end
    end
    return B

end


FEpreLinOp = (x::Array{Float64}, y::Array{Float64}) -> (U,B,alpha,beta) -> begin
    
    N = length(x) - 1
    M = length(y) - 1
    Cijkl = zeros(Float64,(3,3))

    @. B = 0.0

    for j in 1:M-1 #do not do y boundaries
        jj = j+1
        for i in 1:N-1 #do not do x boundaries
            ii = i+1
            n = ii + (jj-1)*(N+1)
            Cijkl .= StencilCoefficients(i, j, x, y)
            for l in -1:1
                for k in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        #eqn 5.96
                        nn = ii+k + (jj+l-1)*(N+1)
                        B[n] = B[n] + Cijkl[k+2,l+2] * U[nn]
                    # end

                end
            end

        end
    end
    #return B

end


#alg 79
#5.2.2.1 pg 186
#SSOR Sweep with Finite ELement Preconditioner using LinearMaps.jl
SSORSweepLinMaps = ( x::Array{Float64}, y::Array{Float64}, omega::Float64) -> LinearMap((length(x))*(length(y)); issymmetric=true,ismutating=true) do Z,R
    
    N = length(x)-1
    M = length(y)-1
    
    @. Z = 0.0
    s = 0.0
    Cijkl = zeros(Float64,(3,3))

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii = i+1
            n = ii + (jj-1)*(N+1)
            #Z[n] = 0.0
            s = 0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for k in -1:1
                for l in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        nn = (i+k)+1 + (j+l+1-1)*(N+1)#ii+k + (jj+l-1)*(N+1)
                        s = s + Cijkl[k+2,l+2]*Z[nn]
                    # end
                end
            end
            Z[n] = Z[n] + omega * (R[n] - s)/Cijkl[0+2,0+2]
        end
    end

    for j in M-1:-1:1
        jj=j+1
        for i in N-1:-1:1
            ii = i+1
            n = ii + (jj-1)*(N+1)
            s=0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for k in -1:1
                for l in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        nn = (i+k)+1 + (j+l+1-1)*(N+1) #ii+k + (jj+l-1)*(N+1)
                        s = s + Cijkl[k+2,l+2]*Z[nn]
                    # end
                end
            end
            Z[n] = Z[n] + omega * (R[n] - s)/Cijkl[0+2,0+2]
        end
    end

   
    return Z   


end


#alg 79
#5.2.2.1 pg 186
#SSOR Sweep with Finite ELement Preconditioner using LinearOperators.jl
SSORSweepLinOp = ( x::Array{Float64}, y::Array{Float64}, omega::Float64) ->(Z,R,alpha,beta) -> begin
    
    N = length(x)-1
    M = length(y)-1
    
    @. Z = 0.0
    s = 0.0
    Cijkl = zeros(Float64,(3,3))

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii = i+1
            n = ii + (jj-1)*(N+1)
            #Z[n] = 0.0
            s = 0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for k in -1:1
                for l in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        nn = (i+k)+1 + (j+l+1-1)*(N+1)#ii+k + (jj+l-1)*(N+1)
                        s = s + Cijkl[k+2,l+2]*Z[nn]
                    # end
                end
            end
            Z[n] = Z[n] + omega * (R[n] - s)/Cijkl[0+2,0+2]
        end
    end

    for j in M-1:-1:1
        jj=j+1
        for i in N-1:-1:1
            ii = i+1
            n = ii + (jj-1)*(N+1)
            s=0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for k in -1:1
                for l in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        nn = (i+k)+1 + (j+l+1-1)*(N+1) #ii+k + (jj+l-1)*(N+1)
                        s = s + Cijkl[k+2,l+2]*Z[nn]
                    # end
                end
            end
            Z[n] = Z[n] + omega * (R[n] - s)/Cijkl[0+2,0+2]
        end
    end


end


#Algorithm 35
function CoarseToFineInterpolation2Da(x::Array{Float64},  y::Array{Float64}, f::Matrix{Float64}, xx::Array{Float64},  yy::Array{Float64},xtol::Float64)

    wx = BarycentricWeights(x)
    wy = BarycentricWeights(y)
    #Tnewold = PolynomialInterpolationMatrix(x, xx, w,xtol)
    
    Nold = length(x)
    Mold = length(y)
    Nnew = length(xx)
    Mnew = length(yy)

    #for j in 1:Mold
    #end
    #just need matrix multiplication by vector
    #Fnewold = [Tnewold * f[:,j] for j in 1:Mold]
    # println("Phi : ",size(f))
    T = PolynomialInterpolationMatrix(x, xx, wx,xtol)
    # println("T : ",size(T))
    # println("phi[:,j] : ", size(f[:,1]))
    F1 = zeros(Float64, (Nnew,Mold))
    # println("F1 : ", size(F1))
    # println("F1[:,j] : ", size(F1[:,1]))
    for j in 1:Mold
        F1[:,j] .= T * f[:,j]
    end
    # println("F : ", size(F1))

    #w = BarycentricWeights(xx)
    T = PolynomialInterpolationMatrix(y, yy, wy,xtol)
    F2 = zeros(Float64, (Nnew,Mnew))
    for i in 1:Nnew
        F2[i,:] .= T * F1[i,:] 
    end

    return F2
    

end


#collocation Transport Operator
ColTranOp = (x::Array{Float64}, y::Array{Float64}, wbx::Array{Float64}, wby::Array{Float64}, u::Function, v::Function) -> LinearMap((length(x))*(length(y)); ismutating=true) do qdPhi,Phi
# function LaplaceOperator!(Phi::Array{Float64}, dPhi::Array{Float64}, x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64},  wbx::Array{Float64}, wby::Array{Float64}, DirBound1::Function, DirBound2::Function, DirBound3::Function,  DirBound4::Function)
    #uses barycentric weights
    N = length(x)-1
    M = length(y)-1

    #set all to zeros
    @. qdPhi = 0.0
    qdphi = 0.0 #dummy variable to store summation

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii=i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            qdphi = 0.0
            for k in 0:N #sum over x
                kk = k+1
                nx = kk + (jj-1)*(N+1)
                qdphi = qdphi + PolynomialDerivativeMatrixElement(x, wbx,i,k)*Phi[nx]
            end
            #add x summation to qdPhi and multiply by u velocity field
            qdPhi[nn] = qdPhi[nn] + u(x[ii])*qdphi
            #reset dummy variable
            qdphi = 0.0 
            for k in 0:M #sum over y
                kk = k+1
                ny = ii + (kk-1)*(N+1)
                qdphi = qdphi + PolynomialDerivativeMatrixElement(y, wby,j,k)*Phi[ny]
            end
            #add v velocity field times y summation to qdPhi
            qdPhi[nn]= qdPhi[nn] + v(y[jj])*qdphi
            #finally multiply by weights for galerkin only
            #qdPhi[nn] = qdPhi[nn]*wx[ii]*wy[jj]

        end
    end


    return qdPhi


end


GalTranOp = (x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64}, wbx::Array{Float64}, wby::Array{Float64}, u::Function, v::Function) -> LinearMap((length(x))*(length(y)); ismutating=true) do qdPhi,Phi

    #uses barycentric weights
    N = length(x)-1
    M = length(y)-1

    #set all to zeros
    @. qdPhi = 0.0
    #println(length(dPhi))
    qdphi = 0.0 #dummy variable to store summation

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii=i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            qdphi = 0.0
            for k in 0:N #sum over x
                kk = k+1
                nx = kk + (jj-1)*(N+1)
                qdphi = qdphi + PolynomialDerivativeMatrixElement(x, wbx,i,k)*Phi[nx]
            end
            #add x summation to qdPhi and multiply by u velocity field
            qdPhi[nn] = qdPhi[nn] + u(x[ii])*qdphi
            #reset dummy variable
            qdphi = 0.0 
            for k in 0:M #sum over y
                kk = k+1
                ny = ii + (kk-1)*(N+1)
                qdphi = qdphi + PolynomialDerivativeMatrixElement(y, wby,j,k)*Phi[ny]
            end
            #add v velocity field times y summation to qdPhi
            qdPhi[nn]= qdPhi[nn] + v(y[jj])*qdphi
            #finally multiply by weights
            qdPhi[nn] = qdPhi[nn]*wx[ii]*wy[jj]

        end
    end


    return qdPhi


end


#LHS operator for Collocation Method with LinearMaps
LHSCollocationLinMaps = (N::Int64, M::Int64, nu::Float64, dt::Float64, LapD::FunctionMap{Float64}) -> LinearMap((N+1)*(M+1); ismutating=true) do dPhi,Phi

    dtnu611 = -6.0*dt * nu/11.0
    
    LPhi = LapD*Phi 
    
    @. dPhi = 0.0
    @. dPhi = Phi + dtnu611*LPhi
    # for j in 0:M
    #     jj=j+1
    #     for i in 0:N
    #         ii=i+1
    #         nn = ii + (jj-1)*(N+1)
    #         dPhi[nn] = Phi[nn] - dtnu611 * LPhi[nn]
    #     end
    # end


    return dPhi


end


#with LinearOperators
LHSCollocationLinOp = (N::Int64, M::Int64, nu::Float64, dt::Float64, LapD::FunctionMap{Float64}) -> (dPhi,Phi,alpha,beta) -> begin

    dtnu611 = -6.0*dt * nu/11.0
    
    LPhi = LapD*Phi 
    
    @. dPhi = 0.0
    @. dPhi = Phi + dtnu611*LPhi
    # for j in 0:M
    #     jj=j+1
    #     for i in 0:N
    #         ii=i+1
    #         nn = ii + (jj-1)*(N+1)
    #         dPhi[nn] = Phi[nn] - dtnu611 * LPhi[nn]
    #     end
    # end


end


# #LHS operator Galerkin Method with LinearMaps
# LHSGalerkinCollMaps = (N::Float64, M::Float64, wx::Array{Float64}, wy::Array{Float64}, nu::Float64, dt::Float64, LapD::FunctionMap{Float64}) -> LinearMap((length(x))*(length(y)); ismutating=true) do dPhi,Phi

#     dtnu611 = -6.0*dt * nu/11.0
#     @. dPhi = 0.0
#     dPhi .= dtnu611 .* LapD*Phi 
    
#     for j in 0:M
#         jj=j+1
#         for i in 0:N
#             ii=i+1
#             nn = ii + (jj-1)*(N+1)
#             dPhi[nn] = x[ii]*wy[jj]*Phi[nn] + dPhi[nn] 
#         end
#     end


#     return dPhi


# end


# #LHS operator Galerkin Method with LinearOperators
# LHSGalerkinLinOp = (N::Float64, M::Float64, wx::Array{Float64}, wy::Array{Float64}, nu::Float64, dt::Float64, LapD::FunctionMap{Float64}) -> (dPhi,Phi,alpha,beta) -> begin

#     dtnu611 = -6.0*dt * nu/11.0
#     @. dPhi = 0.0
#     dPhi .= dtnu611 .* LapD*Phi 
    
#     for j in 0:M
#         jj=j+1
#         for i in 0:N
#             ii=i+1
#             nn = ii + (jj-1)*(N+1)
#             dPhi[nn] = x[ii]*wy[jj]*Phi[nn] + dPhi[nn] 
#         end
#     end


#     #return dPhi


# end


#eqn 5.37 & 5.38
#Eqns 5.12-5.122

# hatAij(dt::Float64, nu::Float64, dxi::Float64, dxip1::Float64,dyj::Float64, dyjp1::Float64) = 1.0 - 6.0*nu*dt/11.0 * Aij(dxi, dxip1,dyj, dyjp1)#-2*(1/(dxi*dxip1) + 1/(dyj*dyjp1))
# hatBij(dt::Float64, nu::Float64,dxi::Float64, dxip1::Float64) = - 6.0*nu*dt/11.0 *Bij(dxi, dxip1)#2*1/(dxi*(dxi+dxip1))
# hatCij(dt::Float64, nu::Float64,dyj::Float64, dyjp1::Float64) = - 6.0*nu*dt/11.0 *Cij(dyj, dyjp1)#2*1/(dyj*(dyj+dyjp1))
# hatEij(dt::Float64, nu::Float64,dxi::Float64, dxip1::Float64) = - 6.0*nu*dt/11.0 *Eij(dxi, dxip1)#2*1/(dxip1*(dxi+dxip1))
# hatFij(dt::Float64, nu::Float64,dyj::Float64, dyjp1::Float64) = - 6.0*nu*dt/11.0 *Fij(dyj, dyjp1)#2*1/(dyjp1*(dyj+dyjp1))


IMEX3FDPreCon = (x::Array{Float64}, y::Array{Float64}, dt::Float64, nu::Float64) -> LinearMap((length(x))*(length(y)); ismutating=true) do Z,U

    N = length(x) - 1
    M = length(y) - 1
    @. Z = 0.0
    dtnu611 = -6.0*dt*nu/11.0

    for j in 1:M-1
        jj = j+1
        for i in 1:N-1
            ii = i+1
            nn = ii + (jj-1)*(N+1)
            nnim1 = ii-1 + (jj-1)*(N+1)
            nnjm1 = ii + (jj-1-1)*(N+1)
            nnip1 = ii+1 + (jj-1)*(N+1)
            nnjp1 = ii + (jj-1+1)*(N+1)
            
            Z[nn] = (1.0 + dtnu611*Aij( (x[ii]-x[ii-1]) , (x[ii+1]-x[ii]), (y[jj]-y[jj-1]), (y[jj+1]-y[jj]) ) )* U[nn] + dtnu611*Bij( (x[ii]-x[ii-1]) , (x[ii+1]-x[ii]) ) * U[nnim1] +  dtnu611*Cij( (y[jj]-y[jj-1]), (y[jj+1]-y[jj]) ) * U[nnjm1] + dtnu611*Eij( (x[ii]-x[ii-1]) , (x[ii+1]-x[ii])  ) * U[nnip1] + dtnu611*Fij( (y[jj]-y[jj-1]), (y[jj+1]-y[jj])  ) * U[nnjp1]

        end
    end

    return Z
    

end


#Transport Operator with Linear Maps
GalTranOpLinMaps = (x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64}, wbx::Array{Float64}, wby::Array{Float64}, u::Function, v::Function) -> LinearMap((length(x))*(length(y)); ismutating=true) do qdPhi,Phi

    #uses barycentric weights
    N = length(x)-1
    M = length(y)-1

    #set all to zeros
    @. qdPhi = 0.0
    #println(length(dPhi))
    qdphi = 0.0 #dummy variable to store summation

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii=i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            qdphi = 0.0
            for k in 0:N #sum over x
                kk = k+1
                nx = kk + (jj-1)*(N+1)
                qdphi = qdphi + PolynomialDerivativeMatrixElement(x, wbx,i,k)*Phi[nx]
            end
            #add x summation to qdPhi and multiply by u velocity field
            qdPhi[nn] = qdPhi[nn] + u(x[ii])*qdphi
            #reset dummy variable
            qdphi = 0.0 
            for k in 0:M #sum over y
                kk = k+1
                ny = ii + (kk-1)*(N+1)
                qdphi = qdphi + PolynomialDerivativeMatrixElement(y, wby,j,k)*Phi[ny]
            end
            #add v velocity field times y summation to qdPhi
            qdPhi[nn]= qdPhi[nn] + v(y[jj])*qdphi
            #finally multiply by weights
            qdPhi[nn] = qdPhi[nn]*wx[ii]*wy[jj]

        end
    end


    return qdPhi


end


#Transport Operator with Linear Maps
GalTranOpLinOp = (x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64}, wbx::Array{Float64}, wby::Array{Float64}, u::Function, v::Function) -> (qdPhi,Phi, alpha, beta) -> begin

    #uses barycentric weights
    N = length(x)-1
    M = length(y)-1

    #set all to zeros
    @. qdPhi = 0.0
    #println(length(dPhi))
    qdphi = 0.0 #dummy variable to store summation

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii=i+1
            nn = ii + (jj-1)*(N+1) #columnwise/fortran
            # nn = j + (i-1)*(M-1) #rowwise/c
            qdphi = 0.0
            for k in 0:N #sum over x
                kk = k+1
                nx = kk + (jj-1)*(N+1)
                qdphi = qdphi + PolynomialDerivativeMatrixElement(x, wbx,i,k)*Phi[nx]
            end
            #add x summation to qdPhi and multiply by u velocity field
            qdPhi[nn] = qdPhi[nn] + u(x[ii])*qdphi
            #reset dummy variable
            qdphi = 0.0 
            for k in 0:M #sum over y
                kk = k+1
                ny = ii + (kk-1)*(N+1)
                qdphi = qdphi + PolynomialDerivativeMatrixElement(y, wby,j,k)*Phi[ny]
            end
            #add v velocity field times y summation to qdPhi
            qdPhi[nn]= qdPhi[nn] + v(y[jj])*qdphi
            #finally multiply by weights
            qdPhi[nn] = qdPhi[nn]*wx[ii]*wy[jj]

        end
    end


    #return qdPhi


end


#LHS operator Galerkin Method with LinearMaps
LHSGalerkinMaps = (N::Int64, M::Int64, wx::Array{Float64}, wy::Array{Float64}, nu::Float64, dt::Float64, LapD::FunctionMap{Float64}) -> LinearMap((length(x))*(length(y)); ismutating=true) do dPhi,Phi

    dtnu611 = -6.0*dt * nu/11.0
    @. dPhi = 0.0
    dPhi .= LapD*Phi 
    
    for j in 0:M
        jj=j+1
        for i in 0:N
            ii=i+1
            nn = ii + (jj-1)*(N+1)
            dPhi[nn] = wx[ii]*wy[jj]*Phi[nn] + dtnu611*dPhi[nn] 
        end
    end


    return dPhi


end

#LHS operator Galerkin Method with LinearOperators
#LHSGalerkinLinOp
LHSGalerkinLinOp = (N::Int64, M::Int64, wx::Array{Float64}, wy::Array{Float64}, nu::Float64, dt::Float64, LapD::LinearOperator{Float64}) -> (dPhi,Phi,alpha,beta) -> begin

    dtnu611 = -6.0*dt * nu/11.0
    @. dPhi = 0.0
    dPhi .=  LapD*Phi 
    
    for j in 0:M
        jj=j+1
        for i in 0:N
            ii=i+1
            nn = ii + (jj-1)*(N+1)
            dPhi[nn] = wx[ii]*wy[jj]*Phi[nn] + dtnu611*dPhi[nn] 
        end
    end


    #return dPhi


end


#function FiniteElementPreconditioner!(X,B)
FEPreTAdvDiffLinMaps = ( x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64}, dt::Float64, nu::Float64) -> LinearMap((length(x))*(length(y)); issymmetric=true, ismutating=true) do Z,R
    
    N = length(x) - 1
    M = length(y) - 1
    Cijkl = zeros(Float64,(3,3))
    z = 0.0
    dtnu611 = - 6.0*dt*nu/11.0
    #one11 = 1.0/11.0

    @. Z = 0.0
    
    for j in 1:M-1 #do not do y boundaries
        jj = j+1
        for i in 1:N-1 #do not do x boundaries
            ii = i+1
            n = ii + (jj-1)*(N+1)
            z = 0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for l in -1:1
                ll=l+2
                for k in -1:1
                    #eqn 5.96
                    kk = k+2
                    nn = ii+k + (jj+l-1)*(N+1)
                    z = z + Cijkl[kk,ll] * R[nn]
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        
                        #A[n,nnn] = A[n,nnn] +  Cijkl[kk,ll]
                        #B[n] = B[n] + Cijkl[kk,ll] * U[nn]
                    # end

                end
            end
            # B[n] = B[n] + Cijkl[1,1] * U[ii-1 + (jj-1-1)*(N+1)] + Cijkl[2,1] + U[ii-0 + (jj-1-1)*(N+1)]+ Cijkl[3,1] * U[ii+1 + (jj-1-1) * (N+1)] + Cijkl[1,2] * U[ii-1 + (jj-0-1)*(N+1)] + Cijkl[2,2] * U[ii-0 + (jj-0-1)*(N+1)] + Cijkl[3,2] * U[ii+1 + (jj-0-1)*(N+1)] + Cijkl[1,3] * U[ii-1 + (jj+1-1)*(N+1)] + Cijkl[2,3] * U[ii-0 + (jj+1-1)*(N+1)] + Cijkl[3,3] * U[ii+1 + (jj+1-1)*(N+1)]
            Z[n] = wx[ii]*wy[jj]*R[n] + dtnu611*z

        end
    end
    return Z

end


#FE Preconditioner with LinearOperators.jl
FEpreTAdvDiffLinOp = ( x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64}, dt::Float64, nu::Float64) -> (Z,R, alpha, beta) -> begin
    
    N = length(x) - 1
    M = length(y) - 1
    Cijkl = zeros(Float64,(3,3))
    z = 0.0
    dtnu611 = - 6.0*dt*nu/11.0
    #one11 = 1.0/11.0

    @. Z = 0.0
    
    for j in 1:M-1 #do not do y boundaries
        jj = j+1
        for i in 1:N-1 #do not do x boundaries
            ii = i+1
            n = ii + (jj-1)*(N+1)
            z = 0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for l in -1:1
                ll=l+2
                for k in -1:1
                    #eqn 5.96
                    kk = k+2
                    nn = ii+k + (jj+l-1)*(N+1)
                    z = z + Cijkl[kk,ll] * R[nn]
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        
                        #A[n,nnn] = A[n,nnn] +  Cijkl[kk,ll]
                        #B[n] = B[n] + Cijkl[kk,ll] * U[nn]
                    # end

                end
            end
            # B[n] = B[n] + Cijkl[1,1] * U[ii-1 + (jj-1-1)*(N+1)] + Cijkl[2,1] + U[ii-0 + (jj-1-1)*(N+1)]+ Cijkl[3,1] * U[ii+1 + (jj-1-1) * (N+1)] + Cijkl[1,2] * U[ii-1 + (jj-0-1)*(N+1)] + Cijkl[2,2] * U[ii-0 + (jj-0-1)*(N+1)] + Cijkl[3,2] * U[ii+1 + (jj-0-1)*(N+1)] + Cijkl[1,3] * U[ii-1 + (jj+1-1)*(N+1)] + Cijkl[2,3] * U[ii-0 + (jj+1-1)*(N+1)] + Cijkl[3,3] * U[ii+1 + (jj+1-1)*(N+1)]
            Z[n] = wx[ii]*wy[jj]*R[n] + dtnu611*z

        end
    end
    #return Z

end


#alg 79
#looks like we have to make this into a linear map
#function SSORSweep!(z::Array{Float64},r::Array{Float64}, omega::Float64, x::Array{Float64}, y::Array{Float64})
SSORSweepTAdvDiffLinMaps = ( x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64}, dt::Float64, nu::Float64, omega::Float64) -> LinearMap((length(x))*(length(y)); issymmetric=true,ismutating=true) do Z,R
    
    N = length(x)-1
    M = length(y)-1
    
    @. Z = 0.0
    s = 0.0
    Cijkl = zeros(Float64,(3,3))
    dtnu611 = -6.0*dt*nu/11.0

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii = i+1
            n = ii + (jj-1)*(N+1)
            #Z[n] = 0.0
            s = 0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for k in -1:1
                for l in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        nn = (i+k)+1 + (j+l+1-1)*(N+1)#ii+k + (jj+l-1)*(N+1)
                        s = s + Cijkl[k+2,l+2]*Z[nn]
                    # end
                end
            end
            s = wx[ii]*wy[jj]*Z[n] + dtnu611*s
            Z[n] = Z[n] + omega * (R[n] - s)/(wx[ii]*wy[jj] + dtnu611*Cijkl[0+2,0+2])
        end
    end

    for j in M-1:-1:1
        jj=j+1
        for i in N-1:-1:1
            ii = i+1
            n = ii + (jj-1)*(N+1)
            s=0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for k in -1:1
                for l in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        nn = (i+k)+1 + (j+l+1-1)*(N+1) #ii+k + (jj+l-1)*(N+1)
                        s = s + Cijkl[k+2,l+2]*Z[nn]
                    # end
                end
            end
            s = wx[ii]*wy[jj]*Z[n] + dtnu611*s
            Z[n] = Z[n] + omega * (R[n] - s)/(wx[ii]*wy[jj] + dtnu611*Cijkl[0+2,0+2])
        end
    end


    return Z   


end


#alg 79
#looks like we have to make this into a linear map
#function SSORSweep!(z::Array{Float64},r::Array{Float64}, omega::Float64, x::Array{Float64}, y::Array{Float64})
SSORSweepTAdvDiffLinOp = ( x::Array{Float64}, y::Array{Float64}, wx::Array{Float64}, wy::Array{Float64}, dt::Float64, nu::Float64, omega::Float64)->(Z,R,alpha,beta) -> begin
    
    N = length(x)-1
    M = length(y)-1
    
    @. Z = 0.0
    s = 0.0
    Cijkl = zeros(Float64,(3,3))
    dtnu611 = -6.0 * dt*nu/11.0

    for j in 1:M-1
        jj=j+1
        for i in 1:N-1
            ii = i+1
            n = ii + (jj-1)*(N+1)
            #Z[n] = 0.0
            s = 0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for k in -1:1
                for l in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        nn = (i+k)+1 + (j+l+1-1)*(N+1)#ii+k + (jj+l-1)*(N+1)
                        s = s + Cijkl[k+2,l+2]*Z[nn]
                    # end
                end
            end
            s = wx[ii]*wy[jj]*Z[n] + dtnu611*s
            Z[n] = Z[n] + omega * (R[n] - s)/(wx[ii]*wy[jj] + dtnu611*Cijkl[0+2,0+2])
        end
    end

    for j in M-1:-1:1
        jj=j+1
        for i in N-1:-1:1
            ii = i+1
            n = ii + (jj-1)*(N+1)
            s=0.0
            Cijkl .= StencilCoefficients(i, j, x, y)
            for k in -1:1
                for l in -1:1
                    # if i+k>=1 && j+l >=1 && i+k<=N-1 && j+l<=M-1
                        nn = (i+k)+1 + (j+l+1-1)*(N+1) #ii+k + (jj+l-1)*(N+1)
                        s = s + Cijkl[k+2,l+2]*Z[nn]
                    # end
                end
            end
            s = wx[ii]*wy[jj]*Z[n] + dtnu611*s
            Z[n] = Z[n] + omega * (R[n] - s)/(wx[ii]*wy[jj] + dtnu611*Cijkl[0+2,0+2])
        end
    end



end