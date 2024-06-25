module Part2a

#this module imports the functions from Part 1 Chapters 1-3 and Part 2 Chapter 4

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


include("Part1.jl")
include("Part2Ch4.jl")

export FourierCollocationTimeDerivative, AdvectionDiffusionTimeDerivative!, EvaluateFourierGalerkinSolution!, FastConvolutionSum!, mthFastChebyshevDerivative, mthFastChebyshevDerivativeC, LegendreTransformNaive, InverseLegendreTransformNaive, FastLegendreTransform, FastInverseLegendreTransform, FastTransformsLegendre, FastInverseTransformsLegendre, ModifiedLegendreBasis, EvaluateLegendreGalerkinSolution, alpha, gamma, mu, beta, initTMatrix, initTMatrixUpperBidiagonal, ModifiedCoefsFromLegendreCoefs, ModifiedCoefsFromLegendreCoefsMatrix, CGDerivativeMatrix, NodalDiscontinuousGalerkinConstruct, ComputeDGDerivative, InterpolateToBoundar, DGTimeDerivative

export FourierInterpolantFromModes, FourierInterpolantFromNodes!, LegendreDerivativeCoefficients!, ChebyshevDerivativeCoefficients!, FourierDerivativeByFFT, FourierDerivativeMatrix, LegendrePolynomial, LegendreSum, ChebyshevPolynomialTrig, ChebyshevPolynomialIter, SPLegendre, LegendrePolynomialandDerivative, LegendreGaussNodesAndWeights, qAndLEvaluation, LegendreGaussLobattoNodesAndWeights, ChebyshevGaussNodesAndWeights, ChebyshevGaussLobattoNodesAndWeights, FastChebyshevTransform1, InvFastChebyshevTransform1, FastChebyshevTransform2, InvFastChebyshevTransform2, mthChebyshevDerivativeCoefficients, BarycentricWeights, LagrangeInterpolation, PolynomialInterpolationMatrix, LagrangeInterpolatingPolynomials, CoarseToFineInterpolation2D, LagrangeInterpolantDerivative, PolynomialDerivativeMatrix, mthOrderPolynomialDerivativeMatrix, EOMatrixDerivative, FastChebyshevDerivative









end # module Part2a
