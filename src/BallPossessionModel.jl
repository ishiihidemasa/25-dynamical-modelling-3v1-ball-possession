module BallPossessionModel
using LinearAlgebra
import Statistics: mean
using Random
using StaticArrays
using DifferentialEquations
using LaTeXStrings

include("BPMSimulation.jl")
include("BPMAnalysis.jl")
end