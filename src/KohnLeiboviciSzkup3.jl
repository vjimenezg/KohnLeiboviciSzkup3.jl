module KohnLeiboviciSzkup3

## ##################################################
# This code solves the model in
# Financial Development and Trade Liberalization
# Kohn, Leibovici, Szkup 2020
####################################################

clearconsole()

# Dependencies

using LinearAlgebra
using DifferentialEquations, Sundials, SimpleDifferentialOperators, DiffEqCallbacks
using DataFrames, DataFramesMeta, CSV, CSVFiles, JSON # results caching
using Interpolations, QuadGK # integration
using NLsolve # root-finding
using NamedTupleTools, Parameters # named tuples
using Roots

# Model files

include("parameters.jl")


end
