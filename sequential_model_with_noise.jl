using PointPatternStatistics
using SequentialPointProcesses
using Distributions
using Optim
using ComponentArrays
using UnPack
include("readpattern.jl")

# Read point pattern and observation window for one subject
pp = readpattern(20)

# Define model and loglikelihood
M1(R, kappa) = Softcore(d -> (R/d)^(2/kappa))
M2(R, kappa, theta) = Mixture(M1(R, kappa), SequentialPointProcesses.Uniform(), theta)
function log_likelihood(x)
    @unpack R, kappa, theta = x
    R <= 0.0 && return -Inf
    0.0 <= kappa <= 1.0 || return -Inf
    logpdf(M2(R, kappa, theta), pp, (nx=120,))
end

# Maximize loglikelihood
@time o2 = maximize(log_likelihood, ComponentVector(R=70.0, kappa=0.4, theta=0.1))
p2 = Optim.maximizer(o2)

# Generate samples from the Softcore model using the estimated parameters
pps = [rand(M1(p2.R, p2.kappa), pp.window, length(pp.data)) for _ in 1:2500]

plot(plot(pp, labels="", title="Data"), plot(pps[1], labels="", title="Simulation"))

# Compute pair correlation functions for all samples and observed pattern
r = 0:500
gsim = [pcf(pp, r)[2:end] for pp in pps]
gobs = pcf(pp, r)[2:end]
r = 1:500
plot(r, gobs, label="obs")
plot!(r, mean(gsim), label="mean")

# Compute global envelope from simulated pair correlation functions
(hi, lo, alpha_interval) = globalenvelope(gsim)
plot!(r, lo, label="lo")
plot!(r, hi, label="hi")

# Use ggplot2 to plot the envelope
using RCall
using DataFrames
@rlibrary(ggplot2)
ggplot(DataFrame(r=1:500, hi=hi, lo=lo, obs=gobs), aes(x=:r, ymax=:hi, y=:obs, ymin=:lo)) +
    geom_ribbon(fill="grey70") + geom_line() + geom_hline(yintercept=1) +
    ggtitle("Global 95% envelope for pair-correlation function")

# Plot simulated pattern using base R (requires spatstat to be installed)
R"plot"(pps[1], main="Simulated pattern")
