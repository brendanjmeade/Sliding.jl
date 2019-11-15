using DifferentialEquations
using PyCall
using PyPlot
using Infiltrator

function aginglaw(v, θ, dc)
    return 1 - θ * v / dc
end

function sliplaw(v, θ, dc)
    return -v * θ / dc * log(v * θ / dc)
end

function calcdadv()
    dadv = 0
    return dadv
end

function calcdvθclassic!(dvθ, vθ, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ, statelaw = p
    θ = vθ[1]
    v = vθ[2]
    dvθ[1] = statelaw(v, θ, dc)
    dvθ[2] = 1 / (η / σn + a / v) * (μ * (Vp - v) / (L * σn) - b * dvθ[1] / θ)
    return nothing
end


function calcdvθ!(dvθ, vθ, p, t)
    dc, η, σn, a, b, μ, vp, L, ρ, statelaw, fstar, vstar, θstar = p
    θ = vθ[1]
    v = vθ[2]
    dθprev = dvθ[1]
    dvprev = dvθ[2]
    dvθ[1] = statelaw(v, θ, dc)
    dadt = 0 # calcdadv(a, v) * dvprev
    dbdt = 0
    numclassic = (μ * (vp - v) / (L * σn) - b * dvθ[1] / θ)
    numab = -dadt * log(v / vstar) - dbdt * log(θ / θstar)
    numσn = 0
    denom = (η / σn + a / v)
    dvθ[2] = (numclassic + numab + numσn) / denom
    return nothing
end

function plottimeseries(sol, siay, titlelabel)
    t1 = [x / siay for x in sol.t]
    θ1 = [x[1] for x in sol.u]
    v1 = [x[2] for x in sol.u]

    figure(figsize = (12, 6))
    subplot(2, 2, 1)
    plot(t1, θ1, "-b", linewidth = 0.5)
    yscale("log")
    ylabel(L"\theta")

    subplot(2, 2, 2)
    plot(1:1:length(t1), θ1, "-b", linewidth = 0.5)
    yscale("log")
    ylabel(L"\theta")

    subplot(2, 2, 3)
    plot(t1, v1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (years)")
    ylabel("v (m/s)")

    subplot(2, 2, 4)
    plot(1:1:length(t1), v1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (step)")
    ylabel("v (m/s)")

    suptitle(titlelabel)

    return nothing
end

function sliding()
    # Model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 2000.0)
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
    L = 60 * 1e3
    a = 0.015
    b = 0.02
    vp = 1e-9
    σn = 30e6
    dc = 0.2
    abstol = 1e-4
    reltol = 1e-4

    # Parameters need for d(a, b σn)/dt
    fstar = 0.6 # ???
    vstar = vp # ???
    θstar = 1e9 # ???
    
    # Time integrate - classic
    icsclassic = [1e8; vp / 1000]
    pclassic = (dc, η, σn, a, b, μ, vp, L, ρ, aginglaw, fstar, vstar, θstar)
    probclassic = ODEProblem(calcdvθclassic!, icsclassic, tspan, pclassic)
    solclassic = solve(probclassic, RK4(), abstol = abstol, reltol = reltol)
    plottimeseries(solclassic, siay, "RSF classic")

    # Time integrate - a, b can evolve
    icsab = [1e8 ; vp / 1000; a ; b]
    # p = (dc, η, σn, a, b, μ, vp, L, ρ, aginglaw, fstar, vstar, θstar)
    # prob = ODEProblem(calcdvθclassic!, ics, tspan, p)
    # sol = solve(prob, RK4(), abstol = abstol, reltol = reltol)
    # plottimeseries(sol, siay, "RSF classic")

    
    return nothing
end
sliding()
