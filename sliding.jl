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

function calcdadv(v)
    amin = 0.045
    ahighdiff =  0.0319 - amin
    δ = 1e-2 # m/s
    # a = amin + ahighdiff * 0.5 * (1 + tanh((v - 5e-3) / δ))
    dadv = ahighdiff * 0.5 * (sech((v - 5e-3) / δ))^2 / δ
    return dadv
end

function calcdbdv(v)
    bmin = 0.085
    bhighdiff =  0.2811 - bmin
    δ = 1e-2 # m/s
    # b = bmin + bhighdiff * 0.5 * (1 + tanh((v - 5e-3) / δ))
    dbdv = bhighdiff * 0.5 * (sech((v - 5e-3) / δ))^2 / δ
    return dbdv
end

function calcdvθclassic!(du, u, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ, statelaw = p
    θ = u[1]
    v = u[2]
    du[1] = statelaw(v, θ, dc)
    du[2] = 1 / (η / σn + a / v) * (μ * (Vp - v) / (L * σn) - b * du[1] / θ)
    return nothing
end

function calcdvθab!(du, u, p, t)
    dc, η, σn, μ, vp, L, ρ, statelaw, fstar, vstar, θstar = p
    θ = u[1]
    v = u[2]
    a = u[3]
    b = u[4]
    dθ = du[1]
    dv = du[2]
    du[1] = statelaw(v, θ, dc)
    da = calcdadv(v) * dv^2 # Pretty sure this square is wrong
    db = calcdbdv(v) * dv^2 # Pretty sure this square is wrong
    numclassic = (μ * (vp - v) / (L * σn) - b * du[1] / θ)
    numab = -da * log(abs(v) / vstar) - db * log(abs(θ) / θstar)
    numσn = 0
    denom = (η / σn + a / v)
    du[2] = (numclassic + numab + numσn) / denom
    du[3] = da
    du[4] = db
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
    tspan = (0.0, siay * 10000.0)
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
    L = 60 * 1e3
    a = 0.015
    b = 0.02
    vp = 1e-9
    σn = 30e6
    dc = 0.1
    abstol = 1e-4
    reltol = 1e-4

    # Parameters need for d(a, b σn)/dt
    fstar = 0.6 # ???
    vstar = vp # ???
    θstar = 1e9 # ???
    
    # Time integrate - classic
    icsclassic = [1e8; vp / 1000]
    pclassic = (dc, η, σn, a, b, μ, vp, L, ρ, aginglaw)
    probclassic = ODEProblem(calcdvθclassic!, icsclassic, tspan, pclassic)
    solclassic = solve(probclassic, RK4(), abstol = abstol, reltol = reltol)
    plottimeseries(solclassic, siay, "RSF (classic)")

    # Time integrate - a, b can evolve
    icsab = [1e8 ; vp / 1000 ; a ; b]
    pab = (dc, η, σn, μ, vp, L, ρ, aginglaw, fstar, vstar, θstar)
    probab = ODEProblem(calcdvθab!, icsab, tspan, pab)
    solab = solve(probab, RK4(), abstol = abstol, reltol = reltol)
    plottimeseries(solab, siay, "RSF (ab)")
    
    return nothing
end
sliding()
