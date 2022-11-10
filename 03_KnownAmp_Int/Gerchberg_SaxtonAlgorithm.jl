using FFTW, LinearAlgebra, Plots, ColorSchemes, StatsBase

## Reference:
#=
    Phase retrival algorithms: a comparison
        J.R. Fienup
        1 Aug 1982
        Applied Optics
        Vol. 21
        No. 15
=#

## Setting up system and directory

cd(".../03_KnownAmp-Int")

include(".../01_SettingUpTestObject_Probe/CreatingLorentzian_TestObject_and_Probe.jl")

#=
imports:    ω: freq domain
            Δω: resolution of ω
            t: time domain
            Δt: resolution of t

            object: object ∈ ℂ
            obj_amp: object Amplitude
            obj_phase: object phase
            Probe_set : set of dimension 1:2001 × -2:2
=#

## Setting up readings

# Measurements: abs_f = obj_amp and abs_F = |Fourier[f(x)]|

abs_f = obj_amp;
abs_F = abs.(fftshift(fft(object)));
plot((abs_F)./maximum(abs_F))
plot!((obj_amp)./maximum(obj_amp))
N = length(ω);
## Gerchberg - Saxton algorithm
#=
N : grid length
Steps:
    1. Fourier transform an estimate of the object
    2. Replace the modulus of the resulting computed Fourier transform with the measured Fourier modulus to form an estimate of the Fourier transform
    3. Inverse Fourier transform the estimate of the Fourier transform
    4. Replace the modulus of the resulting computed image with the measured object modulus to form a new estimate of the object.

Quantity    Estimate of
    gₖ          f   : Object
    θₖ          η   : Phase of object
    G'ₖ         F   : FFT(Object)
    ϕₖ          ψ   : Phase of FFT(Object)

Constraints:
Fourier domain constraint: Bₖ = (1/√N) Σ|Gₖ - F|²
Function domain contraint: Eₖ = Σ|f - g'ₖ|²
=#
