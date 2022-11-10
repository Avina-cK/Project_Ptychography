## Naive reconstruction
using FFTW, LinearAlgebra, Plots, ColorSchemes, CSV, Tables, DataFrames
# using ToeplitzMatrices;

cd(".../02_ForwardModel")
include("ForwardModel.jl")

#=
imports:    ω: freq domain
            Δω: resolution of ω
            t: time domain
            Δt: resolution of t

            object: object ∈ ℂ
            obj_amp: object Amplitude
            obj_phase: object phase
            Probe_set : set of dimension 1:2001 × -2:2
            Int_Measured : set of measured intensities
=#

## Forward model
#= |F[Pₘ(ω)∘O(ω)]|² = Iₘ(t)
 where  ∘ is the Hadamard (or Schur) Product or element wise multiplication
        Pₘ(ω) is the probe
        O(ω) is the object
=#
# Int_Measured was imported
#take one intensity
(m,n) = size(Int_Measured);
n_mid = Int64(floor(n/2));
n_range =-n_mid:n_mid;

Results = zeros(m,n) .+ (im.*zeros(m,n));
for i=1:n
    intensity = Int_Measured[:,i];
    probe1 = parent(Probe_set)[:,i];
    sqrt_int = sqrt.(intensity);
    ifft_sqrt_int = (ifftshift(ifft(sqrt_int)));
    result = ifft_sqrt_int./probe1;
    global Results[:,i]=result;
end

plot()
colors =range(HSV(0,1,1), stop=HSV(-360,1,1), length=n)
for i=1:n
    plot!(ω, normalize(angle.(Results[:,i])), alpha=0.2
    , marker=true, markeralpha=0.3, markersize=2
    , markercolor=colors[i], markerstrokewidth=0
    , markerstrokecolor=colors[i], linecolor=colors[i]
    , label="NR, probe shift "*string(n_range[i]))
end
plot!(ω, normalize(angle.(object)), marker=true
    , markeralpha=0.6, markersize=2, markerstrokewidth=0
    , markerstrokecolor=:orange, markercolor=:orange
    ,linecolor=:orange, label="Actual Object")
plot!(legend=:topleft, legendtitle = "Naïve reconstruction (NR)"
    , size=(1400,800)  , title="Comparing Phase", xlabel="ω")
#savefig("ComparingPhase.png")

plot()
for i=1:n
    plot!(ω, abs.(normalize(Results[:,i]))
    , alpha=0.3, linewidth=2
    , marker=true, markeralpha=0.3, markersize=2
    , markercolor=colors[i], markerstrokewidth=0
    , markerstrokecolor=colors[i], linecolor=colors[i]
    , label="NR, probe shift "*string(n_range[i]))
end
plot!(ω, abs.(normalize(object)), linealpha=0.4, marker=true, markeralpha=0.8, markersize=2, label="Object",markerstrokecolor=:orange, markercolor=:orange, color=:orange, linewidth=2)
plot!(legend=:top, title="Comparing Absolute Value", size=(1400,800), legendtitle="Naïve reconstruction(NR)", xlabel="ω")

savefig("ComparingAbsValue.png")


# Next naïve step would be to take the average of all the results
Avg_Result = vec(mean(Results, dims=2));
plot(ω, abs.(normalize(Avg_Result)), linewidth=2
    , label="NR, Avg(Result)")
plot!(ω, abs.(normalize(object)), linewidth=2, label="Object")
plot!(legend=:top, title="Comparing Absolute Value"
    ,size=(1400,800), legendtitle="Naïve reconstruction(NR)"
    , xlabel="ω", ylabel="|O(ω)|")
savefig("ComparingAbsValue_AveragevsTrue.png")

plot(ω, (normalize(angle.(Avg_Result)))
    , label="NR, Avg(Result)", linealpha=0.3
    , marker=true, markeralpha=0.4
    , markersize=2, markerstrokewidth=0)
plot!(ω, (normalize(angle.(object))), label="Object"
    , linewidth=2)
plot!(legend=:topleft
    , legendtitle = "Naïve reconstruction (NR)"
    , size=(1400,800)  , title="Comparing Phase"
    , xlabel="ω", ylabel ="Phase(Object)")
savefig("ComparingPhase_AveragevsTrue.png")
