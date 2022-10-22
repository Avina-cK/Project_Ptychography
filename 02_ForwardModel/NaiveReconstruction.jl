## Naive reconstruction
using FFTW, LinearAlgebra, Plots, ColorSchemes, CSV, Tables, DataFrames
# using ToeplitzMatrices;

cd("C://Users/avina/Documents/2020-2022_MSc_MathMods/401_Thesis/JuliaCode/02_ForwardModel")
#=
imports:    x: domain
            object: object ∈ ℂ
            Probe_set : set of dimension 1:2001 × -2:2
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
    plot!(x, normalize(angle.(Results[:,i])), alpha=0.2
    , marker=true, markeralpha=0.3, markersize=2
    , markercolor=colors[i], markerstrokewidth=0
    , markerstrokecolor=colors[i], linecolor=colors[i]
    , label="NR, probe shift "*string(n_range[i]))
end
plot!(x, normalize(angle.(object)), marker=true
    , markeralpha=0.6, markersize=2, markerstrokewidth=0
    , markerstrokecolor=:orange, markercolor=:orange
    ,linecolor=:orange, label="Actual Object")
plot!(legend=:topleft, legendtitle = "Naïve reconstruction (NR)"
    , size=(1400,800)  , title="Comparing Phase")

#savefig("ComparingPhase.png")

plot(abs.(result./norm(result)), linealpha=0.8, marker=true, markeralpha=0.2, markersize=2, label="Absolute Value (Naïve reconstruction)",markerstrokecolor=:blue, markercolor=:blue, linewidth=2)
plot!(abs.(object./norm(object)), linealpha=0.4, marker=true, markeralpha=0.8, markersize=2, label="Absolute Value (Object)",markerstrokecolor=:orange, markercolor=:orange, linewidth=2)
plot!(legend=:top, ylabel="Absolute Value", size=(1400,800))

#savefig("ComparingAbsValue.png")
