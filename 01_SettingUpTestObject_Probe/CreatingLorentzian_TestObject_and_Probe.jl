using CSV
using LinearAlgebra
using OffsetArrays
using Plots
using Tables

include(".../00_BasicFunctions/BasicFunctionsAndPlots.jl")
cd(".../01_SettingUpTestObject_Probe")

ω_end = 0.5;
ω_res = 0.0005;
ω = -ω_end:ω_res:ω_end;
ω_l = (length(ω) - 1)/2;
Δω = ω_res;
Δt = (1/Δω)*(2*π/(ω_l));
t = Δt .* (-ω_l:ω_l);

## Object Amplitude
#= creating 4 Lorentzian peaks:
    .            .
   . .  .    .  . .
...   .. .... ..   ...
4 peaks, symmetric along the y-axis
=#

no_of_peaks = 4;
no_of_centres = no_of_peaks+2;
if mod(no_of_peaks,2)==0
    global no_of_centres = no_of_centres+1;
end

div_cen =  (2*ω_end)/no_of_centres;
centres = div_cen .* (-no_of_peaks/2:no_of_peaks/2);

# first peak
c₋₁ = centres[1]; # centre
Γ₋₁ = 0.02; # width

# second peak
c₋₂ = centres[2];
Γ₋₂ = 0.05;

#third peak
c₁ = centres[4];   # centre
Γ₁ = Γ₋₂;   # width

# fourth peak
c₂ = centres[5];   # centre
Γ₂ = Γ₋₁;   # width

# Lorentzian Function
"""
    Lorentzian(x, x₀, Γ)
    
    x ∈ ℝⁿ is the domain
    x₀ ∈ ℝ is the centre of the curve
    Γ ∈ ℝ⁺ is the width of the curve
"""
function Lorentzian(x, x₀, Γ)
    (1/π).*((0.5 .* Γ) ./ (((x .- x₀).^2) .+ (0.5 .* Γ).^2))
end
peak₋₂ = normalize(Lorentzian(ω, c₋₂, Γ₋₂));

# peak -1
peak₋₁ = normalize(Lorentzian(ω, c₋₁, Γ₋₁));

# peak 1
peak₁ = normalize(Lorentzian(ω, c₁, Γ₁));

peak₂ = normalize(Lorentzian(ω, c₂, Γ₂));

all_peaks = peak₋₂ .+ peak₋₁ .+ peak₁ .+ peak₂;


objectAmplitude = normalize(all_peaks);    # Object Amplitude

## Object Phase

"""
    step_func(x, centres)

creating a step function
                ....
            ....
        ....
    ....
"""
function step_func(x, centres)
    n = length(x);
    y = zeros(n);
    for i=1:n
        if x[i]<=centres[1]
            y[i] = -2.0;
        else
            y[i] = -1.0;
        end
        if x[i]>centres[2]
            y[i] = 0.0;
        end
        if x[i]>centres[4]
            y[i] = 1.0;
        end
        if x[i]>centres[5]
            y[i] = 2.0;
        end
    end

    return y;
end

a = 1.0;    # height of peak
b = 0.05;    # centre of peak
c² = 0.02;   # c = standard deviation : controls width

#phase = normalize(Gaussian(ω, a, b.+(Δω), c²));
#a = 10.0; phase = sinc_probe(ω, a);
phase = step_func(ω, centres);

## Complete Complex Object
object = objectAmplitude .* exp.(im .*phase);

## Probe
# Parameters
a = 1.0;    # height of peak
b = 0.0;    # centre of peak
c² = 0.02;   # c = standard deviation : controls width


"""
Function to define probe. 
    m is the index of offset.
"""
function probe(m, Δω)
    #normalize(Gaussian(ω, a, b.+(m*Δω), c²));  #old
    sinc_probe(ω.+(m*Δω), 10.0);
end

# 
"""
    Probe_set_function(M_max)
Function to create probe set of size M_max 
M_max = no. of probe readings. ideally odd.
"""
function Probe_set_function(M_max)
    m_end = Int64(floor(M_max/2));
    Δω = ω_end/((M_max-1));         # Old version: ω_end/(2*(M_max-1))

    Probe_set = zeros(length(ω), M_max) + im.*zeros(length(ω), M_max);

    Probe_set = OffsetArray(Probe_set, 1:length(ω), -m_end:m_end);
    for i=-m_end:m_end
        global Probe_set[:,i] = (probe(i, Δω)) .* exp.(im.*probe(i, Δω));
    end
    return Probe_set
end


M_max = 5;  #no. of probe readings. ideally an odd number
m_end = Int64(floor(M_max/2));
Probe_set = Probe_set_function(M_max);

##

obj_amp = objectAmplitude;
obj_phase = phase;

function plottingprobeset(probeset)
    M_max=size(probeset)[2];
    m_end = Int64(floor(M_max/2));
    colors = range(HSV(0,1,1), stop=HSV(-300,1,1), length=M_max);
    plot(size=(1800, 800),
        xlabel="ω", ylabel="|Pm̃(ω)|");
    for i=-m_end:m_end
    plot!(ω, abs.(probeset[:,i]), label=i, linewidth=2 , color=colors[i+m_mid+1]);
    end
    p_Probeset = plot!(legendtitle="m̃", legendfontsize=10, legendtitlefontsize=15);
    headerprobe = string.(-m_end:1:m_end);
    return  plot(p_Probeset, title="|Pm̃(ω)| set", titlefontsize=15, left_margin = 5mm, bottom_margin = 5mm)
end

## Writing and plotting

#=
using Tables, CSV

CSV.write("TestDomain.csv", Tables.table(ω), header=false)

p_objAmp = plot(ω, objectAmplitude, title="Object Amplitude"
        , xlabel ="ω", ylabel="|O(ω)|"
        , linewidth=2.5, label=false, color=:green, titlefontsize=15)
#CSV.write("TestObjectAmplitude.csv", Tables.table(objectAmplitude), header=false)

p_phase = plot(ω, phase, title="Object Phase"
        ,  xlabel ="ω", ylabel="θₒ(ω)"
        , linewidth=2.5, color=:darkgreen, titlefontsize=15);
#CSV.write("TestPhase.csv", Tables.table(phase), header=false)

plot(real(object), imag(object));
#CSV.write("TestObject.csv", Tables.table(object), header=false)

p_probe=plot(ω, abs.(Probe_set[:,0]), title="Probe", xlabel ="ω", ylabel="Probe(ω)")
#CSV.write("TestProbe.csv", Tables.table(vec(Probe_set[:,0])), header=false)

plot(p_objAmp, p_phase
        #, p_probe
    , layout=(1,2), legend=false, size = (1800,800)
    , left_margin = 5mm
    , bottom_margin = 5mm
    )

savefig("ObjectAmplitude_Phase_Probe.png")

plottingprobeset(Probe_set)
savefig("Probe_set.png")
=#
