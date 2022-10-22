using LinearAlgebra, OffsetArrays
using Plots
using CSV, Tables

cd("C://Users/avina/Documents/2020-2022_MSc_MathMods/401_Thesis/JuliaCode/01_SettingUpTestObject_Probe")

x_end = 0.5
x_res = 0.0005
x = -x_end:x_res:x_end;

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

div_cen =  (2*x_end)/no_of_centres;
centres = div_cen .* (-no_of_peaks/2:no_of_peaks/2)

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
Lorentzian(x, x₀, Γ) = (1/π).*((0.5 .* Γ) ./ (((x .- x₀).^2) .+ (0.5 .* Γ).^2));

# peak -2
peak₋₂ = normalize(Lorentzian(x, c₋₂, Γ₋₂));

# peak -1
peak₋₁ = normalize(Lorentzian(x, c₋₁, Γ₋₁));

# peak 1
peak₁ = normalize(Lorentzian(x, c₁, Γ₁));

peak₂ = normalize(Lorentzian(x, c₂, Γ₂));

all_peaks = peak₋₂ .+ peak₋₁ .+ peak₁ .+ peak₂;


objectAmplitude = normalize(all_peaks);    # Object Amplitude

## Object Phase
#=
creating a step function
                ....
            ....
        ....
    ....
=#

function step_func(x, centres)
    n = length(x);
    y = zeros(n);
    for i=1:n
        if x[i]<=centres[1]
            y[i] = -2.0
        else
            y[i] = -1.0
        end
        if x[i]>centres[2]
            y[i] = 0.0
        end
        if x[i]>centres[4]
            y[i] = 1.0
        end
        if x[i]>centres[5]
            y[i] = 2.0
        end
    end
    return y
end

phase = normalize(step_func(x, centres));

## Complete Complex Object

object = objectAmplitude .* exp.(im .*phase);

## Probe
# Parameters
a = 1.0;    # height of peak
b = 0.0;    # centre of peak
c² = 0.02;   # c = standard deviation : controls width

Gaussian(x, a, b, c²) = a .* exp.(-((x .- b).^2)/(2.0 .* c²))

M_max = 5;  #odd number
m_end = Int64(floor(M_max/2));
Δω = x_end/(2*(M_max-1));
probe(m, Δω) = normalize(Gaussian(x, a, b.+(m*Δω), c²));

Probe_set = zeros(length(x), M_max)

Probe_set = OffsetArray(Probe_set, 1:length(x), -m_end:m_end)
for i=-m_end:m_end
    global Probe_set[:,i] = probe(i, Δω)
end

## Writing and plotting

#=

CSV.write("TestDomain.csv", Tables.table(x), header=false)

p_objAmp = plot(x, objectAmplitude, title="Object Amplitude", xlabel ="x", ylabel="Obj(Amp)(ω)")
CSV.write("TestObjectAmplitude.csv", Tables.table(objectAmplitude), header=false)

p_phase = plot(x, phase, title="Object Phase",  xlabel ="x", ylabel="Obj(Phase)(ω)");
CSV.write("TestPhase.csv", Tables.table(phase), header=false)

plot(real(object), imag(object));
CSV.write("TestObject.csv", Tables.table(object), header=false)

p_probe=plot(x, probe(0.0,0.0), title="Probe", xlabel ="x", ylabel="Probe(ω)")
CSV.write("TestProbe.csv", Tables.table(probe), header=false)

plot(p_objAmp, p_phase, p_probe
        ,layout=(1,3), legend=false, size = (1800, 600)
        ,xlabel=" "
    )

savefig("ObjectAmplitude_Phase_Probe.png")

plot(size=(1400, 1000));
for i=-m_end:m_end
    plot!(Probe_set[:,i], label=i, linewidth=2);
end
p_Probeset = plot!(legendtitle="m");
headerprobe = string.(-m_end:1:m_end)
CSV.write("TestProbeSet.csv", Tables.table(parent(Probe_set)), header=headerprobe)
plot(p_Probeset)
savefig("Probe_set.png")
=#
