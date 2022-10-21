using LinearAlgebra, Plots, FFTW, CSV, Tables

cd("...")

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
A₋₁ = 0.4;  # Amplitude

# second peak
c₋₂ = centres[2];
Γ₋₂ = 0.05;
A₋₂ = 1.0;

#third peak
c₁ = centres[4];   # centre
Γ₁ = Γ₋₂;   # width
A₁ = A₋₂;  # Amplitude

# fourth peak
c₂ = centres[5];   # centre
Γ₂ = Γ₋₁;   # width
A₂ = A₋₁;  # Amplitude

# Lorentzian Function
Lorentzian(x, x₀, Γ) = (1/π).*((0.5 .* Γ) ./ (((x .- x₀).^2) .+ (0.5 .* Γ).^2));

# peak -2
peak₋₂ = Lorentzian(x, c₋₂, Γ₋₂);
plot(x, peak₋₂);

# peak -1
peak₋₁ = Lorentzian(x, c₋₁, Γ₋₁);

# peak 1
peak₁ = Lorentzian(x, c₁, Γ₁);

peak₂ = Lorentzian(x, c₂, Γ₂);

all_peaks = peak₋₂ .+ peak₋₁ .+ peak₁ .+ peak₂;


objectAmplitude = (all_peaks);    # Object Amplitude
p_objAmp = plot(x, objectAmplitude, title="Object Amplitude", xlabel ="x", ylabel="Obj(Amp)(ω)");
CSV.write("TestObjectAmplitude.csv", Tables.table(objectAmplitude), header=false)

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

phase = (step_func(x, centres));
p_phase = plot(x, phase, title="Object Phase",  xlabel ="x", ylabel="Obj(Phase)(ω)");

CSV.write("TestPhase.csv", Tables.table(phase), header=false)

## Complete Complex Object

object = objectAmplitude .* exp.(im .*phase);
plot(real(object), imag(object));
CSV.write("TestObject.csv", Tables.table(object), header=false)

## Probe
# Parameters
a = 1.0;    # height of peak
b = 0.0;    # centre of peak
c² = 0.02;   # c = standard deviation : controls width

Gaussian(x, a, b, c²) = a .* exp.(-((x .- b).^2)/(2.0 .* c²))

probe = normalize(Gaussian(x, a, b, c²));
p_probe=plot(x, probe, title="Probe", xlabel ="x", ylabel="Probe(ω)")

CSV.write("TestProbe.csv", Tables.table(probe), header=false)

plot(p_objAmp, p_phase, p_probe
        ,layout=(1,3), legend=false, size = (1200, 400)
        ,xlabel=" "
    )

savefig("ObjectAmplitude_Phase_Probe.png")
