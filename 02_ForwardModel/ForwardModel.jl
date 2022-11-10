using FFTW, LinearAlgebra, Plots, ColorSchemes, CSV, Tables, DataFrames, OffsetArrays
# using ToeplitzMatrices;

cd(".../02_ForwardModel")

include(".../01_SettingUpTestObject_Probe/CreatingLorentzian_TestObject_and_Probe.jl")

obj_amp = objectAmplitude;
obj_phase = phase;
object = obj_amp .* exp.(im .* obj_phase);
probe0 = Probe_set[:,0]

## Forward model
#= |F[Pₘ(ω)∘O(ω)]|² = Iₘ(t)
 where  ∘ is the Hadamard (or Schur) Product or element wise multiplication
        Pₘ(ω) is the probe
        O(ω) is the object
=#
function Intensity_measurements(object, probeset)
    (n,m)=size(probeset);
    IntMeas = zeros(n,m);
    probeset = parent(probeset);
    for i=1:m
        prod_P_O = probeset[:,i] .* object;
        FFT_prod = fftshift(fft(prod_P_O));
        intensity = (abs.(FFT_prod)).^2;
        IntMeas[:,i] = intensity
    end
    return IntMeas
end # function
Int_Measured = Intensity_measurements(object, Probe_set)
mid = Int64(floor(length(x)/2));
int_plot_half = 35;
int_plot = mid-int_plot_half:mid+int_plot_half

(n,m)=size(Int_Measured);
m_mid = Int64(floor(m/2));
m_range =-m_mid:m_mid;

colors =range(HSV(0,1,1), stop=HSV(-360,1,1), length=m)
plot()
for i=1:m_mid+1
    plot!(int_plot, Int_Measured[int_plot,i], linestyle=:solid, linewidth=2, color=colors[i], label=m_range[i])
end
plot!()
for i=m_mid+2:m
    plot!(int_plot, Int_Measured[int_plot,i], marker=true, linealpha=0, markersize=4, markerstrokewidth=0, markercolor=colors[i], label=m_range[i])
end
plot!(size=(1090,727), legendtitle="m", title="Iₘ(t)")
