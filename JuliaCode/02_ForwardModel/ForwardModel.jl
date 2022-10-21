using FFTW, LinearAlgebra, Plots, ColorSchemes, CSV, Tables, DataFrames
# using ToeplitzMatrices;

cd("...")

test_object_amp = CSV.read("TestObjectAmplitude.csv", DataFrame, header=false);
test_object_phase = CSV.read("TestPhase.csv", DataFrame, header=false);
test_probe = CSV.read("TestProbe.csv", DataFrame, header=false);

obj_amp = Matrix(test_object_amp);
obj_phase = Matrix(test_object_phase);
object = obj_amp .* exp.(im .* obj_phase);
probe = Matrix(test_probe);

## Forward model
#= |F[Pₘ(ω)∘O(ω)]|² = Iₘ(t)
 where  ∘ is the Hadamard (or Schur) Product or element wise multiplication
        Pₘ(ω) is the probe
        O(ω) is the object
=#
prod_P_O = probe .* object;
intensity = (abs.(fftshift(fft(prod_P_O)))).^2
