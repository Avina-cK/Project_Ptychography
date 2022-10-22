using LinearAlgebra, Plots, FFTW

cd("C://Users/avina/Documents/2020-2022_MSc_MathMods/401_Thesis/JuliaCode/00_BasicFunctions")

x_end = 0.5
x_res = 0.0005
x = -x_end:x_res:x_end;

## Lorentzian
# Lorentzian Function
Lorentzian(x, x₀, Γ, A) = ((1/π).*((0.5 .* Γ) ./ (((x .- x₀).^2) .+ (0.5 .* Γ).^2)));

# Parameters
c1 = 0.1;   # centre
Γ1 = 0.1;   # width
Lorentziancurve1 = normalize(Lorentzian(x, c1, Γ1));

c2 = 0.0;   # centre
Γ2 = 0.2;   # width
Lorentziancurve2 = normalize(Lorentzian(x, c2, Γ2));

c3 = -0.10;   # centre
Γ3 = 0.3;   # width
Lorentziancurve3 = normalize(Lorentzian(x, c3, Γ3));

#Plot
pL = plot(x, Lorentziancurve1, title="Lorentzian function", xlabel ="x", ylabel="L(x)", label="x₀=0.1, Γ=0.1", linewidth=2)
plot!(x, Lorentziancurve2, label="x₀=0, Γ=0.2", linewidth=2)
plot!(x, Lorentziancurve3, label="x₀=-0.1, Γ=0.3", size=(1200,700), linewidth=2)
savefig("Lorentzian function.png")


FTL(k,x₀,Γ) = exp.((-2.0 .*π .* im .* k .* x₀) .- (Γ.*π.*abs.(k)));

## Gaussian
Gaussian(x, a, b, c²) = a .* exp.(-((x .- b).^2)/(2.0 .* c²))

# Parameters
a1 = 1.0;   #height
b1 = 0.1;   # centre
c²1 = 0.01;   # width
Gaussiancurve1 = Gaussian(x, a1, b1, c²1);

a2 = 0.8;
b2 = 0.0;   # centre
c²2 = 0.02;   # width
Gaussiancurve2 = Gaussian(x, a2, b2, c²2);

a3 = 0.6
b3 = -0.10;   # centre
c²3 = 0.03;   # width
Gaussiancurve3 = Gaussian(x, a3, b3, c²3);

plot(x, Gaussiancurve1, title="Gaussian", xlabel ="x", ylabel="G(x)", size=(1200, 700), label="a = 1, b=0.1, c²=0.01", linewidth=2)
plot!(x, Gaussiancurve2, label="a = 0.8, b=0, c²=0.02", linewidth=2)
plot!(x, Gaussiancurve3, label="a = 1, b=-0.1, c²=0.03", linewidth=2)

savefig("GaussianCurve.png")
