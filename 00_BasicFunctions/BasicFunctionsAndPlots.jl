using LinearAlgebra, Plots, FFTW

cd(".../00_BasicFunctions")

x_end = 0.5
x_res = 0.0005
x = -x_end:x_res:x_end;

"""
    Lorentzian(x_vec, x₀, Γ)

Defining Lorentzian curve with:
    x_vec ∈ ℝⁿ is the domain
    x₀ ∈ ℝ is the centre of the curve
    Γ ∈ ℝ is the width of the curve 
"""
function Lorentzian(x_vec, x₀, Γ)
    ((1/π).*((0.5 .* Γ) ./ (((x_vec .- x₀).^2) .+ (0.5 .* Γ).^2)))
end

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
#savefig("Lorentzian function.png")

"""
    FTL(k,x₀,Γ)

Fourier transform of the Lorentzian?
"""
function FTL(k,x₀,Γ)
    exp.((-2.0 .*π .* im .* k .* x₀) .- (Γ.*π.*abs.(k)))
end

"""
    Gaussian(x, a, b, c²)
Create a Gaussian curve with the following inputs:
    x ∈ ℝⁿ is the domain
    a ∈ ℝ is the height of the peak's curve
    b ∈ ℝ is the position of the central peak
    c ∈ ℝ⁺ is the standard deviation
"""
function Gaussian(x, a, b, c²)
    a .* exp.(-((x .- b).^2)/(2.0 .* c²))
end
a1 = 1.0;   #height
b1 = 0.1;   # centre
c²1 = 0.01;   # width
Gaussiancurve1 = Gaussian(x, a1, b1, c²1);

a2 = 0.8;
b2 = 0.0;   # centre
c²2 = 0.02;   # width
Gaussiancurve2 = Gaussian(x, a2, b2, c²2);

a3 = 0.6;
b3 = -0.10;   # centre
c²3 = 0.03;   # width
Gaussiancurve3 = Gaussian(x, a3, b3, c²3);

plot(x, Gaussiancurve1, title="Gaussian", xlabel ="x", ylabel="G(x)", size=(1200, 700), label="a = 1, b=0.1, c²=0.01", linewidth=2)
plot!(x, Gaussiancurve2, label="a = 0.8, b=0, c²=0.02", linewidth=2)
plot!(x, Gaussiancurve3, label="a = 1, b=-0.1, c²=0.03", linewidth=2)
#savefig("GaussianCurve.png")

"""
    sinc probe    

"""
function sinc_probe(x, a)
    abs.(sinc.(a.*x))
end

plot(title="|sinc(a.x)|", legendtitle="a=")
plot!(sinc_probe(x,10), label="10")
plot!(sinc_probe(x,20), label="20")
plot!(sinc_probe(x,30), label="30")
