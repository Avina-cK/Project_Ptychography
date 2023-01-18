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

cd("C://Users/avina/Documents/2020-2022_MSc_MathMods/401_Thesis/JuliaCode/03_KnownAmp-Int")

include("C://Users/avina/Documents/2020-2022_MSc_MathMods/401_Thesis/JuliaCode/01_SettingUpTestObject_Probe/CreatingLorentzian_TestObject_and_Probe.jl")

#=
imports:    ω: freq domain
            Δω: resolution of ω
            t: time domain
            Δt: resolution of t

            object: object ∈ ℂ
            obj_amp: object Amplitude
            obj_phase: object phase
                        phase_gaus if gaussian phase
            Probe_set : set of dimension 1:2001 × -2:2
=#

## Setting up readings
# Measurements: abs_f = obj_amp and abs_F = |Fourier[f(x)]|

abs_f = obj_amp;

# Let us simplify the problem and make the max value of both 1
abs_f = (abs_f);
fft_object = fftshift(fft(object));
abs_F = sqrt.(abs.(fft_object));

plot(abs_F)
plot!(abs_f)

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

# Fourier domain constraint:
function B(G,F)
    N = length(G)
    b = (abs.(G).-abs.(F)).^2
    B = (1/√N)*sum(b)
    return B
end

# Function domain contraint:
function E(f,g)
    e = (abs.(f).-abs.(g)).^2
    E = sum(e)
    return E
end

# Create initial guess of object
phase_0 = ω;
abs_obj_0 = abs_f;
N = length(ω);

i = 1;
i_f = 500;

phase_int = zeros(N,i_f);
phase_int[:,1] = phase_0;

obj_int = zeros(N,i_f) + im.*zeros(N,i_f);
obj_int[:,1] = abs_obj_0.*exp.(im.*phase_int[:,1]);

B_it = zeros(i_f-1) ;
E_it = zeros(i_f-1) ;

rand_range = -1.0 : 0.00000001 : 1.0;

# Tolerance for condition
ϵ = 0.1;
abs_obj = obj_int[:,1];

while i<i_f
    est_obj = obj_int[:,i];
    # step 1: FT the estimate of the object
    ft_estobj = fftshift(fft(est_obj));
    ft_estobj_phase = angle.(ft_estobj);

    # Does it satisfy the Fourier constraints?
    global B_it[i] = B(ft_estobj,fft_object);
    if B_it[i]<ϵ
        println("B<ϵ")
    else #println("B>=ϵ ", B_it[i])
    end

    # step 2: Replace the modulus of the computed FT with the measure FT
    est_ft_estobj = abs_F.*exp.(im.*ft_estobj_phase);
    # step 3: inverse FT the estimate
    ifft_estft = ifftshift(ifft(est_ft_estobj));
    ifft_estft_phase = angle.(ifft_estft)

    # Does it satisfy the function constraints?
    global E_it[i] = E(object,ifft_estft);
    if E_it[i]<ϵ
        println("E<ϵ")
    else #println("B>=ϵ ", B_it[i])
    end

    # step 4: replace modulus of computed image with measured obj modulus
    new_est_obj = abs_obj.*exp.(im.*ifft_estft_phase);

    # Imposing constraints to update
    constraint1 = zeros(N);
    for j=1:N
        if angle.(new_est_obj)[j]<0.0
            constraint1[j] = -1.0
        else constraint1[j] = 1.0
        end
    end
    for k=1:Int8(floor(N/2))
        if constraint1[k] == -1.0
            global obj_int[k,i+1] = new_est_obj[k];
        else
            global obj_int[k,i+1] = obj_int[k,i];
        end
    end
    for k=(Int8(floor(N/2))+1):N
        if constraint1[k] == 1.0
            global obj_int[k,i+1] = new_est_obj[k];
        else
            global obj_int[k,i+1] = obj_int[k,i];
        end
    end
    # global obj_int[:,i+1] = new_est_obj;
    
    global i = i+1;
    
end
plot(B_it)
plot(E_it.-minimum(E_it))
