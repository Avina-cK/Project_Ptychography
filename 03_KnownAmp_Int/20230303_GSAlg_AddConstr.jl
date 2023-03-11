using FFTW, LinearAlgebra, Plots, ColorSchemes, StatsBase
using Plots.PlotMeasures

## Reference:
#=

    A practical algorithm for the determination
    of phase from image and diffraction plane pictures
        R. W. Gerchberg and W. O. Saxton,
        1972
        Optik
        Vol. 35 
        pp. 237–247.

and

    Phase retrival algorithms: a comparison
        J.R. Fienup
        1 Aug 1982
        Applied Optics
        Vol. 21
        No. 15
=#

## Setting up system and directory
include(".../01_SettingUpTestObject_Probe/CreatingTestObject_and_Probe.jl")

cd(".../03_KnownAmp-Int")


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

# auto-correlation function
ift_fft_object = ifft(ifftshift(fft_object));

# diameter of real and imaginary parts of the auto-correlation function
max_real = maximum(real.(ift_fft_object));  # maximum of real part
min_real = minimum(real.(ift_fft_object));  # minimum of real part
max_imag = maximum(imag.(ift_fft_object));  # maximum of imaginary part
min_imag = minimum(imag.(ift_fft_object));  # minimum of imaginary part
diam_real = max_real - min_real;    # diameter of real part
diam_imag = max_imag - min_imag ;   # diameter of imaginary part

abs_F = sqrt.(abs.(fft_object));

plot(t, abs_F
    , size = (1450,600)                        
    , ylabel = "I(x)", xlabel = "x"
    , title = "Intensity measurement"
    , legend = false
    , color =:blue
    , titlefontsize = 15, guidefontsize = 12
    )
#savefig("IntensityMeasurement_GBAlg.png")

plot(ω ,abs_f
    , size = (1450,600)  
    , ylabel = "|O(ω)|", xlabel = "ω"
    , title = "Object magnitude"
    , legend = false
    , color =:green
    , titlefontsize = 15, guidefontsize = 12
    )
#savefig("ObjectMagnitude_GBAlg.png")

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
i_f = 5000;

phase_int = zeros(N,i_f);   # matrix to store changes in phase
phase_int[:,1] = phase_0;

# plotting initial phase guess
plot(ω , phase_0
    , size = (1450,600)  
    , ylabel = "Phase(O(g,1)(ω))", xlabel = "ω"
    , title = "Initial guess of object phase"
    , legend = false
    , color =:darkgreen
    , titlefontsize = 15, guidefontsize = 12
    )
# savefig("InitialObjGuessPhase_GSAlg.png")


obj_int = zeros(N,i_f) + im.*zeros(N,i_f);  # matrix to store changes in object
obj_int[:,1] = abs_obj_0.*exp.(im.*phase_int[:,1]); # initial object guess

B_it = zeros(i_f-1) ;   # vector to store Fourier-domain error
E_it = zeros(i_f-1) ;   # vector to store Object-domain error

#rand_range = -1.0 : 0.00000001 : 1.0;

# Tolerance for condition
ϵ = 0.05;

abs_obj = obj_int[:,1];

# maximum number of iterations
i_end = i_f;

while i<i_f
    # current object guess
    est_obj = obj_int[:,i];

    # step 1: FT the estimate of the object
    ft_estobj = fftshift(fft(est_obj));
    ft_estobj_phase = angle.(ft_estobj);

    # Does it satisfy the Fourier constraints?
    global B_it[i] = B(ft_estobj,fft_object);
    if B_it[i]<ϵ
        println("B<ϵ")
    # else println("B>=ϵ ", B_it[i])
    end

    # step 2: Replace the modulus of the computed FT with the measure FT
    est_ft_estobj = abs_F.*exp.(im.*ft_estobj_phase);
    
    # step 3: inverse FT the estimate
    ifft_estft = ifftshift(ifft(est_ft_estobj));
    ifft_estft_phase = angle.(ifft_estft)

    # Does it satisfy the object function constraints?
    global E_it[i] = E(object,ifft_estft);  #update current function constraint
    # check if the object function constraint is satisfied
    if E_it[i]<ϵ
        println("E<ϵ")
        # global i_end = i;
        global i = i_f+10;  
    else #println("E>=ϵ ", E_it[i])
    end

    # step 4: replace modulus of computed image with measured obj modulus
    new_est_obj = abs_obj.*exp.(im.*ifft_estft_phase);
    real_new_est_obj = real.(new_est_obj);
    imag_new_est_obj = imag.(new_est_obj);
    
    # create a vector to store constraint
    indices = zeros(N); 
    # Constraint: is the guessed object within the estimated diameter
    for j in 1:N
        if real_new_est_obj[j]<=max_real && imag_new_est_obj[j]<=max_imag
            if real_new_est_obj[j]>=min_real && imag_new_est_obj[j]>=min_imag
                indices[j] = 1;
            end
        end
    end    

    global obj_int[:,i+1] = new_est_obj.*indices;
    
    if mod(i,100)==0
        print(i," ")
    end
    
    global i = i+1;

end

plot(B_it, title="Fourier domain constraint")

E_mins_index = findall(x->x<=0.85, E_it); 
plot_ind_i= collect(1:20:length(E_it));
plot_ind = sort([plot_ind_i; E_mins_index]);

plot(size = (1450,600))
plot!(plot_ind, E_it[plot_ind]
    , title="Object domain constraint", label=false
    , alpha = 0.9, marker = true
    , xlabel = "Iteration [n]", ylabel = "SSE₀(n)"
    , color =:purple
    , titlefontsize = 15, guidefontsize = 12
    , bottom_margin = 5mm
    , left_margin = 6mm
    ) 
# savefig("ObjDomainConstraint_GS.png")

E_min_index = findall(x->x==minimum(E_it), E_it); 
plot(ω, angle.(obj_int[:,E_min_index])
    , label = "Phase of object with min(SSE₀)"
    , xlabel = "ω"
    , ylabel = "Phase(O(ω))"
    , size = (1450,600)
    )
plot!(ω, obj_phase
    , label = "True phase")
