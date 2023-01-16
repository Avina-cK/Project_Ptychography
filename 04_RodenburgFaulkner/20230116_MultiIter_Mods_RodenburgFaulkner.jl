using FFTW, LinearAlgebra, Plots, ColorSchemes, StatsBase

## Reference:
#=
    A Phase retrival algorithm for shifting illumination
        J. M. Rodenburg and H. M. L. Faulkner
        15 November 2004
        Applied Physics Letters
        Vol. 85
        No. 20
        https://aip.scitation.org/doi/10.1063/1.1823034
=#

# Setting up system and directory
cd("C://Users/avina/Documents/2020-2022_MSc_MathMods/401_Thesis/JuliaCode/04_RodenburgFaulkner")

include("C://Users/avina/Documents/2020-2022_MSc_MathMods/401_Thesis/JuliaCode/02_ForwardModel/ForwardModel.jl");
#= imports:    
    ω: freq domain
    Δω: resolution of ω
    t: time domain
    Δt: resolution of t

    object: object ∈ ℂ
    obj_amp: object Amplitude
    obj_phase: object phase
    Probe_set : set of dimension 1:2001 × -m_end:m_end
    Int_Measured : Set of internsity measurements 1:2001 × -m_end:m_end
=#

#=
Algorithm
   7. Move to the next position of 𝐑, for which the illumination in part overlaps
        with that of the previous position.
    8. Repeat (2) - (7) until the sum squared error (SSE) is sufficiently small. 
        The SSE is measured in the diffraction plane as
            SSE = (|Ψ(𝐤,𝐑)|²-|Ψᵧ,ₙ(𝐤,𝐑)|²)²/N
        where N is the number of pixels in the array representing the wave function.  
=#

# Initial vectors and variables
O_i_amp = peak₋₂ .+ peak₋₁ .+ peak₂; #initial guess of object
O_i_amp = normalize(O_i_amp);
a = 1.0;    # height of peak
b = 0.04;    # centre of peak
c² = 0.01;   # c = standard deviation : controls width

phase_i = normalize(Gaussian(ω, a, b.+(m*Δω), c²));

N = length(ω);  # number of pixels in the array

#= Step 1:
    Start with an initial guess of the object function. 𝑂_1 ∈ ℝ?  
=#
O_1 = O_i_amp.*exp.(im.*phase_i);  

probe_offset_0 = -m_end; #Offset of Probe
probe_offset_end = m_end;

# create an object matrix to store updated objects
Object = zeros(N, M_max) + im .*zeros(N, M_max);
Object = OffsetArray(Object,1:N, -m_end: m_end);
Object[:, probe_offset_0] = O_1;

# update function for object
function O_new(O_old, P, α, β, ψ_c, ψ_g)
    P = vec(P);
    O_old = vec(O_old);
    ψ_c = vec(ψ_c);
    ψ_g = vec(ψ_g);
    a = abs.(P)/(maximum(abs.(P)));
    b = (1/(maximum(abs.(P))+α)).*conj.(P);
    c = β.*(ψ_c-ψ_g);
    return O_old + a.*b.*c
end

# Sum squared error
function SSqE(Ψ, Ψ_g, len_N)
    Ψ = vec(Ψ);
    Ψ_g = vec(Ψ_g);
    return (norm(Ψ)^2 - norm(Ψ_g)^2)^2/(len_N);   
end
# vector storing SSE
SSqE_it = zeros(size(Probe_set)[2],2);
SSqE_it = OffsetArray(SSqE_it,-6:6,0:1);

i= probe_offset_0;

# setting parameters
alpha = 0.0001;
beta = 01.;
ϵ = 2.1e-7;
i_end = probe_offset_end;

while i<=probe_offset_end
    probe_offset = i;
    # current probe 
    Probe_i = Probe_set[:,probe_offset];
    # Associated sqrt of intensity reading
    mod_read_i = sqrt.(Int_Measured[:, probe_offset]);
    # Step 1: Current guess for object
    O_i = Object[:,i];
    #= Step 2:
    Produce the guessed exit wave function for position 𝐑,
        ψᵧ,ₙ(𝐫,𝐑) = 𝑂ᵧ,ₙ(𝐫).𝑃(𝐫-𝐑) =#
    ψ_i = O_i .* Probe_i;
    #= Step 3:
    Transform ψᵧ,ₙ(𝐫,𝐑) to get the wave function in the diffraction space plane
        Ψᵧ,ₙ(𝐤,𝐑) = ℱ[ψᵧ,ₙ(𝐫,𝐑)]
    Ψ is a "guessed" version of the actual wave function in diffraction space. =#
    Ψ_γ_i = fftshift(fft(ψ_i));
    # Decompose into Amplitude and Phase 
    abs_Ψ_γ_i = abs.(Ψ_γ_i);        # Amplitude of guess
    phase_Ψ_γ_i = angle.(Ψ_γ_i);    # phase of guess
    #= Step 4:
    Correct the intensities of the guessed wave function to known values,
        Ψc,n(𝐤,𝐑) = |Ψ(𝐤,𝐑)|exp(iθᵧ,ₙ(𝐤,𝐑))
    where |Ψ(𝐤,𝐑)| is known modulus. =#
    Ψ_c_i = mod_read_i.*exp.(im*phase_Ψ_γ_i);
    #= Step 5
    Inverse transform back to real space to obtain the new and improved guess at 
    the exit wave function
        ψc,n(𝐫,𝐑) = ℱ ⁻¹[Ψc,n(𝐤,𝐑)] =#
    ψ_c_i = ifftshift(ifft(Ψ_c_i));
    #= Step 6
    Update the guessed object wave function in the area covered by the probe, using the update function
        𝑂ᵧ,ₙ₊₁(𝐫) = 𝑂ᵧ,ₙ(𝐫) + (|P(𝐫-𝐑)|/|Pₘₐₓ(𝐫-𝐑)|)(P*(𝐫-𝐑)/(|P(𝐫-𝐑)|²+α)) 
                    × β(ψc,n(𝐫,𝐑)- ψᵧ,ₙ(𝐤,𝐑))
    =#
    O_iplus1 = O_new(O_i, Probe_i, alpha, beta, ψ_c_i, ψ_i);
    SSE = SSqE(mod_read_i, Ψ_γ_i, N);
    if i<probe_offset_end
        global Object[:,i+1] = O_iplus1;
        global SSqE_it[i,1] = SSE;
    end
    if SSE<ϵ
        global i_end = i;
        global i = probe_offset_end+10;
    else
        global i=i+1;
    end
end
print("Loop Over")
plot((angle.(Object[:,i_end]))./maximum(angle.(Object[:,i_end])), label="Object final angle")
plot!((angle.(object))./maximum(angle.(object)), label="Object true angle")
plot!((angle.(Object[:,-6])), label = "Object initial angle guess")
SSqE_it = vec(SSqE_it[:,1]);
plot(SSqE_it[probe_offset_0:i_end-1], title="SSE") 
