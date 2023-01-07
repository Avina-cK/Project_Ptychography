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
cd(".../04_RodenburgFaulkner")

include(".../02_ForwardModel/ForwardModel.jl");
#= imports:    
    ω: freq domain
    Δω: resolution of ω
    t: time domain
    Δt: resolution of t

    object: object ∈ ℂ
    obj_amp: object Amplitude
    obj_phase: object phase
    Probe_set : set of dimension 1:2001 × -M_max:M_max
    Int_Measured : Set of internsity measurements 2001 × M_max
=#

# 2001×M_max OffsetArray with indices 1:2001 × -M:M
Int_Measured = OffsetArray(Int_Measured, 1:length(ω), -m_end:m_end);

#=
Algorithm
    1. Start with a guess of the object function 𝑂ᵧ,ₙ (𝐫). γ,n represents a 
        guessed function at the 𝑛th iteration of the algorithm. 𝑂 ∈ ℝ.   
    2. Multiply current guess at the object function by the illumination function
        at the current position 𝐑, 𝑃(𝐫-𝐑). This produces the guessed exit wave 
        function for position 𝐑,
            ψᵧ,ₙ(𝐫,𝐑) = 𝑂ᵧ,ₙ(𝐫).𝑃(𝐫-𝐑)
    3. Transform ψᵧ,ₙ(𝐫,𝐑) to obtain the corresponding wave function in the 
        diffraction space plane, for that position 𝐑,
            Ψᵧ,ₙ(𝐤,𝐑) = ℱ[ψᵧ,ₙ(𝐫,𝐑)]
        𝐤 is the usual reciprocal space coordinate. Ψ is a "guessed" version of the
        actual wave function in diffraction space. 
        Decompose into Amplitude and Phase:
            Ψᵧ,ₙ(𝐤,𝐑) = |Ψᵧ,ₙ(𝐤,𝐑)|exp(iθᵧ,ₙ(𝐤,𝐑))    
    4. Correct the intensities of the guessed diffraction space wave function to 
        known values,
            Ψc,n(𝐤,𝐑) = |Ψ(𝐤,𝐑)|exp(iθᵧ,ₙ(𝐤,𝐑))
        where |Ψ(𝐤,𝐑)| is known modulus.
    5. Inverse transform back to real space to obtain the new and improved guess at 
        the exit wave function
            ψc,n(𝐫,𝐑) = ℱ ⁻¹[Ψc,n(𝐤,𝐑)]
    6. Update the guessed object wave function in the area covered by the aperture 
        or probe, using the update function
            𝑂ᵧ,ₙ₊₁(𝐫) = 𝑂ᵧ,ₙ(𝐫) + (|P(𝐫-𝐑)|/|Pₘₐₓ(𝐫-𝐑)|)(P*(𝐫-𝐑)/(|P(𝐫-𝐑)|²+α)) 
                        × β(ψc,n(𝐫,𝐑)- ψᵧ,ₙ(𝐤,𝐑))
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
phase_0 = ω;
N = length(ω);  # number of pixels in the array
probe_offset = -2; #Offset of Probe
P_1 = Probe_set[:,probe_offset];

#= Step 1:
    Start with an initial guess of the object function. 𝑂_1 ∈ ℝ?  
=#
O_1 = O_i_amp;  

#= Step 2:
    Produce the guessed exit wave function for position 𝐑,
        ψᵧ,ₙ(𝐫,𝐑) = 𝑂ᵧ,ₙ(𝐫).𝑃(𝐫-𝐑)
=#
ψ_1 = O_1.*P_1;

#= Step 3:
    Transform ψᵧ,ₙ(𝐫,𝐑) to get the wave function in the diffraction space plane
        Ψᵧ,ₙ(𝐤,𝐑) = ℱ[ψᵧ,ₙ(𝐫,𝐑)]
    Ψ is a "guessed" version of the actual wave function in diffraction space. 
    Decompose into Amplitude and Phase  
=#

Ψ_γ_1 = fftshift(fft(ψ_1));   # guessed version of wave function
abs_Ψ_γ_1 = abs.(Ψ_γ_1);        # Amplitude of guess
phase_Ψ_γ_1 = angle.(Ψ_γ_1);    # phase of guess

#=
4. Correct the intensities of the guessed diffraction space wave function to 
    known values,
        Ψc,n(𝐤,𝐑) = |Ψ(𝐤,𝐑)|exp(iθᵧ,ₙ(𝐤,𝐑))
    where |Ψ(𝐤,𝐑)| is known modulus.
=#
Ψ_c_1 = Int_Measured[:,probe_offset].*exp.(im.*phase_Ψ_γ_1);

#= Step 5
 Inverse transform back to real space to obtain the new and improved guess at 
the exit wave function
    ψc,n(𝐫,𝐑) = ℱ ⁻¹[Ψc,n(𝐤,𝐑)]
=#
ψ_c_1 = ifft(ifftshift(Ψ_c_1));

#= Step 6
6. Update the guessed object wave function in the area covered by the aperture 
or probe, using the update function
    𝑂ᵧ,ₙ₊₁(𝐫) = 𝑂ᵧ,ₙ(𝐫) + (|P(𝐫-𝐑)|/|Pₘₐₓ(𝐫-𝐑)|)(P*(𝐫-𝐑)/(|P(𝐫-𝐑)|²+α)) 
                × β(ψc,n(𝐫,𝐑)- ψᵧ,ₙ(𝐤,𝐑))
=#
# update function
function O_new(O_old, P, α, β, ψ_c, ψ_g)
    P = vec(P);
    O_old = vec(O_old);
    ψ_c = vec(ψ_c);
    ψ_g = vec(ψ_g);
    return O_old .+ ((norm(P)./abs(maximum(abs.(P)))).*(vec(conj.(P))./((norm(P))^2 .+ α))) .* (β.*(ψ_c-ψ_g))
end

O_2 = O_new(O_1, P_1, 0.0, 0.1, ψ_c_1, Ψ_γ_1);

