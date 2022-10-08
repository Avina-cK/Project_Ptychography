using FFTW, Plots, LinearAlgebra, ColorSchemes, CSV, Tables;

# 1 D object
#resolution of HD object : res_obj
res_obj = 256;
xfin = 5;

# Function for object
function re_test_object(x,x_fin)
    xdiv = (x_fin/x);
    xvec = xdiv:xdiv:x_fin;
    ans_reobj = sin.(xvec.^2) .* cos.(6 .*xvec) .* sin.(xvec);
    return ans_reobj
end # function for real part of test object

function im_test_object(x,x_fin)
    xdiv = (x_fin/x);
    xvec = xdiv:xdiv:x_fin;
    ans_imobj = sin.(2 .*xvec) .* cos.(xvec.^2) .* sin.(xvec);
end # function for imaginary part of test object

# real part of test object
re_obj = re_test_object(res_obj, xfin);
# imaginary part of test object
im_obj = im_test_object(res_obj,xfin);
# complete test object
object_i = re_obj + (im_obj)*1im;

# test object into a sparse diagonal matrix
objectM = Diagonal(object_i);

obj_Intensity = abs.(objectM);
obj_Amplitute = sqrt.(obj_Intensity);
obj_Amplitute_1d = diag(obj_Amplitute);
plot(obj_Amplitute_1d, legend=false, title="Input Object", ylabel="Amplitude", xlabel="Index")

CSV.write("C://.../TestObjectAmplitude.csv", Tables.table(obj_Amplitute_1d'), writeheader=false)
