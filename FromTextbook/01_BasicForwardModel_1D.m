%% 1 D object
%plot grid
m_grid = 2; n_grid=3;
%resolution of HD object : res_obj
res_obj = 256;
xfin = 5;

% real part of test object
re_obj = re_test_object(res_obj, xfin);
% imaginary part of test object
im_obj = im_test_object(res_obj,xfin);
% complete test object
object_i = complex(re_obj,im_obj);

% test object into a diagonal matrix
objectM = diag(object_i);

obj_Intensity = abs(objectM);
obj_AmplituteM = sqrt(obj_Intensity);
obj_Amplitute = diag(obj_AmplituteM);
subplot(m_grid,n_grid,1);
plot(obj_Amplitute);
title('Input Object');


%%
% set up the coherent imaging system (parameters)
lambda = 0.5e-6;            % wavelength
k0 = 2*pi/lambda;    
PS = 0.5e-6;                % Pixel size of image sensor
NA = 0.5;                   % Numerical aperture of objective lens
CutOffFreq = NA*k0;

% simulate low pass filtering process of imaging system
FT_obj_Intensity = fftshift(fft2(obj_AmplituteM));
[m n] = size(obj_AmplituteM);
k_x = (-pi/PS):(2*pi/(PS*(n-1))):(pi/PS);
k_y = (-pi/PS):(2*pi/(PS*(n-1))):(pi/PS);
[k_xm k_ym] = meshgrid(k_x, k_y);
CohTransFunc = (k_ym.^2 + k_xm.^2)<CutOffFreq^2;
subplot(m_grid,n_grid,2);
imshow(CohTransFunc, []);
title('Coherent transfer function(CTF) in Fourier domain');

% set up incoherent transfer function
cohPSF = fftshift(ifft2(ifftshift(CohTransFunc)));
incohPSF = (abs(cohPSF)).^2;
incohTransFunc = abs(fftshift(fft2(ifftshift(incohPSF))));  %incoherent transfer function
incohTransFunc = incohTransFunc./max(max(incohTransFunc));
subplot(m_grid,n_grid,3);
imshow(abs(incohTransFunc), []);
title('Incoherent transfer function (ITF)');

% perform low-pass filtering and generate output itensity image
output_FT = incohTransFunc.*FT_obj_Intensity;
subplot(m_grid,n_grid,4);
imshow(log(abs(output_FT)),[]);
title('Output Image in Fourier domain (ITF)');
output_Int = ifft2(ifftshift(output_FT));
subplot(m_grid,n_grid,5);
plot(diag(output_Int));
title('Output Image (ITF)');

% error (|final - initial| image)
subplot(m_grid,n_grid, 6);
scatter(obj_Amplitute,diag(output_Int));
xlabel('Input'); 
ylabel('Output')
%% Functions
% real part of test object
function ans_reobj = re_test_object(x,x_fin)
    xdiv = (x_fin/x);
    xvec = xdiv:xdiv:x_fin;
    ans_reobj = sin(xvec.^2) .* cos(6.*xvec) .* sin(xvec);
end

% imaginary part of test object
function ans_imobj = im_test_object(x,x_fin)
    xdiv = (x_fin/x);
    xvec = xdiv:xdiv:x_fin;
    ans_imobj = sin(2 .*xvec) .* cos(xvec.^2) .* sin(xvec);
end
