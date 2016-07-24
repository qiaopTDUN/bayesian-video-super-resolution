function [y] = conv2d(I, kernel)
    fk = psf2otf(kernel, size(I));
    fI = fft2(I);
    fy = fk .* fI;
    y  = ifft2(fy);
end