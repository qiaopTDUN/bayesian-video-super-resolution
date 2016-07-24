function [y] = conv2dt(I, kernel)
    fk = conj(psf2otf(kernel, size(I)));
    fI = fft2(I);
    fy = fk .* fI;
    y  = ifft2(fy);
end