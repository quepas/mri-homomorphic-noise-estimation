from scipy import ndimage, special, fftpack
import numpy
import time
from numpy import exp, genfromtxt, array, real

__author__ = 'Quepas'

def filter2B(image, mask):
    (Mx, My) = mask.shape
    # check if Mrow and Mcol are odd
    if Mx % 2 == 0 or My % 2 == 0:
        print "Mask width and height must be odd"
        return array([])
    Nx = (Mx - 1) / 2
    Ny = (My - 1) / 2
    filtered_image = ndimage.filters.convolve(image, mask, mode='nearest')
    return filtered_image

def approxI1_I0(image):
    count = numpy.sum(image < 1.5)
    image_oct = 8 * image
    Mn = 1. - 3. / image_oct - 15. / 2. / (image_oct ** 2) - (3. * 5. * 21.) / 6. / (image_oct ** 3)
    Md = 1. + 1. / image_oct + 9. / 2. / (image_oct ** 2) + (25. * 9.) / 6. / (image_oct ** 3)
    M = Mn / Md
    if (count > 1):
        K = numpy.flatnonzero(image < 1.5)
        M.ravel()[K] = special.iv(1, image.ravel()[K]) / special.iv(0, image.ravel()[K])
    K = numpy.flatnonzero(image == 0)
    M.ravel()[K] = 0
    return M

def em_ml_rice2D(rice_image, iterations, window_size):
    rice_image = rice_image.astype(numpy.double)
    mask = numpy.ones(window_size) / numpy.prod(window_size)

    diff_left = 2 * filter2B(rice_image ** 2, mask) ** 2
    diff_right = filter2B(rice_image ** 4, mask)

    a_k = numpy.sqrt(numpy.sqrt(numpy.maximum(diff_left - diff_right, 0)))
    sigma_k2 = 0.5 * numpy.maximum(filter2B(rice_image ** 2, mask) - a_k ** 2, 0.01)

    for i in range(iterations):
        approx_image = approxI1_I0(a_k * rice_image / sigma_k2) * rice_image
        a_k = numpy.maximum(filter2B(approx_image, mask), 0)
        image_abs = abs(rice_image) ** 2
        sigma_k2 = numpy.maximum(0.5 * filter2B(image_abs, mask) - a_k ** 2 / 2.0, 0.01)

    return (a_k, numpy.sqrt(sigma_k2))

def matlab_style_gauss2D(shape=(3, 3), sigma=0.5):
    """
    2D gaussian mask - should give the same result as MATLAB's
    fspecial('gaussian',[shape],[sigma])
    """
    m, n = [(ss - 1.) / 2. for ss in shape]
    y, x = numpy.ogrid[-m:m + 1, -n:n + 1]
    h = numpy.exp(-(x * x + y * y) / (2. * sigma * sigma))
    h[h < numpy.finfo(h.dtype).eps * h.max()] = 0
    sumh = h.sum()
    if sumh != 0:
        h /= sumh
    return h

def lpf(image, sigma, mode=2):
    (Mx, My) = image.shape
    if mode == 1:
        kernel = matlab_style_gauss2D(image.shape, sigma)
        kernel /= numpy.max(kernel)

        if Mx == 1 or My == 1:
            fft = numpy.fft.fft(image)
            fft = numpy.fft.fftshift(fft)
            fft *= kernel
            result = numpy.real(numpy.fft.ifft(numpy.fft.fftshift(fft)))
        else:
            fft = numpy.fft.fftshift(numpy.fft.fft2(image))
            fft *= kernel
            result = numpy.real(numpy.fft.ifft2(numpy.fft.fftshift(fft)))
    elif mode == 2:
        new_dim = 2 * array(image.shape)
        kernel = matlab_style_gauss2D((new_dim[0], new_dim[1]), sigma * 2)
        kernel /= numpy.max(kernel)
        kernel = kernel[Mx:, My:]

        image = image.astype(numpy.double)
        if Mx == 1 or My == 1:
            dct = fftpack.dct(image, type=1)
            dct *= kernel
            result = numpy.real(fftpack.idct(dct, type=1))
        else:
            dct = fftpack.dct(fftpack.dct(image.T, type=2, norm='ortho').T, type=2, norm='ortho')
            dct *= kernel
            result = numpy.real(fftpack.idct(fftpack.idct(dct.T, type=2, norm='ortho').T, type=2, norm='ortho'))
    return result

def correct_rice_gauss(SNR):
    c = [-0.289549906258443,  -0.038892257560633,    0.409867108141953,
         -0.355237628488567,  0.149328280945610,     -0.0357861117942093,
         0.00497952893859122, -0.000374756374477592, 1.18020229140092e-05]
    rice_gauss = c[0] + \
                 c[1] * SNR + \
                 c[2] * (SNR ** 2) + \
                 c[3] * (SNR ** 3) + \
                 c[4] * (SNR ** 4) + \
                 c[5] * (SNR ** 5) + \
                 c[6] * (SNR ** 6) + \
                 c[7] * (SNR ** 7) + \
                 c[8] * (SNR ** 8)
    return rice_gauss * (SNR <= 7)

def rice_homomorf_est(image, SNR = 0, LPF = 4.8, mode = 2):
    (M2, Sigma_n) = em_ml_rice2D(image, 10, [3, 3]);
    Sigma_n2 = lpf(Sigma_n, 1.2)
    M1 = filter2B(image, numpy.ones((5, 5)) / 25)

    if (SNR.shape[0] == 1) and SNR == 0:
        SNR = M2 / Sigma_n

    Rn = abs(image - M1)
    lRn = numpy.log(Rn * (Rn != 0) + 0.001 * (Rn == 0))
    LPF2 = lpf(lRn, LPF)
    Mapa2 = numpy.exp(LPF2)
    MapaG = Mapa2 * 2 / numpy.sqrt(2) * numpy.exp(-special.psi(1)/2.)

    LocalMean = 0
    if mode == 1:
        LocalMean = M1
    elif mode == 2:
        LocalMean = M2

    Rn = numpy.abs(image - LocalMean)
    lRn = numpy.log(Rn * (Rn != 0) + 0.001 * (Rn == 0))
    LPF2 = lpf(lRn, LPF)
    Fc1 = correct_rice_gauss(SNR)
    LPF1 = LPF2 - Fc1
    LPF1 = lpf(LPF1, LPF + 2, 2.)
    Mapa1 = exp(LPF1)
    MapaR = Mapa1*2/numpy.sqrt(2)*numpy.exp(-special.psi(1)/2.)
    return (MapaR, MapaG)

def run_example():
    MR_SNR = genfromtxt('MR_SNR.csv', delimiter=',')
    MR_noisy = genfromtxt('MR_noisy.csv', delimiter=',')
    # estymacja EM przy znanym SNR
    now = time.clock()
    (MapaR_EM_SNR, MapaG_EM_SNR) = rice_homomorf_est(MR_noisy, MR_SNR, 3.4, 2);
    numpy.savetxt("MapaR_EM_SNR.csv", MapaR_EM_SNR, delimiter=',')
    numpy.savetxt("MapaG_EM_SNR.csv", MapaG_EM_SNR, delimiter=',')
    print "EM_SNR time : " + str(time.clock() - now) + "s"

    # estymacja EM przy nieznanym SNR
    now = time.clock()
    (MapaR_EM, MapaG_EM) = rice_homomorf_est(MR_noisy, array([0]), 3.4, 2);
    numpy.savetxt("MapaR_EM.csv", MapaR_EM, delimiter=',')
    numpy.savetxt("MapaG_EM.csv", MapaG_EM, delimiter=',')
    print "EM time : " + str(time.clock() - now) + "s"

    # estymacja local mean przy znanym SNR
    now = time.clock()
    (MapaR_LM_SNR, MapaG_LM_SNR) = rice_homomorf_est(MR_noisy, MR_SNR, 3.4, 1);
    numpy.savetxt("MapaR_LM_SNR.csv", MapaR_LM_SNR, delimiter=',')
    numpy.savetxt("MapaG_LM_SNR.csv", MapaG_LM_SNR, delimiter=',')
    print "LM SNR time : " + str(time.clock() - now) + "s"

    # estymacja local mean przy nieznanym SNR
    now = time.clock()
    (MapaR_LM, MapaG_LM) = rice_homomorf_est(MR_noisy, array([0]), 3.4, 1);
    numpy.savetxt("MapaR_LM.csv", MapaR_LM, delimiter=',')
    numpy.savetxt("MapaG_LM.csv", MapaG_LM, delimiter=',')
    print "LM time : " + str(time.clock() - now) + "s"

run_example()
