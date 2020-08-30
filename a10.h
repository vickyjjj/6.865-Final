#ifndef A10_H_PHUDVTKB
#define A10_H_PHUDVTKB

#include "Image.h"
#include "pixel.h"
#include "filtering.h"
#include "basicImageManipulation.h"

// TEXTURE SYNTHESIS FUNCTIONS
Image createNewImage(int new_width, int new_height, Image input, bool is_hole);
std::vector<Pixel> getUnfilledNeighbors(Image input);
Image getNeighborhoodWindow(Pixel p, Image input, int window_size);
Image getValidMask(Image input);
Pixel findMatch(Image temp, Image input, int window_size, Pixel p, bool is_hole);
Image textureSynthesis(Image input, int new_width, int new_height, int window_size, bool is_hole);

// INPAINTING FUNCTIONS
Image invert(Image input, int niter);
Image computeTensor(Image im, Image mask, float sigmaG = 1, float factorSigma = 0.5);
vector<Image> createImageAndMasks(Image input);
Image createFlatImage(Image input);
float dotProdIm(Image im1, Image im2);
Image poisson(Image bg, Image fg, Image mask, Image mask_inv, int niter = 200);
Image windowDiffsToCenter(Image input, Image tensor_original, Image tensor_guide, int x, int y, int window_size = 3);
Image diffusionCoefficients(Image window, float constant = 5.0f);
Image anisotropicDiffusion(Image input, Image tensor_original, Image tensor_guide, float step_size = 2.0f);
Image inpaint(Image im, int niter = 200);


#endif /* end of include guard: A10_H_PHUDVTKB */

