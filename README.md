# 6.865 Final: Texture Synthesis and Image Inpainting

I compared the motivation and results of texture synthesis and image inpainting for my 6.865 Computational Photography final project. For texture synthesis and hole-filling, I coded the method in "Texture Synthesis by Non-Parametric Sampling" (Efros and Leung, 1999) over images of both stochastic and tiled patterns, by implementing a per-pixel window sampling algorithm Gaussian blur and a normalized sum of squared differences. For image inpainting, I implemented the iterative method in "Image Inpainting" (Bertalmio, Sapiro, Caselles and Ballester, 2000), which aims to connect isophite lines at region boundaries using the structure tensor calculated from the luminance gradients and Poisson-based gradient descent, and then estimates coloration using anisotropic diffusion.

Note: Some files have been redacted as they are part of class files. 

## Texture Synthesis Methods

### Image createNewImage(int new_width, int new_height, Image input, bool is_hole);
Given a new width and new height, return an image of that size with the input image centered within.

### Image adjustInput(Image im)
Increases all black pixels to EMPTY_THRESHOLD, since we treat black pixels as parts that need to be filled.

### std::vector\<Pixel\> getUnfilledNeighbors(Image input);
Returns all pixels that are not filled in input image with filled pixels as neighbors.

### Image getNeighborhoodWindow(Pixel p, Image input, int window_size);
Returns an image with the window of size window_size of neighboring pixels to p in input..

### Image getValidMask(Image input);
Returns 1's where image is filled, 0 otherwise.

### Pixel findMatch(Image temp, Image input, int window_size, Pixel p, bool is_hole);
Finds suitable pixel match and returns it. If no match found, returns Pixel with value (-1, -1).

### Image textureSynthesis(Image input, int new_width, int new_height, int window_size, bool is_hole);

## Image Inpainting Methods

### Image invert(Image input, int niter);
Inverts black and white image pixels.

### Image computeTensor(Image im, Image mask, float sigmaG = 1, float factorSigma = 0.5);

### vector\<Image\> createImageAndMasks(Image input);
Returns a 3-element vector with the new_image, where masked areas are empty (have value 0), the mask, where masked areas have value 1 and non-masked areas have value 0, and the mask inverse. 

### Image createFlatImage(Image input);
Returns a "flat" white image.

### float dotProdIm(Image im1, Image im2);
Returns the dot product of an image by multiplying two images element-wise, then summing the values. 

### Image poisson(Image bg, Image fg, Image mask, Image mask_inv, int niter = 200);

### Image windowDiffsToCenter(Image input, Image tensor_original, Image tensor_guide, int x, int y, int window_size = 3);
Calculates differences between pixels in window from center pixel value. Original and guide should both be structure tensors.

### Image diffusionCoefficients(Image window, float constant = 5.0f);
Calculates positional diffusion coefficients using formual e^(-(gradI/k)^2)

### Image anisotropicDiffusion(Image input, Image tensor_original, Image tensor_guide, float step_size = 2.0f);
Computes anistropic diffusion following instructions from http://www.cs.utah.edu/~manasi/coursework/cs7960/p2/project2.html.

### Image inpaint(Image im, int niter = 200);
