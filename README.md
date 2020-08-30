# 6.865 Final: Texture Synthesis and Image Inpainting

I compared the motivation and results of texture synthesis and image inpainting for my 6.865 Computational Photography final project. For texture synthesis and hole-filling, I coded the method in "Texture Synthesis by Non-Parametric Sampling" (Efros and Leung, 1999) over images of both stochastic and tiled patterns, by implementing a per-pixel window sampling algorithm Gaussian blur and a normalized sum of squared differences. For image inpainting, I implemented the iterative method in "Image Inpainting" (Bertalmio, Sapiro, Caselles and Ballester, 2000), which aims to connect isophite lines at region boundaries using the structure tensor calculated from the luminance gradients and Poisson-based gradient descent, and then estimates coloration using anisotropic diffusion.

Note: Some files have been redacted as they are part of class files. 

## Texture Synthesis Methods

### Image createNewImage(int new_width, int new_height, Image input, bool is_hole);

### std::vector<Pixel> getUnfilledNeighbors(Image input);

### Image getNeighborhoodWindow(Pixel p, Image input, int window_size);

### Image getValidMask(Image input);

### Pixel findMatch(Image temp, Image input, int window_size, Pixel p, bool is_hole);

### Image textureSynthesis(Image input, int new_width, int new_height, int window_size, bool is_hole);

## Image Inpainting Methods

### Image invert(Image input, int niter);

### Image computeTensor(Image im, Image mask, float sigmaG = 1, float factorSigma = 0.5);

### vector<Image> createImageAndMasks(Image input);

### Image createFlatImage(Image input);

### float dotProdIm(Image im1, Image im2);

### Image poisson(Image bg, Image fg, Image mask, Image mask_inv, int niter = 200);

### Image windowDiffsToCenter(Image input, Image tensor_original, Image tensor_guide, int x, int y, int window_size = 3);

### Image diffusionCoefficients(Image window, float constant = 5.0f);

### Image anisotropicDiffusion(Image input, Image tensor_original, Image tensor_guide, float step_size = 2.0f);

### Image inpaint(Image im, int niter = 200);
