#include "a10.h"

#define ERR_THRESHOLD 0.1
#define MAX_ERR_THRESHOLD 0.3
#define EMPTY_THRESHOLD 0.01

using namespace std;

// TEXTURE SYNTHESIS FUNCTIONS

// Given a new width and new height, return an image of that size
// with the input image centered within.
Image createNewImage(int new_width, int new_height, Image input, bool is_hole) {
	Image output(new_width, new_height, 1); 

	// For all pixels with value 0 and not hole filling, 
	// add slight error because we detect holes as purely black spaces 
	int i, j;
	for (i = 0; i < input.width(); i++) {
		for (j = 0; j < input.height(); j++) {
			if (!is_hole && input(i, j) == 0) {
				output(i + new_width/2 - input.width()/2, j + new_height/2 - input.height()/2) = input(i, j) + EMPTY_THRESHOLD;
			} else {
				output(i + new_width/2 - input.width()/2, j + new_height/2 - input.height()/2) = input(i, j);
			}
		}
	}

	return output;
}

// Increases all black pixels to EMPTY_THRESHOLD, since we treat black pixels as parts that need to be filled
Image adjustInput(Image im) {
	Image output(im.width(), im.height(), 1);

	int i, j;
	for (i = 0; i < im.width(); i++) {
		for (j = 0; j < im.height(); j++) {
			if (im(i, j) == 0) {
				output(i, j) = im(i, j) + EMPTY_THRESHOLD;
			} else {
				output(i, j) = im(i, j);
			}
		}
	}
	return output;
}

// Returns all pixels that are not filled in input image with filled pixels as neighbors
std::vector<Pixel> getUnfilledNeighbors(Image input) {
	// TODO permute and sort 
	std::vector<Pixel> output;

	Image dilated = dilate_filterClass(input, 3, true);

	// perform "subtraction" and get points that are diff
	int i, j;
	for (i = 0; i < input.width(); i++) {
		for (j = 0; j < input.height(); j++) {
			if (dilated(i, j) > 0) {
				if (input(i, j) == 0) {
					Pixel p(i, j);
					output.push_back(p);
				}
			}
		}
	}

	return output;
}

// Returns an image with the window of size window_size of neighboring pixels to p in input.
Image getNeighborhoodWindow(Pixel p, Image input, int window_size) {
	Image output(window_size, window_size, input.channels());

	int i, j;
	for (i = 0; i < window_size; i++) {
		for (j = 0; j < window_size; j++) {
			output(i, j) = input.smartAccessor(p.x + (i - floor(window_size/2)), p.y + (j - floor(window_size/2)), 0, true);
		}
	}

	return output;
}

// Returns 1's where image is filled, 0 otherwise
Image getValidMask(Image input) {
	Image output(input.width(), input.height(), 1);

	int i, j;
	for (i = 0; i < input.width(); i++) {
		for (j = 0; j < input.height(); j++) {
			if (input(i, j) > 0) {
				output(i, j) = 1.0f; 
			}
		}
	}

	return output;
}
// Finds suitable pixel match and returns it.
// If no match found, returns Pixel with value (-1, -1).
Pixel findMatch(Image temp, Image input, int window_size, Pixel p, bool is_hole) {

	Image valid_mask = getValidMask(temp);
	Image gaussian_mask = gaussian2DWindow(window_size, window_size/6.4);

	// Find total weight for normalization
	Image valid_gaussian_mask = valid_mask * gaussian_mask;
	float total_weight = valid_gaussian_mask.sum();

	Image ssd(input.width(), input.height(), 1);
	int ii, jj, i, j;
	for (i = 0; i < input.width(); i++) {
		for (j = 0; j < input.height(); j++) {
			if (!(i == p.x && j == p.y) && (input(i, j) != 0)) {
				for (ii = 0; ii < temp.width(); ii++) {
					for (jj = 0; jj < temp.height(); jj++) {
						// find squared distance between temp and input
						float curr = input.smartAccessor(i+(ii - temp.width()/2), j+(jj - temp.height()/2), 0, true);
						float dist = pow(temp(ii, jj) - curr, 2);
						ssd(i, j) += dist * valid_gaussian_mask(ii, jj);
					}
				}
				ssd(i, j) = ssd(i, j) / total_weight;
			} else {
				// if same pixel or pixel is 0, set ssd to inf so that it's not small (if 0, will affect ssd min)
				ssd(i, j) = INFINITY;
			}
		}
	}

	std::vector<Pixel> best_matches;
	for (i = 0; i < ssd.width(); i++) {
		for (j = 0; j < ssd.height(); j++) {
			if (ssd(i, j) <= ssd.min() * (1 + ERR_THRESHOLD) && ssd(i, j) < MAX_ERR_THRESHOLD) {
				Pixel pp(i, j);
				best_matches.push_back(pp);
			}
		}
	}

	// randomly select a best match pixel
	if (best_matches.size() > 0) {
		int random_index = rand() % best_matches.size();
		return best_matches[random_index];
	} else {
		Pixel no_match(-1, -1);
		return no_match;
	}
}

Image textureSynthesis(Image input, int new_width, int new_height, int window_size, bool is_hole) {
    Image output = createNewImage(new_width, new_height, input, is_hole);
    Image adjusted_input = input;
    if (!is_hole) {
    	adjusted_input = adjustInput(input);
    }
	
	int counter = 0;
   	while (true) {
   		// get pixel locations that are not filled with filled pixels as neighbors
		std::vector<Pixel> pixel_list = getUnfilledNeighbors(output);
		// if no more unfilled neighbors, return. image finished filling
		if (pixel_list.size() == 0) {
			break;
		}
   		// iterate through pixels
   		for (int p = 0; p < pixel_list.size(); p++) {
   			Pixel pixel = pixel_list[p];
   			// find neighborhood of pixel
   			Image temp = getNeighborhoodWindow(pixel, output, window_size);
   			// find matches to neighborhood in input
   			// TODO lighten input
   			Pixel best_match = findMatch(temp, adjusted_input, window_size, pixel, is_hole);
   			// randomly select a best match pixel
   			if (best_match.x != -1 && best_match.y != -1) {
   				// add slight threshold if empty pixel, which happens only when expanding
   				output(pixel.x, pixel.y) = adjusted_input(best_match.x, best_match.y);
   			}
   			if (counter % 100 == 0) {
   				output.write("./Output/intermed_" + std::to_string(counter) + ".png");
   			}
   			counter++;
   		}
   	}
   	return output;
}

Image invert(Image input) {
	Image output(input.width(), input.height(), input.channels());
	int i, j, k;
	for (i = 0; i < input.width(); i++) {
		for (j = 0; j < input.height(); j++) {
			for (k = 0; k < input.channels(); k++) {
				if (input(i, j, k) == 0) {
					output(i, j, k) = 1.0f;
				} else {
					output(i, j, k) = 0.0f;
				}
			}
		}
	}
	return output;
}

// INAPINTING FUNCTIONS

// Returns a 3-element vector with the new_image, where masked areas are empty (have value 0),
// the mask, where masked areas have value 1 and non-masked areas have value 0, and the mask inverse. 
// TODO add more color options
vector<Image> createImageAndMasks(Image input) {
	Image new_image(input.width(), input.height(), 1); 
	Image mask(input.width(), input.height(), 3);

	int i, j;
	for (i = 0; i < input.width(); i++) {
		for (j = 0; j < input.height(); j++) {
			// if masked by white infill
			if (input(i, j) == 1.0f) {
				new_image(i, j) = 0.0f;
				mask(i, j, 0) = 1.0f;
				mask(i, j, 1) = 1.0f;
				mask(i, j, 2) = 1.0f;
			} else {
				new_image(i, j) = input(i, j);
				mask(i, j, 0) = 0.0f;
				mask(i, j, 1) = 0.0f;
				mask(i, j, 2) = 0.0f;
			}
		}
	}

	// dilate the mask so that it covers a little more than region requested
	mask = dilate_filterClass(mask, 15); // TODO make this value user inputtable
	// normalize the filter values so that they are only 0 or 1
	mask.binary();

	// invert dilation
	Image mask_inv = invert(mask);

	vector<Image> output;
	output.push_back(new_image);
	output.push_back(mask);
	output.push_back(mask_inv);
	return output;
}

// Returns a "flat" white image.
Image createFlatImage(Image input) {
	Image output(input.width(), input.height(), 3);
	output.set_color(1.0f);
	return output;
}

// Returns the dot product of an image by multiplying two images element-wise, then summing the values. 
float dotProdIm(Image im1, Image im2) {
	Image product = im1 * im2;
	return product.sum();
}

// From PSET 7: panorama.cpp
Image computeTensor(Image im, Image mask, float sigmaG, float factorSigma) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Compute xx/xy/yy Tensor of an image. (stored in that order)
  // luminance 
	Image luminance = lumiChromi(im)[0];
	// blur
	Image blurredLuminance = gaussianBlur_separable(luminance, sigmaG);
	// luminance gradients
	Image gradX = gradientX(blurredLuminance);
	Image gradY = gradientY(blurredLuminance);
	Image output(im.width(), im.height(), 3);
	int i, j;
	for (i = 0; i < im.width(); i++) {
		for (j = 0; j < im.height(); j++) {
		  output(i, j, 0) = gradX(i, j, 0) * gradX(i, j, 0);
		  output(i, j, 1) = gradX(i, j, 0) * gradY(i, j, 0);
		  output(i, j, 2) = gradY(i, j, 0) * gradY(i, j, 0);
		}
	}
	// weighting
	return gaussianBlur_separableMask(output, mask, sigmaG*factorSigma);
}

// Follows a8 instructions given in handout
Image poisson(Image bg, Image fg, Image mask, Image mask_inv, int niter) {
	// apply laplacian operator to foreground
	Image b = laplacian_filterClass(fg);
	// masked areas black, non masked areas from original image
	Image x = bg * mask_inv;
	int i;
	for (i = 0; i < niter; i++) {
		// residual: apply laplacian to estimate, subtract from b
		// multiply residual by mask as soon as computed so that updates will only consider pixels inside mask
		Image r = (b - laplacian_filterClass(x)) * mask;
		// solve for alpha using dot product for images
		float alpha = dotProdIm(r, r) / dotProdIm(r, laplacian_filterClass(r));
		// update estimate of x
		x = x + alpha * r;
	}
	return x;
}

// Calculates differences between pixels in window from center pixel value. 
// Original and guide should both be structure tensors 
Image windowDiffsToCenter(Image input, Image tensor_original, Image tensor_guide, int x, int y, int window_size) {
	Image output(window_size, window_size, 3);

	int i, j, k;
	int radius = (int) floor(window_size/2);
	for (i = x - radius; i < x + radius + 1; i++) {
		for (j = y - radius; j < y + radius + 1; j++) {
			for (k = 0; k < 3; k++) {
				// gradient = value of pixel in original - value of pixel in guide
				// TODO: assumes input is one channel
				output(i - x + radius, j - y + radius, k) = input.smartAccessor(i, j, 0, true) * 1 - (abs(tensor_original.smartAccessor(i, j, k, true) - tensor_guide.smartAccessor(x, y, k, true)));
			}
		}
	}
	return output;
}

// Calculates positional diffusion coefficients using formual e^(-(gradI/k)^2)
Image diffusionCoefficients(Image window, float constant) {
	Image output(window.width(), window.height(), 3);
	int i, j, k;
	for (i = 0; i < window.width(); i++) {
		for (j = 0; j < window.height(); j++) {
			for (k = 0; k < 3; k++) {
				output(i, j) = exp(-pow(window(i, j, k)/constant, 2));
			}
		}
	}
	return output;
}

// Helpful: http://www.cs.utah.edu/~manasi/coursework/cs7960/p2/project2.html
Image anisotropicDiffusion(Image input, Image tensor_original, Image tensor_guide, float step_size) {
	Image output(input.width(), input.height(), input.channels());
	int i, j;
	for (i = 0; i < input.width(); i++) {
		for (j = 0; j < input.height(); j++) {
			// check if mask value (we only want to operate on masked areas)
			if (input(i, j) == 0.0f) {
				// calculate gradients using guide 
				Image gradient = windowDiffsToCenter(input, tensor_original, tensor_guide, i, j);
				Image coeff = diffusionCoefficients(gradient);
				// new value = old input + step size * mean of gradient * coeff elementwise multiplication
				float new_val = input(i, j) + step_size * (gradient*coeff).sum() / 9;
				printf("New val at location: %f, (%d, %d)\n", new_val, i, j);
				output(i, j) = new_val;
			} else {
				// Directly copy non-masked values
				output(i, j) = input(i, j);
			}
		}
	}
	return output;
}

Image inpaint(Image im, int niter) {
	// mask pixels in actual image, i.e. make empty
	vector<Image> imAndMasks = createImageAndMasks(im);
	Image zeroed = imAndMasks[0];
	Image mask = imAndMasks[1];
	Image mask_inv = imAndMasks[2];
	zeroed.write("./Output/zeroed.png");
	mask.write("./Output/mask.png");
	mask_inv.write("./Output/mask_inv.png");
	printf("inpaint1\n");

	// Compute flat white image
	Image flat = createFlatImage(im);
	printf("inpaint2\n");

	// compute structure tensor
	Image tensor = computeTensor(zeroed, mask);
	tensor.write("./Output/tensor.png");
	printf("inpaint3\n");

	// run Poisson with a flat source to interpolate the structure tensor inside the region.
	Image interpolated_struc_tensor = poisson(tensor, flat, mask, mask_inv);
	// Image interpolated_struc_tensor("./Input/in/poisson.png");
	interpolated_struc_tensor.write("./Output/interpolated.png");
	printf("inpaint4\n");

	// Interpolate color values using anisotropic diffusion
	Image output = anisotropicDiffusion(zeroed, tensor, interpolated_struc_tensor);
	printf("inpaint5\n");

	return output;
}

