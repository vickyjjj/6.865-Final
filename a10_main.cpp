#include <iostream>

#include "a10.h"

using namespace std;

void testDilateFilter() {
	Image cloth("./Input/ts/cloth.png");
	Image larger_cloth = createNewImage(cloth.width() * 2, cloth.height() * 2, cloth, true);
	Image output = dilate_filterClass(larger_cloth, 3, true);
	output.write("./Output/cloth_dilateFilter.png");
}

void testCreateNewImage() {
	Image cloth("./Input/ts/cloth.png");
	Image output = createNewImage(cloth.width() * 2, cloth.height() * 2, cloth, true);
	output.write("./Output/cloth_createNewImage.png");
}

void testGetUnfilledNeighbors() {
	Image bread("./Input/ts/bread.png");
	Image larger_cloth = createNewImage(bread.width(), bread.height(), bread, false);
	std::vector<Pixel> unfilledNeighbors = getUnfilledNeighbors(larger_cloth);

	// Highlight unfilled neighbor pixels in red 
	int i;
	for (i = 0; i < unfilledNeighbors.size(); i++) {
		Pixel p = unfilledNeighbors[i];
		printf("Pixel: (%d, %d)\n", p.x, p.y);
		larger_cloth(p.x, p.y, 0) = 1;
	}
	larger_cloth.write("./Output/cloth_getUnfilledNeighbors.png");
}

void testGetNeighborhoodWindow() {
	Image input("./Input/ts/bread.png");
	Image larger_input = createNewImage(input.width(), input.height(), input, false);
	Pixel p(55, 31); // test with top left corner
	Image window = getNeighborhoodWindow(p, larger_input, 9);
	window.write("./Output/bread_getNeighborhoodWindow.png");
}

void testGetValidMask() {
	Image cloth("./Input/ts/bread.png");
	Image larger_cloth = createNewImage(cloth.width(), cloth.height(), cloth, true);
	Image output = getValidMask(larger_cloth);
	output.write("./Output/cloth_getValidMask.png");
}

void testGaussian2DWindow() {
	Image output = gaussian2DWindow(15, 15.0f/6.4f);
	output.write("./Output/gaussianWindow_5.png");
	// print values
	for (int i = 0; i < 15; i++) {
		for (int j = 0; j < 15; j++) {
			printf("%f ", output(i, j));
		}
		printf("\n");
	}
}

void testFindMatch() {
	Image cloth("./Input/ts/bread.png");
	Image larger_cloth = createNewImage(cloth.width(), cloth.height(), cloth, true);
	Pixel p(55, 31); // test with top left corner
	Image temp = getNeighborhoodWindow(p, larger_cloth, 13);
	Pixel match = findMatch(temp, larger_cloth, 13, p, true);

	Image output(larger_cloth.width(), larger_cloth.height(), 3);
	// Highlight match in red
	output(match.x, match.y, 0) = 1;
	output(match.x, match.y, 1) = 0;
	output(match.x, match.y, 2) = 0;

	output.write("./Output/cloth_findMatch.png");
}

void testtextureSynthesisBW() {
	// Image grid("./Input/ts/grid.png");
	// Image output_grid = textureSynthesis(grid, grid.width() * 2, grid.height() * 2, 13, false);
	// output_grid.write("./Output/grid_grow.png");

	Image cloth("./Input/ts/cloth.png");
	Image output_cloth = textureSynthesis(cloth, cloth.width() * 2, cloth.height() * 2, 13, false);
	output_cloth.write("./Output/cloth_grow.png");

	Image text("./Input/ts/text.png");
	Image output_text = textureSynthesis(text, text.width() * 2, text.height() * 2, 25, false);
	output_text.write("./Output/text_grow.png");
}

void testtextureSynthesisColor() {
	// Image brick("./Input/ts/brick-color.png");
	// Image output_brick = textureSynthesis(brick, brick.width() * 2, brick.height() * 2, 13, false);
	// output_brick.write("./Output/brick_grow.png");

	Image wood("./Input/ts/wood-color.png");
	Image output_wood = textureSynthesis(wood, wood.width() * 2, wood.height() * 2, 13, false);
	output_wood.write("./Output/wood_grow.png");
}

void testHoleFill() {
	Image bread("./Input/ts/bread.png");
	Image output_bread = textureSynthesis(bread, bread.width(), bread.height(), 13, true);
	output_bread.write("./Output/bread_grow.png");
}

// INPAINTING TESTS

void testLaplacian() {
	Image cambridge("./Input/in/cambridge.png");
	Image output = laplacian_filterClass(cambridge);
	output.write("./Output/laplacian.png");
}

void testInvert() {
	// TODO
}

void testBinary() {
	// TODO 
}

void testComputeTensor() {
	Image cambridge("./Input/in/cambridge.png");
	vector<Image> intermed = createImageAndMasks(cambridge);
	Image tensor = computeTensor(intermed[0], intermed[1]);
	tensor.write("./Output/masked.png");
}

void testCreateImageAndMasks() {
	Image cambridge("./Input/in/cambridge.png");
	vector<Image> output = createImageAndMasks(cambridge);
	output[0].write("./Output/masked_input.png");
	output[1].write("./Output/mask.png");
	output[2].write("./Output/mask_inverse.png");
}

void testCreateFlatImage() {
	Image cambridge("./Input/in/cambridge.png");
	Image flatIm = createFlatImage(cambridge);
	flatIm.write("./Output/flatImage.png");
}

void testDotProdIm() {
	Image empty(10, 10, 1);
	Image flatIm = createFlatImage(empty);
	float output = dotProdIm(flatIm, flatIm);
	cout << "Expected output: 100" << endl;
	cout << "Observed output: " << output << endl;
}

void testPoisson() {
	Image cambridge("./Input/in/cambridge.png");
	// mask pixels in actual image, i.e. make empty
	vector<Image> imAndMasks = createImageAndMasks(cambridge);
	Image zeroed = imAndMasks[0];
	Image mask = imAndMasks[1];
	Image mask_inv = imAndMasks[2];

	// compute structure tensor
	Image tensor = computeTensor(zeroed, mask);
	tensor.write("./Output/masked.png");

	// Compute flat white image
	// Image flat(cambridge.width(), cambridge.height(), cambridge.channels());
	Image flat = createFlatImage(cambridge);

	Image output = poisson(tensor, flat, mask, mask_inv, 10000);
	output.write("./Output/poisson.png");
}

void testInpaintBW() {
	Image cambridge("./Input/in/cambridge.png");
	Image output = inpaint(cambridge);
	output.write("./Output/inpaint.png");
}

void testWindowDiffsToCenter() {
	Image cambridge("./Input/in/cambridge.png");
	// TODO 
}

void testDiffusionCoefficients() {
	Image cambridge("./Input/in/cambridge.png");
	// TODO 
}

void testGaussianBlur_separableMask() {
	Image cambridge("./Input/in/cambridge.png");
	vector<Image> imAndMasks = createImageAndMasks(cambridge);
	Image zeroed = imAndMasks[0];
	Image mask = imAndMasks[1];
	Image mask_inv = imAndMasks[2];
	Image output = gaussianBlur_separableMask(cambridge, mask);
	output.write("./Output/cambridge-gaussianWithMask.png");
}

void testFilterConvolveWithoutMasked() {
	Image cambridge("./Input/in/cambridge.png");
	// TODO 
}

int main()
{
	// TEXTURE SYNTHESIS TESTS
    // testCreateNewImage();
	// testDilateFilter();
    // testGetUnfilledNeighbors();
    // testGetNeighborhoodWindow();
    // testGetValidMask();
    // testGaussian2DWindow();
    // testFindMatch();
    // testtextureSynthesisBW();
    // testtextureSynthesisColor();
    testHoleFill();

    // INPAINTING TESTS
    // testLaplacian();
    // testInvert(); 
    // testBinary(); 
    // testComputeTensor();
    // testCreateImageAndMasks();
    // testCreateFlatImage();
    // testDotProdIm();
    // testWindowDiffsToCenter();
    // testDiffusionCoefficients();
    // testPoisson();
    // testInpaintBW();
    // testGaussianBlur_separableMask();
    // testFilterConvolveWithoutMasked();
    
    return EXIT_SUCCESS;
}
