/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:36
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/noisetexturegenerator.h>
#include <labutils/rgbaimage.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo NoiseTextureGenerator::processorInfo_{
    "org.inviwo.NoiseTextureGenerator",  // Class identifier
    "Noise Texture Generator",           // Display name
    "KTH Labs",                          // Category
    CodeState::Experimental,             // Code state
    Tags::None,                          // Tags
};

const ProcessorInfo NoiseTextureGenerator::getProcessorInfo() const { return processorInfo_; }

NoiseTextureGenerator::NoiseTextureGenerator()
    : Processor()
    , texOut_("texOut")
    , texSize_("texSize", "Texture Size", vec2(512, 512), vec2(1, 1), vec2(2048, 2048), vec2(1, 1))
// TODO: Register additional properties
    , propTexture("Texture", "Texture")
    , randomSeed_("randomSeed", "Random Seed", 1, 1, 100000)
{
    // Register ports
    addPort(texOut_);

    // Register properties
    addProperty(texSize_);
    addProperty(randomSeed_);

    // TODO: Register additional properties
    addProperty(propTexture);
    propTexture.addOption("gray", "Gray", 0);
    propTexture.addOption("blackAndWhite", "Black & White", 1);
    //propTexture.addOption("color", "Color", 2); //not meaningful
}

void NoiseTextureGenerator::process() {
    // The output of the generation process is an Image
    // With the given dimensions
    // With the data format DataVec4UInt8, this means values for RGB-alpha range between 0 and 255
    auto outImage =
        std::make_shared<Image>(size2_t(texSize_.get().x, texSize_.get().y), DataVec4UInt8::get());

    // Similar to ScalarField and VectorField, the RGBAImage has some methods to sample from and set
    // values
    RGBAImage noiseTexture(outImage);

    // Setting pixels in the image/texture
    // setPixelGrayScale will set the value to (val,val,val,255) at the pixel with indices (i,j)
    int val = 4;
    noiseTexture.setPixelGrayScale(size2_t(0, 0), val);
    // setPixel allows to set all color components (red,green,blue,alpha) at the pixel with indices
    // (i,j)
    noiseTexture.setPixel(size2_t(0, 0), vec4(val, val, val, 255));

    // Reading from the image
    // readPixelGrayScale returns the averge of the three colors (red+green+blue)/3 at the pixel
    // with indices (i,j)
    double value = noiseTexture.readPixelGrayScale(size2_t(0, 0));
    // readPixel returns all color components (red,green,blue,alpha) at the pixel with indices (i,j)
    dvec4 color = noiseTexture.readPixel(size2_t(0, 0));
    LogProcessorInfo("The color at index (0,0) is " << color << " with grayscale value " << value
                                                    << ".");
    // sample peforms bilinear interpolation. For (0.5,0.5) this would involve the values at pixels
    // (0,0), (1,0), (0,1), and (1,1)
    color = noiseTexture.sample(dvec2(0.5, 0.5));
    // The grayscale version again does the same but returns an average of the three color values
    value = noiseTexture.sampleGrayScale(dvec2(0.5, 0.5));
    LogProcessorInfo("The interpolated color at (0.5,0.5) is " << color << " with grayscale value "
                                                               << value << ".");


    //Initialize random number generator

    randGenerator.seed(static_cast<std::mt19937::result_type>(randomSeed_.get()));
    for (int j = 0; j < texSize_.get().y; j++) {
        for (int i = 0; i < texSize_.get().x; i++) {
            // TODO: Randomly sample values for the texture, this produces the same gray value for
            // all pixels
            // A value within the ouput image is set by specifying pixel position and color
            //noiseTexture.setPixelGrayScale(size2_t(i, j), val);
            // Alternatively, the entire color can be specified
            int randomRed = (int) randomValue(0, 255);
            int randomGreen = (int) randomValue(0, 255);
            int randomBlue = (int) randomValue(0, 255);
            if(propTexture == 0) { //if texture is Gray
                noiseTexture.setPixelGrayScale(size2_t(i, j), randomRed); //randomRed used as the randdom gray value
            }
            else if (propTexture == 1) { //if texture is Black and White
                if(randomRed > 127) {
                    randomRed = 255;
                }
                else {
                    randomRed = 0;
                }
                noiseTexture.setPixelGrayScale(size2_t(i, j), randomRed);
            }
            //else { //if texture is Color
            //    noiseTexture.setPixel(size2_t(i, j), vec4(randomRed, randomGreen, randomBlue, 255));
            //}
        }
    }
    texOut_.setData(outImage);
}

float NoiseTextureGenerator::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}

}  // namespace inviwo
