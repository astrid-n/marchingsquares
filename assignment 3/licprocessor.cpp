/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
// TODO: Register additional properties
    , kernelSize("kernelSize", "Kernel Size", 1, 1, 100)
    , contrastEnhancement("contrastEnhancement", "Contrast Enhancement")
    , mean("mean", "Mean", 0.5, 0, 1)
    , standardDeviation("standardDeviation", "Standard Deviation", 0.1, 0.001, 1)
    , propLIC("LIC", "LIC")
    , propTexture("Texture", "Texture")
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    // TODO: Register additional properties
    addProperty(kernelSize);
    addProperty(contrastEnhancement);
    addProperty(mean);
    addProperty(standardDeviation);
    addProperty(propLIC);
    propLIC.addOption("normalLIC", "Normal LIC", 0);
    propLIC.addOption("fastLIC", "Fast LIC", 1);
    addProperty(propTexture);
    propTexture.addOption("gray", "Gray", 0);
    propTexture.addOption("color", "Color", 1);

    util::hide(mean, standardDeviation);
    contrastEnhancement.onChange([this]() {
        if (contrastEnhancement.get() == 0) {
            util::hide(mean, standardDeviation);
        } else {
            util::show(mean, standardDeviation);
        }
    });
}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();

    double value = texture.readPixelGrayScale(size2_t(0, 0));

    LogProcessorInfo(value);

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 0.0));

    auto mesh = std::make_shared<BasicMesh>();
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    std::vector<BasicMesh::Vertex> vertices;

    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    // TODO: Implement LIC and FastLIC*
    if (propLIC == 0) {
        LIC(licImage, vectorField, texture);
    } 
    else if(propLIC==1) { //if fast LIC
        fastLIC(licImage, vectorField, texture);       
    }


    if (contrastEnhancement.get() == 1) { //apply contrast enhancement with given mean and standard deviation value
        enhanceContrast(licImage);
    }
    

    /*
    // This code instead sets all pixels to the texture value
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            //int val = int(licTexture[i][j]);
            dvec4 vec = texture.readPixel(size2_t(i, j)); 
            licImage.setPixel(size2_t(i, j), vec); //uncomment to view the random color texture
            // or
            //int val = int(texture.readPixelGrayScale(size2_t(i, j))); 
            //licImage.setPixelGrayScale(size2_t(i, j), val); //uncomment to view the random gray texture
        }
    }
    */

    licOut_.setData(outImage);
}

dvec2 LICProcessor::vectorFieldToTextureCoords(dvec2 const position) {
    return vec2((position[0] - BBoxMin_[0]) * texDims_.x / (BBoxMax_[0] - BBoxMin_[0]), (position[1] - BBoxMin_[1]) * texDims_.y / (BBoxMax_[1] - BBoxMin_[1]));
}

dvec2 LICProcessor::textureToVectorFieldCoords(dvec2 const coords) {
    return vec2((coords[0] / texDims_.x) * (BBoxMax_[0] - BBoxMin_[0]) + BBoxMin_[0], (coords[1] / texDims_.y) * (BBoxMax_[1] - BBoxMin_[1]) + BBoxMin_[1]);
}

void LICProcessor::LIC(RGBAImage& licImage, const VectorField2& vectorField, const RGBAImage texture) {
    double stepSize = std::min((BBoxMax_[0] - BBoxMin_[0]) / (double)texDims_.x, (BBoxMax_[1] - BBoxMin_[1]) / (double)texDims_.y);
    double maxMagnitude = 0.0; 
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            //take the stream line containing the point
            dvec2 textureCoords_ij = vec2(i + 0.5 / (double)texDims_.x, j + 0.5 / (double)texDims_.y); //coordinates at the center of the pixel (i,j) in the texture
            dvec2 startPoint = textureToVectorFieldCoords(textureCoords_ij); //corresponding coordinates for the vector field
            std::vector<dvec2> streamline;
    
            streamline = Integrator::computeEquidistantStreamline(startPoint, vectorField, stepSize, kernelSize.get());
            
            double value = 0;
            for(int k = 0 ; k < streamline.size() ; k++) {
                dvec2 textureCoords_k = vectorFieldToTextureCoords(streamline[k]);
                value += texture.sampleGrayScale(textureCoords_k);
            }
            //value /= streamline.size();
            value /= (double)kernelSize.get();
            licImage.setPixelGrayScale(size2_t(i, j), value);

            if(propTexture == 1) {
                if(maxMagnitude < glm::length(vectorField.interpolate(startPoint))) {
                    maxMagnitude = glm::length(vectorField.interpolate(startPoint));
                }
            }
        }
    }

    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            if(propTexture == 1) {
                dvec2 textureCoords_ij = vec2(i + 0.5 / (double)texDims_.x, j + 0.5 / (double)texDims_.y); //coordinates at the center of the pixel (i,j) in the texture
                dvec2 position_ij = textureToVectorFieldCoords(textureCoords_ij); //corresponding coordinates for the vector field
                double magnitude = glm::length(vectorField.interpolate(position_ij)) / maxMagnitude;
                dvec4 color = licImage.readPixelGrayScale(size2_t(i, j)) * vec4(magnitude, 0, 1-magnitude, 255);
                //set pixel intensity for (i,j)
                licImage.setPixel(size2_t(i, j), color);
            }
        }
    }
}

void LICProcessor::fastLIC(RGBAImage& licImage, const VectorField2& vectorField, const RGBAImage texture) {
    // Hint: Output an image showing which pixels you have visited for debugging
    std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));
    double stepSize = std::min((BBoxMax_[0] - BBoxMin_[0]) / (double)texDims_.x, (BBoxMax_[1] - BBoxMin_[1]) / (double)texDims_.y);
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            if (visited[i][j] == 0) { //we only compute the streamline if the point has not been visited yet
                //take the stream line containing the point
                dvec2 textureCoords_ij = vec2(i + 0.5 / (double)texDims_.x, j + 0.5 / (double)texDims_.y); //coordinates at the center of the pixel (i,j) in the texture
                dvec2 startPoint = textureToVectorFieldCoords(textureCoords_ij); //corresponding coordinates for the vector field
                std::deque<dvec2> streamline = Integrator::computeEquidistantMaxStreamline(startPoint, vectorField, stepSize);
                //initialize color for the first pixel
                double value = 0;
                for (int k = 0; k < kernelSize.get()/2.0 ; k++) {
                    dvec2 textureCoords_k = vectorFieldToTextureCoords(streamline[k]);
                    value += texture.sampleGrayScale(textureCoords_k);
                }
                //for the next pixels, just shift the value by adding/removing the right texture samples
                for (int k = 0; k < streamline.size(); k++) {
                    //remove first pixel in the convolution
                    if (k > kernelSize.get()/2.0) {
                        dvec2 textureCoords_first = vectorFieldToTextureCoords(streamline[k - kernelSize.get()/2.0 - 1]);
                        value -= texture.sampleGrayScale(textureCoords_first);
                    }
                    //add last pixel in the convolution
                    if( k < kernelSize.get()/2.0) {
                        dvec2 textureCoords_last = vectorFieldToTextureCoords(streamline[k + kernelSize.get()/2.0]);
                        value += texture.sampleGrayScale(textureCoords_last);                        
                    }
                    //get the average
                    dvec2 textureCoords_k = vectorFieldToTextureCoords(streamline[k]);
                    //get the index (i_k, j_k) of the pixel corresponding to point k
                    int i_k = std::floor(textureCoords_k[0]);
                    int j_k = std::floor(textureCoords_k[1]);
                    if ((i_k >= 701) || (j_k >= 401)) {
                        std::cout << "texSize = " << texDims_ << " vs (i_k, j_k) = " << i_k << ", " << j_k << std::endl;
                        std::cout << "bbox = " << BBoxMin_ << " , " << BBoxMax_ << " vs coords = " << streamline[k] << std::endl;
                    }

                    //update texture and mark pixel as visited
                    if ((i_k < texDims_.x) && (j_k < texDims_.y)) {
                        licImage.setPixelGrayScale(size2_t(i_k, j_k), value / (double)kernelSize.get());
                        visited[i_k][j_k] = 1;                       
                    }
                }
            }
        }
    }
}

void LICProcessor::enhanceContrast(RGBAImage& licImage) {
    //compute mean and standard deviation of non-zero pixels
    double sumPixels = 0;
    double sumPixelsSquared = 0;
    int numberNonBlackPixels = 0;
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            double pixel_ij = licImage.readPixelGrayScale(size2_t(i,j)); 
            if (pixel_ij != 0) {
                sumPixels += pixel_ij;
                sumPixelsSquared += pixel_ij * pixel_ij;
                numberNonBlackPixels++;
            }
        }
    }
    double meanNonZero = sumPixels / (double)numberNonBlackPixels;
    double standardDeviationNonZero = std::sqrt((sumPixelsSquared - numberNonBlackPixels * meanNonZero * meanNonZero) / (double)(numberNonBlackPixels - 1));
    //compute desired pixel values
    double stretchingFactor = std::min(255 * standardDeviation.get() / standardDeviationNonZero, 3.0); //we restrict to the maximum value
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            double pixel_ij = licImage.readPixelGrayScale(size2_t(i,j)); 
            double newPixel_ij = 255 * mean + stretchingFactor * (pixel_ij - meanNonZero); 
            licImage.setPixelGrayScale(size2_t(i, j), newPixel_ij);
        }
    }
}

}  // namespace inviwo
