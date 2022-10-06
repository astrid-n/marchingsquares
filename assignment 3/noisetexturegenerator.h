/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:36
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/optionproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <lablic/lablicmoduledefine.h>
#include <random>

namespace inviwo {

/** \docpage{org.inviwo.NoiseTextureGenerator, Noise Texture Generator}
    ![](org.inviwo.NoiseTextureGenerator.png?classIdentifier=org.inviwo.NoiseTextureGenerator)

    Generates a noise texture with a given number of pixels in x- and y-direction.

    ### Outports
      * __outImage__ Generated texture.

    ### Properties
      * __texSize__ Size of the texture to be generated.
*/
class IVW_MODULE_LABLIC_API NoiseTextureGenerator : public Processor {
    // Friends
    // Types
public:
    // Construction / Deconstruction
public:
    NoiseTextureGenerator();
    virtual ~NoiseTextureGenerator() = default;

    // Methods
public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

protected:
    /// Our main computation function
    virtual void process() override;

    float randomValue(const float min, const float max) const;

    // Ports
public:
    ImageOutport texOut_;

    // Properties
public:
    IntVec2Property texSize_;
    IntProperty randomSeed_;
    // TODO: Declare additional properties
    TemplateOptionProperty<int> propTexture;

    // IntProperty properyName;
    // TemplateOptionProperty<int> propertyName2;

    // Attributes
private:
    mutable std::mt19937 randGenerator;
    mutable std::uniform_real_distribution<float> uniformReal;
};

}  // namespace inviwo
