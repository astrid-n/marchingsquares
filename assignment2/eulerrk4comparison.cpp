/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp
 *  Init    : Tuesday, September 19, 2017 - 15:08:24
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/algorithm/boundingbox.h>
#include <inviwo/core/interaction/events/mouseevent.h>
#include <labstreamlines/eulerrk4comparison.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo EulerRK4Comparison::processorInfo_{
    "org.inviwo.EulerRK4Comparison",  // Class identifier
    "Euler RK4 Comparison",           // Display name
    "KTH Lab",                        // Category
    CodeState::Experimental,          // Code state
    Tags::None,                       // Tags
};

const ProcessorInfo EulerRK4Comparison::getProcessorInfo() const { return processorInfo_; }

EulerRK4Comparison::EulerRK4Comparison()
    : Processor()
    , inData("inData")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
          MouseButton::Left, MouseState::Press | MouseState::Move)
    , maxNumberIntegrationStepsEuler("maxNumberIntegrationStepsEuler", "Max Integration Steps for Euler", 10)
    , maxNumberIntegrationStepsRK4("maxNumberIntegrationStepsRK4", "Max Integration Steps for RK4", 10)
    , propIntegrationMethod("integrationMethod", "Integration Method")
    , stepSizeEuler("stepSizeEuler", "Step Size for Euler", 1.0f)
    , stepSizeRK4("stepSizeRK4", "Step Size for RK4", 1.0f)
    , propIntegrationDirection("integrationDirection", "Integration Direction")
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(meshOut);
    addPort(meshBBoxOut);
    addPort(inData);

    // Register Properties
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);

    addProperty(propIntegrationMethod);
    propIntegrationMethod.addOption("euler", "Euler", 0);
    propIntegrationMethod.addOption("RK4", "Runge Kutta 4", 1);
    propIntegrationMethod.addOption("both", "Both", 2);

    addProperty(maxNumberIntegrationStepsEuler);
    addProperty(maxNumberIntegrationStepsRK4);
    addProperty(stepSizeEuler);
    addProperty(stepSizeRK4);

    addProperty(propIntegrationDirection);
    propIntegrationDirection.addOption("forward", "Forward", 0);
    propIntegrationDirection.addOption("backward", "Backward", 1);
    
    // TODO: Register additional properties
    // addProperty(propertyName);

}

void EulerRK4Comparison::eventMoveStart(Event* event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to bounding box range
    mousePos[0] *= static_cast<float>(BBoxMax_[0] - BBoxMin_[0]);
    mousePos[1] *= static_cast<float>(BBoxMax_[1] - BBoxMin_[1]);
    mousePos += static_cast<vec2>(BBoxMin_);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void EulerRK4Comparison::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    // The start point should be inside the volume (set maximum to the upper right corner)
    propStartPoint.setMinValue(BBoxMin_ - dvec2(1, 1));
    propStartPoint.setMaxValue(BBoxMax_ + dvec2(1, 1));

    // Initialize mesh, vertices and index buffers for the two streamlines and the points
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin_[0], BBoxMax_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax_[0], BBoxMin_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    dvec2 startPoint = propStartPoint.get();

    // TODO: Implement the Euler and Runge-Kutta of 4th order integration schemes
    // and then integrate forward for a specified number of integration steps and a given stepsize
    // (these should be additional properties of the processor)
    vec4 red = vec4(1, 0, 0, 1);
    vec4 blue = vec4(0, 0, 1, 1);

    int integrationDirection = 1; //1 if we integrate forward, -1 if we integrate backwards
    if (propIntegrationDirection.get() == 1) { //then we integrate backwards instead of forward
        integrationDirection = -1;
    }

    if(propIntegrationMethod==0 or propIntegrationMethod==2) {
        //Euler integration
        Integrator::drawNextPointInPolyline(startPoint, red, indexBufferEuler.get(), vertices);
        dvec2 currentPoint = startPoint;
        for (int i = 0; i < maxNumberIntegrationStepsEuler; i++) {
            currentPoint = Integrator::Euler(vectorField, currentPoint, integrationDirection * stepSizeEuler);
            Integrator::drawNextPointInPolyline(currentPoint, red, indexBufferEuler.get(), vertices);
            Integrator::drawPoint(currentPoint, red, indexBufferPoints.get(), vertices);
            
        }
    }
    
    if(propIntegrationMethod==1 or propIntegrationMethod==2) {
        //RK4 integration
        Integrator::drawNextPointInPolyline(startPoint, blue, indexBufferRK.get(), vertices);
        dvec2 currentPoint = startPoint;
        for (int i = 0; i < maxNumberIntegrationStepsRK4; i++) {
            currentPoint = Integrator::RK4(vectorField, currentPoint, integrationDirection * stepSizeRK4);
            Integrator::drawNextPointInPolyline(currentPoint, blue, indexBufferRK.get(), vertices);
            Integrator::drawPoint(currentPoint, blue, indexBufferPoints.get(), vertices);
        }
    }

    // Draw start point
    Integrator::drawPoint(startPoint, black, indexBufferPoints.get(), vertices);

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

}  // namespace inviwo
