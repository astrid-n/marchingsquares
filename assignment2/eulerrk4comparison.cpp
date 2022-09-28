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
    , maxNumberIntegrationSteps("maxNumberIntegrationSteps", "Max Integration Steps", 0)
    , propIntegrationMethod("integrationMethod", "Integration Method")
    , stepSize("stepSize", "Step Size", 1.0f)
    , propIntegrationDirection("integrationDirection", "Integration Direction")
    , propIntegrateInDirectionField("integrateInDirectionField", "Integrate in the direction field")
    , maxArcLenght("maxArcLenght", "Max Arc Lenght", 0)
    , stopAtBoundary("stopAtBoundary", "Stop At Boundary")
    , minVelocity("minVelocity", "Min Velocity")
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

    addProperty(maxNumberIntegrationSteps);
    addProperty(stepSize);

    addProperty(propIntegrationDirection);
    propIntegrationDirection.addOption("forward", "Forward", 0);
    propIntegrationDirection.addOption("backward", "Backward", 1);

    addProperty(propIntegrateInDirectionField);
    addProperty(maxArcLenght);
    addProperty(stopAtBoundary);
    addProperty(minVelocity);
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

    // Draw start point
    dvec2 startPoint = propStartPoint.get();
    dvec2 nextPoint;
    Integrator::drawPoint(startPoint, black, indexBufferPoints.get(), vertices);

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
        double arcLenght = 0;
        Integrator::drawNextPointInPolyline(startPoint, black, indexBufferEuler.get(), vertices);
        dvec2 currentPoint = startPoint;
        for (int i = 0; i < maxNumberIntegrationSteps; i++) {
            nextPoint = Integrator::Euler(vectorField, currentPoint, integrationDirection * stepSize);
            arcLenght += sqrt(pow((currentPoint[0]-nextPoint[0]),2) + pow((currentPoint[1]-nextPoint[1]),2));
            if(arcLenght < maxArcLenght.get()) {
                currentPoint = nextPoint;
                Integrator::drawNextPointInPolyline(currentPoint, red, indexBufferEuler.get(), vertices);
                Integrator::drawPoint(currentPoint, red, indexBufferPoints.get(), vertices);
            }
            
        }
    }
    
    if(propIntegrationMethod==1 or propIntegrationMethod==2) {
        //RK4 integration
        double arcLenght = 0;
        Integrator::drawNextPointInPolyline(startPoint, black, indexBufferRK.get(), vertices);
        dvec2 currentPoint = startPoint;
        for (int i = 0; i < maxNumberIntegrationSteps; i++) {
            nextPoint = Integrator::RK4(vectorField, currentPoint, integrationDirection * stepSize);
            arcLenght += sqrt(pow((currentPoint[0]-nextPoint[0]),2) + pow((currentPoint[1]-nextPoint[1]),2));
            if(arcLenght < maxArcLenght.get()) {
                currentPoint = nextPoint;
                Integrator::drawNextPointInPolyline(currentPoint, blue, indexBufferRK.get(), vertices);
                Integrator::drawPoint(currentPoint, blue, indexBufferPoints.get(), vertices);
            }
        }
    }

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

}  // namespace inviwo
