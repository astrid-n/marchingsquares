/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming
// scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propDisplayPoints("displayPoints", "Display Points", true)
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(-1), vec2(1), vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , propNumStepsTaken("numstepstaken", "Number of actual steps", 0, 0, 100000)
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
          MouseButton::Left, MouseState::Press | MouseState::Move)

// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional),
// increment (optional)); propertyIdentifier cannot have spaces
    , maxNumberIntegrationSteps("maxNumberIntegrationSteps", "Max Integration Steps", 10, 0, 100)
    , stepSize("stepSize", "Step Size", 0.1f, 0.0f)
    , propIntegrationDirection("integrationDirection", "Integration Direction")
    , propIntegrateInDirectionField("integrateInDirectionField", "Integrate in the direction field")
    , maxArcLength("maxArcLenght", "Max Arc Lenght", 10, 0, 50)
    , minVelocity("minVelocity", "Min Velocity", 0, 0, 2)
    , streamlineMode("streamlineMode", "Stream Line Mode")
    , numberOfStreamlines("numberOfStreamlines", "Number of Streamlines", 1, 1, 1000)
    , numberGridPoints("numberGridPoints", "Number of Grid Points", vec2(5, 5), vec2(0), vec2(20))

{
    // Register Ports
    addPort(inData);
    addPort(meshOut);
    addPort(meshBBoxOut);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(propDisplayPoints);
    addProperty(propNumStepsTaken);
    propNumStepsTaken.setReadOnly(true);
    propNumStepsTaken.setSemantics(PropertySemantics::Text);
    addProperty(mouseMoveStart);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(maxNumberIntegrationSteps);
    addProperty(stepSize);

    addProperty(propIntegrationDirection);
    propIntegrationDirection.addOption("forward", "Forward", 0);
    propIntegrationDirection.addOption("backward", "Backward", 1);

    addProperty(propIntegrateInDirectionField);
    addProperty(maxArcLength);
    addProperty(minVelocity);

    addProperty(streamlineMode);
    streamlineMode.addOption("random", "Random", 0);
    streamlineMode.addOption("uniformGrid", "Uniform Grid", 1);
    streamlineMode.addOption("magnitude", "Magnitude Distribution", 2);
    addProperty(numberOfStreamlines);
    addProperty(numberGridPoints);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    util::hide(streamlineMode, numberOfStreamlines, numberGridPoints);
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::hide(streamlineMode, numberOfStreamlines, numberGridPoints);
        } else {
            util::hide(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::show(streamlineMode, numberOfStreamlines, numberGridPoints);
        }
    });
    // Show/hide properties associated to stream line modes
    streamlineMode.onChange([this]() {
        if ((streamlineMode.get() == 0) || (streamlineMode.get() == 2)) {
            util::show(numberOfStreamlines);
            util::hide(numberGridPoints);
        } else if (streamlineMode.get() == 1) {
            util::hide(numberOfStreamlines);
            util::show(numberGridPoints);
        }   
    });
}

void StreamlineIntegrator::eventMoveStart(Event* event) {
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

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    //Initialize random number generator
    randGenerator.seed(static_cast<std::mt19937::result_type>(0));

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

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    if (propSeedMode.get() == 0) {
        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
        vec2 startPoint = propStartPoint.get();
        
        //draw the stream line and retrieve the number of steps integrated
        int numberStepsTaken = drawStreamLine(startPoint, vectorField, indexBufferPoints.get(), indexBufferRK.get(), vertices, propDisplayPoints.get() != 0);

        // TODO: Use the propNumStepsTaken property to show how many steps have actually been
        // integrated This could be different from the desired number of steps due to stopping
        // conditions (too slow, boundary, ...)
        propNumStepsTaken.set(numberStepsTaken);

        // Draw start point
        if (propDisplayPoints.get() != 0)
            Integrator::drawPoint(startPoint, vec4(1, 0, 0, 1), indexBufferPoints.get(), vertices);

    } else {        
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        if (streamlineMode.get() == 0) {
            for (int i = 0; i < numberOfStreamlines; i++) {
                auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                float rand_x = randomValue(BBoxMin_[0], BBoxMax_[0]);
                float rand_y = randomValue(BBoxMin_[1], BBoxMax_[1]);
                dvec2 randomStartPoint = vec2(rand_x, rand_y);
                //draw starting point
                if (propDisplayPoints.get() != 0) 
                    Integrator::drawPoint(randomStartPoint, vec4(1, 0, 0, 1), indexBufferPoints.get(), vertices);
                //draw streamline
                drawStreamLine(randomStartPoint, vectorField, indexBufferPoints.get(), indexBufferRK.get(), vertices, propDisplayPoints.get() != 0);
            }
        }
        else if ((streamlineMode.get() == 1)){
            for (int i=0; i < numberGridPoints.get()[0]; i++) {
                for (int j=0; j < numberGridPoints.get()[1]; j++) {
                    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                    dvec2 gridPoint = vec2(BBoxMin_[0] + (i+1)*(BBoxMax_[0]-BBoxMin_[0])/(numberGridPoints.get()[0]+1), BBoxMin_[1] + (j+1)*(BBoxMax_[1]-BBoxMin_[1])/(numberGridPoints.get()[1]+1));
                    //draw starting point
                    if (propDisplayPoints.get() != 0) 
                        Integrator::drawPoint(gridPoint, vec4(1, 0, 0, 1), indexBufferPoints.get(), vertices);
                    //draw streamline
                    drawStreamLine(gridPoint, vectorField, indexBufferPoints.get(), indexBufferRK.get(), vertices, propDisplayPoints.get() != 0);
                    }
            }
        }
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
        else if ((streamlineMode.get() == 2)) {
            //associate to each cell of a uniform grid a probability of being chosen that is proportional to the magnitude at the center of the cell
            int gridResolution = 50; //number of cells in each direction
            double cellSize = (BBoxMax_[0]-BBoxMin_[0])/gridResolution;
            std::vector<std::vector<double> > samplingProbability = computeSamplingProbability(gridResolution, vectorField);
            //use these probabilities to sample the points
            for (int i = 0; i < numberOfStreamlines; i++) {
                auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                ivec2 randomCellIndex = randomCell(samplingProbability);
                float rand_x = randomValue(BBoxMin_[0] + randomCellIndex[0]*cellSize, BBoxMin_[0] + (randomCellIndex[0]+1)*cellSize);
                float rand_y = randomValue(BBoxMin_[1] + randomCellIndex[1]*cellSize, BBoxMin_[1] + (randomCellIndex[1]+1)*cellSize);
                dvec2 randomStartPoint = vec2(rand_x, rand_y);
                //draw starting point
                if (propDisplayPoints.get() != 0) 
                    Integrator::drawPoint(randomStartPoint, vec4(1, 0, 0, 1), indexBufferPoints.get(), vertices);
                //draw streamline
                drawStreamLine(randomStartPoint, vectorField, indexBufferPoints.get(), indexBufferRK.get(), vertices, propDisplayPoints.get() != 0);
            }
        }
    }

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

int StreamlineIntegrator::drawStreamLine(const dvec2 startPoint, const VectorField2 vectorField, IndexBufferRAM* indexBufferPoints, IndexBufferRAM* indexBufferLines, std::vector<BasicMesh::Vertex>& vertices, bool displayPoints) {
    /* Parameters:
        - startPoint: starting point of the stream line
        - vectorField: vector field
        - indexBufferPoints: index buffer to draw the points of the stream line
        - indexBufferLines: index buffer to draw the polyline between the vertices of the stream line
        - vertices: vector of vertices to draw
        - displayPoints: if true, we display the points on the stream line

    Output: number of steps actually integrated
    */
    // TODO: Create one stream line from the given start point
    vec4 black = vec4(0, 0, 0, 1);
    //get integration direction
    int integrationDirection = 1; //1 if we integrate forward, -1 if we integrate backwards
    if (propIntegrationDirection.get() == 1) { //then we integrate backwards instead of forward
        integrationDirection = -1;
    }  
    //we use RK4 integration
    Integrator::drawNextPointInPolyline(startPoint, black, indexBufferLines, vertices);
    dvec2 currentPoint = startPoint;
    double arcLength = 0;
    int numberStepsTaken = 0;
    for (int i = 0; i < maxNumberIntegrationSteps; i++) {
        dvec2 nextPoint = Integrator::RK4(vectorField, currentPoint, integrationDirection * stepSize, propIntegrateInDirectionField.get());
        double velocity = glm::length(currentPoint - nextPoint);
        arcLength += velocity;
        //std::cout << "i=" << i << std::endl;
        //std::cout << "arcLength = " << arcLength << " VS maxArcLength = " << maxArcLength << std::endl;
        //std::cout << "velocity = " << velocity << " VS minVelocity = " << minVelocity << std::endl;
        if ((nextPoint[0] <= BBoxMax_[0]) && (nextPoint[0] >= BBoxMin_[0])&&((nextPoint[1] <= BBoxMax_[1]))&&(nextPoint[1] >= BBoxMin_[1])&&
        (arcLength < maxArcLength.get())&&(velocity != 0)&&(velocity >= minVelocity)) {
            currentPoint = nextPoint;
            numberStepsTaken += 1;
            Integrator::drawNextPointInPolyline(nextPoint, black, indexBufferLines, vertices);
            if (displayPoints){
                Integrator::drawPoint(nextPoint, black, indexBufferPoints, vertices);
            }
        }
        else {
            break;
        }
    }

    return numberStepsTaken;
}

float StreamlineIntegrator::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}

std::vector<std::vector<double> > StreamlineIntegrator::computeSamplingProbability(int gridResolution, const VectorField2 vectorField) {
    //create a uniform grid inside the bounding box and sample the vector field magnitude in the center of each cell
    std::vector<std::vector<double> > samplingProbability(gridResolution, std::vector<double>(gridResolution));
    double sumMagnitudes = 0; //compute the sum of all magnitudes and divide every magnitude at the end to get a probability
    for (int i = 0; i < gridResolution; i++) {
        for (int j=0; j < gridResolution; j++) { //sample vector field in cell (i,j)
            dvec2 centerCell = vec2(BBoxMin_[0] + (i+0.5)*(BBoxMax_[0]-BBoxMin_[0])/gridResolution, BBoxMin_[1] + (j+0.5)*(BBoxMax_[1]-BBoxMin_[1])/gridResolution);
            samplingProbability[i][j] = glm::length(vectorField.interpolate(centerCell));
            sumMagnitudes += samplingProbability[i][j];
        }
    }
    //divide by the sum of magnitudes to get the probability
    for (int i = 0; i < gridResolution; i++) {
        for (int j=0; j < gridResolution; j++) { //sample vector field in cell (i,j)
        samplingProbability[i][j] /= sumMagnitudes;
        }
    }
    return samplingProbability;
}

ivec2 StreamlineIntegrator::randomCell(std::vector<std::vector<double> > samplingProbability) {
    float rand = randomValue(0, 1);
    double cumulatedProbability = 0;
    for (int i = 0; i < samplingProbability.size(); i++) {
        for (int j = 0; j < samplingProbability[0].size(); j++) {
            cumulatedProbability += samplingProbability[i][j];
            if (rand <= cumulatedProbability) 
                return vec2(i, j);
        }
    }
}

}  // namespace inviwo
