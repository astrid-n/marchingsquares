/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>

#include <iostream>
#include <queue> 

namespace inviwo {

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1)     // Center - Green
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",    // Class identifier
    "Vector Field Topology",  // Display name
    "KTH Lab",                // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const { return processorInfo_; }

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , meshBBoxOut("meshBBoxOut")
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
}

void Topology::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Initialize mesh, vertices and index buffers for seperatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    /*
    //change of sign test
    dvec2 v1 = vectorField.interpolate(vec2(BBoxMin[0], BBoxMin[1])); //down left
    dvec2 v2 = vectorField.interpolate(vec2(BBoxMax[0], BBoxMin[1])); //down right
    dvec2 v3 = vectorField.interpolate(vec2(BBoxMax[0], BBoxMax[1])); //up right
    dvec2 v4 = vectorField.interpolate(vec2(BBoxMin[0], BBoxMax[1])); //up left

    double err = pow(10,5)*/


    size2_t dims = vectorField.getNumVerticesPerDim();
    //std::cout << "vector fild dimensions ";
    //std::cout << dims << endl;

    // Looping through all values in the vector field.
    for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t i = 0; i < dims[0]; ++i) {
            dvec2 vectorValue = vectorField.getValueAtVertex(size2_t(i, j));
        }
    }

    // Other helpful functions
    // dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
    // Computing the jacobian at a position
    // dmat2 jacobian = vectorField.derive(pos);
    // Doing the eigen analysis
    // auto eigenResult = util::eigenAnalysis(jacobian);
    // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
    // eigenvectors

    // Accessing the colors
    vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                               IndexBufferRAM* indexBuffer,
                               std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

std::vector<dvec2> Topology::extractCriticalPoints(const VectorField2 vectorField, const double eps) {
    std::vector<dvec2> criticalPoints;
    size2_t dims = vectorField.getNumVerticesPerDim();

    // Looping through all values in the vector field
    for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t i = 0; i < dims[0]; ++i) {
            //for each cell (i,j), we do domain decomposition to check for a critical point
            //we initialize a queue containing the corners of the cell (bottomLeft, topRight)
            //at every step, we dequeue a cell and if it may contain a zero, we add its subcells to the queue
            //we stop when one of the corners is close enough to zero
            std::vector<dvec2> initial_corners = {cellPosition(i,j, vectorField), cellPosition(i+1, j+1, vectorField)};
            std::queue<std::vector<dvec2> > cornersQueue;
            cornersQueue.push(initial_corners);
            while (!cornersQueue.empty()) {
                //extract one subcell from the queue
                std::vector<dvec2> bottomLeftAndTopRight = cornersQueue.front();
                dvec2 bottomLeft = bottomLeftAndTopRight[0];
                dvec2 topRight = bottomLeftAndTopRight[1];
                cornersQueue.pop();
                //get all four corners
                std::vector<dvec2> corners;
                corners.push_back(bottomLeft); //bottom left
                corners.push_back(vec2(topRight[0], bottomLeft[1])); //bottom right
                corners.push_back(vec2(bottomLeft[0], topRight[1])); //top left
                corners.push_back(topRight); //top right
                //check for zero at the corner or for possible zero in one of the subcells
                bool zeroFound = false;
                ivec2 signChanges = vec2(0, 0); //signChanges[i] = 1 iif component i has a sign change
                findZerosAndSignChanges(criticalPoints, zeroFound, signChanges, corners, vectorField, eps);
                //if there is a possible zero, queue all four subcells
                if (!zeroFound) {//if a zero was found, we do not queue anything
                    if ((signChanges[0] == 1) && (signChanges[1] == 1)) {//sign change => possible zeros
                        //add subcells to queue (easier to visualize with a drawing)
                        dvec2 center = 0.5*(bottomLeft+topRight);
                        cornersQueue.push({bottomLeft, center}); //bottom left cell
                        cornersQueue.push({vec2(center[0], bottomLeft[1]), vec2(topRight[0], center[1])}); //bottom right cell
                        cornersQueue.push({vec2(bottomLeft[0], center[1]), vec2(center[0], topRight[1])}); //top left cell
                        cornersQueue.push({center, topRight});
                    }

                }
            }

        }
    }
    return criticalPoints;

}

void Topology::findZerosAndSignChanges(std::vector<dvec2>& criticalPoints, bool& zeroFound, ivec2& signChanges, const std::vector<dvec2> corners, const VectorField2 vectorField, const double eps) {
    ivec2 lastSign = vec2(0, 0); //sign(0) = 0, sign(x < 0) = -1, sign(x > 0) = 1
    for (int i_corner = 0; i_corner < 4; i_corner++) {
        vec2 valueAtCorner = vectorField.interpolate(corners[i_corner]);
        //check if value is close enough to zero
        if (glm::length(valueAtCorner) <= eps) {//zero found
            criticalPoints.push_back(corners[i_corner]);
            zeroFound = true;
            break;
        }
        else { //check for sign change and update related variables
            //get signs of both components of the vector
            ivec2 currentSign = vec2(sgn(valueAtCorner[0]), sgn(valueAtCorner[1]));
            for (int i_component = 0; i_component < 2; i_component++) {
                //check for sign change
                if (currentSign[i_component]*lastSign[i_component] < 0) {//sign change
                    signChanges[i_component] = true;
                }
                //update lastSign[i_component] only if currentSign[i_component] is non-zero
                if (currentSign[i_component] != 0) {
                    lastSign[i_component] = currentSign[i_component];
                }
            }
        }                   
    }
}


vec4 Topology::classifyCriticalPoints(const vec2 pos, const VectorField2 vectorField) {
    dmat2 jacobian = vectorField.derive(pos);
    auto eigenResult = util::eigenAnalysis(jacobian);
    double R1 = eigenResult.eigenvaluesRe[0];
    double R2 = eigenResult.eigenvaluesRe[1];
    double I1 = eigenResult.eigenvaluesIm[0];
    double I2 = eigenResult.eigenvaluesIm[1];
   
    
    // Calculate eigenvalues, real and imaginary parts R1, R2, I1, I2
    if ((I1==I2)&&(I1==0)&&(R1<0)&&(R2>0)){
        return ColorsCP[0];}//'Saddlepoint'
    if ((I1==I2)&&(I1==0)&&(R1>0)&&(R2>0)){
        return ColorsCP[1];} //'Repelling Node'
    if ((I1==I2)&&(I1==0)&&(R1<0)&&(R2<0)){
        return ColorsCP[2];} //'Attracting Node'
    if ((R1==R2)&&(R1==0)&&(I1==-I2)&&(I2!=0)){
        return ColorsCP[5];} //'Center'
    if ((R1==R2)&&(R1>0)&&(I1==-I2)&&(I2!=0)){
        return ColorsCP[3];} //'Attracting Focus'
    if ((R1==R2)&&(R1<0)&&(I1==-I2)&&(I2!=0)){
        return ColorsCP[4];} //'Repelling Focus'
    // Each class should somehow correspond to a color in the ColorsCP array with typeCP enum value
    


}

vec2 Topology::cellPosition(int i, int j, const VectorField2 vectorField) {
    const dvec2 BBoxMin = vectorField.getBBoxMin();
    const dvec2 BBoxMax = vectorField.getBBoxMax();
    return vec2(BBoxMin[0] + (BBoxMax[0] - BBoxMin[0])*i, BBoxMin[1] + (BBoxMax[1] - BBoxMin[1])*j);
}

int Topology::sgn(double val) { //returns the sign (-1, 0 or 1) of a double
    return (0 < val) - (val < 0);
}

}  // namespace inviwo
