/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, Anke Friederici
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_{
    "org.inviwo.MarchingSquares",  // Class identifier
    "Marching Squares",            // Display name
    "KTH Lab",                     // Category
    CodeState::Experimental,       // Code state
    Tags::None,                    // Tags
};

const ProcessorInfo MarchingSquares::getProcessorInfo() const { return processorInfo_; }

MarchingSquares::MarchingSquares()
    : Processor()
    , inData("volumeIn")
    , meshIsoOut("meshIsoOut")
    , meshGridOut("meshGridOut")
    , propShowGrid("showGrid", "Show Grid")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f),
                    vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                    PropertySemantics::Color)
    , propDeciderType("deciderType", "Decider Type")
    , propRandomSeed("seed", "Random Seed", 0, 0, std::mt19937::max())
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f),
                   vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData)
    , propGaussianFilter("gaussianFilter", "Apply Gaussian Filter") {
    // Register ports
    addPort(inData);
    addPort(meshIsoOut);
    addPort(meshGridOut);

    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);

    addProperty(propDeciderType);
    propDeciderType.addOption("asymptotic", "Asymptotic", 0);
    propDeciderType.addOption("random", "Random", 1);

    addProperty(propRandomSeed);
    propRandomSeed.setSemantics(PropertySemantics::Text);

    addProperty(propMultiple);

    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    //we add a checkbox for Gaussian filtering
    addProperty(propGaussianFilter);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.get().add(1.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propRandomSeed, propNumContours, propIsoTransferFunc);

    propDeciderType.onChange([this]() {
        if (propDeciderType.get() == 1) {
            util::show(propRandomSeed);
        } else {
            util::hide(propRandomSeed);
        }
    });

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]() {
        if (propShowGrid.get()) {
            util::show(propGridColor);
        } else {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]() {
        if (propMultiple.get() == 0) {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        } else {
            //util::hide(propIsoValue);
            //util::show(propIsoColor, propNumContours);

            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            util::hide(propIsoValue, propIsoColor);
            util::show(propNumContours, propIsoTransferFunc);
        }
    });
}

void MarchingSquares::process() {
    if (!inData.hasData()) {
        return;
    }

    // Create a structured grid from the input volume
    auto vol = inData.getData();
    auto grid = ScalarField2::createFieldFromVolume(vol);

    // Extract the minimum and maximum value from the input data
    const double minValue = grid.getMinValue();
    const double maxValue = grid.getMaxValue();

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue
                                                                  << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from

    // Get the definition of our structured grid with
    // - number of vertices in each dimension {nx, ny}
    const ivec2 nVertPerDim = grid.getNumVerticesPerDim();
    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = grid.getBBoxMin();
    const dvec2 bBoxMax = grid.getBBoxMax();
    // - cell size {dx, dy}
    const dvec2 cellSize = grid.getCellSize();

    // Values at the vertex positions can be accessed by the indices of the vertex
    // with index i ranging between [0, nx-1] and j in [0, ny-1]
    ivec2 ij = {0, 0};
    double valueAt00 = grid.getValueAtVertex(ij);
    LogProcessorInfo("The value at (0,0) is: " << valueAt00 << ".");

    // Initialize the output: mesh and vertices for the grid and bounding box
    auto gridmesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> gridvertices;

    auto indexBufferBBox = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // bottomLeft to topLeft
    drawLineSegment(bBoxMin, vec2(bBoxMin[0], bBoxMax[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topLeft to topRight
    drawLineSegment(vec2(bBoxMin[0], bBoxMax[1]), bBoxMax, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topRight to bottomRight
    drawLineSegment(bBoxMax, vec2(bBoxMax[0], bBoxMin[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // bottomRight to bottomLeft
    drawLineSegment(vec2(bBoxMax[0], bBoxMin[1]), bBoxMin, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);

    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));
    // You can create a random sample between min and max with
    float minRand = 0.0;
    float maxRand = 1.0;
    float rand = randomValue(minRand, maxRand);
    LogProcessorInfo("The first random sample for seed " << propRandomSeed.get() << " between "
                                                         << minRand << " and " << maxRand << " is "
                                                         << rand << ".");

    // Properties are accessed with propertyName.get()
    if (propShowGrid.get()) {
        // TODO: Add grid lines of the given color

        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up
        // lines/trianges/quads. Here two vertices make up a line segment.
        auto indexBufferGrid = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        //vertical lines
        std::cout << "nVertexPerDim " << nVertPerDim << std::endl;
        std::cout << "bBoxMin " << bBoxMin << std::endl;
        std::cout << "bBoxMax " << bBoxMax << std::endl;
        for (int i = 0; i < nVertPerDim[0]; i++) {
            vec2 v1 = vec2(bBoxMin[0] + i*cellSize[0], bBoxMin[1]);
            vec2 v2 = vec2(bBoxMin[0] + i*cellSize[0], bBoxMax[0]);
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);
        }

        //horizontal lines
        for (int i = 0; i < nVertPerDim[1]; i++) {
            vec2 v1 = vec2(bBoxMin[0], bBoxMin[1] + i*cellSize[1]);
            vec2 v2 = vec2(bBoxMax[0], bBoxMin[1] + i*cellSize[1]);
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);
        }

        // Draw a line segment from v1 to v2 with a the given color for the grid
        // vec2 v1 = vec2(10, 0.5);
        // vec2 v2 = vec2(0.7, 0.7);
        // drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);
    }

    // Set the created grid mesh as output
    gridmesh->addVertices(gridvertices);
    meshGridOut.setData(gridmesh);
    
    // TODO (Bonus) Gaussian filter
    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    // Create an editable structured grid with ScalarField2 smoothedField =
    // ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin); Values can be set with
    // smoothedField.setValueAtVertex({0, 0}, 4.2);
    // and read again in the same way as before
    // smoothedField.getValueAtVertex(ij);

    ScalarField2 smoothedGrid = gaussianFilter(grid, nVertPerDim, bBoxMin, bBoxMax);

    // Initialize the output: mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    // get the data values for the regular grid (either original or smoothed)
    std::vector<std::vector<double> > data(nVertPerDim[0], std::vector<double>(nVertPerDim[1]));
    for (int i = 0; i < nVertPerDim[0]; i++) {
        for (int j = 0; j < nVertPerDim[1]; j++) {
            if(propGaussianFilter.get()) { //use the smoothed values
                data[i][j] = smoothedGrid.getValueAtVertex({i, j});
            }
            else { //use the original values
                data[i][j] = grid.getValueAtVertex({i, j});
            }
        }           
    }

    if (propMultiple.get() == 0) {
        // TODO: Draw a single isoline at the specified isovalue (propIsoValue)
        // and color it with the specified color (propIsoColor)
        auto indexBufferIso = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
        algo(propIsoValue, data, bBoxMin, cellSize, propIsoColor, indexBufferIso.get(), vertices);

    }

    else {
        // TODO: Draw the given number (propNumContours) of isolines between
        // the minimum and maximum value
        std::cout << "minValue: " << minValue << " and maxValue: " << maxValue << std::endl;
        std::cout << "numContours: " << propNumContours << std::endl;
    
        for (int i=0 ; i<propNumContours; i++) {
            auto indexBufferIso = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
            double isovalue = minValue + ((i+1) / (double(propNumContours)+1))*(maxValue - minValue);
            
            //if only one color, uncomment the line below
            //algo(isovalue, data, bBoxMin, cellSize, propIsoColor, indexBufferIso.get(), vertices);

            // TODO (Bonus): Use the transfer function property to assign a color
            // The transfer function normalizes the input data and sampling colors
            // from the transfer function assumes normalized input, that means
            // vec4 color = propIsoTransferFunc.get().sample(0.0f);
            // is the color for the minimum value in the data
            // vec4 color = propIsoTransferFunc.get().sample(1.0f);
            // is the color for the maximum value in the data
            vec4 color = propIsoTransferFunc.get().sample((i+1) / (double(propNumContours)+1));
            algo(isovalue, data, bBoxMin, cellSize, color, indexBufferIso.get(), vertices);
        }
    


        
    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour

    mesh->addVertices(vertices);
    meshIsoOut.setData(mesh);
}

float MarchingSquares::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

void MarchingSquares::algo(const double isovalue, std::vector<std::vector<double> > data, const dvec2 bBoxMin, const dvec2 cellSize, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {

    for (int i = 0; i < data.size() -1; i++) {
        for (int j = 0; j < data[0].size() -1 ; j++) { 
            //get the coordinates of the intersections on the cell {i,j}'s edges
            std::vector<vec2> intersections; 
            if ( (data[i][j]>=isovalue) != (data[i][j+1]>=isovalue) ) { //left vertical edge
                double alpha = (isovalue - data[i][j])/(data[i][j+1] - data[i][j]);
                vec2 corner = vec2(bBoxMin[0]+i*cellSize[0], bBoxMin[1]+j*cellSize[1]);
                intersections.push_back(vec2(corner[0], corner[1] + alpha * cellSize[1]));
            }

            if ( (data[i][j]>=isovalue) != (data[i+1][j]>=isovalue) ) {
                double alpha = (isovalue - data[i][j])/(data[i+1][j] - data[i][j]);
                vec2 corner = vec2(bBoxMin[0]+i*cellSize[0], bBoxMin[1]+j*cellSize[1]);
                intersections.push_back(vec2(corner[0] + alpha * cellSize[0], corner[1]));
            }

            if ( (data[i+1][j]>=isovalue) != (data[i+1][j+1]>=isovalue) ) { //right vertical edge
                double alpha = (isovalue - data[i+1][j])/(data[i+1][j+1] - data[i+1][j]);
                vec2 corner = vec2(bBoxMin[0]+(i+1)*cellSize[0], bBoxMin[1]+j*cellSize[1]);
                intersections.push_back( vec2(corner[0], corner[1] + alpha * cellSize[1]));
            }
            
            if ( (data[i][j+1]>=isovalue) != (data[i+1][j+1]>=isovalue) ) {
                double alpha = (isovalue - data[i][j+1])/(data[i+1][j+1] - data[i][j+1]);
                vec2 corner = vec2(bBoxMin[0]+i*cellSize[0], bBoxMin[1]+(j+1)*cellSize[1]);
                intersections.push_back(vec2(corner[0] + alpha * cellSize[0], corner[1]));
            }


            //draw the isolines
            if (intersections.size() == 2) {
                drawLineSegment(intersections[0], intersections[1], color, indexBuffer, vertices);
            }
            else if (intersections.size() == 4) {
                //asymptotic decider
                //we check which intersection vertex on the horizontal edges is closest to the left vertical edge (where intersections[0] is)
                if(propDeciderType == 1) {//if random option choosed
                    if(randomValue(0, 1) < 0.5) {// this segments \\ cells
                        drawLineSegment(intersections[0], intersections[1], color, indexBuffer, vertices);
                        drawLineSegment(intersections[2], intersections[3], color, indexBuffer, vertices);
                    }
                    else { // this segments // cells
                        drawLineSegment(intersections[0], intersections[3], color, indexBuffer, vertices);
                        drawLineSegment(intersections[2], intersections[1], color, indexBuffer, vertices);
                    }
                }
            
                else {
                    if (intersections[1][0] < intersections[3][0]) { //intersections[1] is closest 
                        drawLineSegment(intersections[0], intersections[1], color, indexBuffer, vertices);
                        drawLineSegment(intersections[2], intersections[3], color, indexBuffer, vertices);
                    }
                    else { //intersections[3] is closest
                        drawLineSegment(intersections[0], intersections[3], color, indexBuffer, vertices);
                        drawLineSegment(intersections[2], intersections[1], color, indexBuffer, vertices);
                    }
                }
            }
        }
    }
}

ScalarField2 MarchingSquares::gaussianFilter(const ScalarField2 field, ivec2 nVertPerDim, dvec2 bBoxMin, dvec2 bBoxMax) {
    //gaussian kernel
    double gaussian_kernel[5][5] = {{2, 4, 5, 4, 2},
                                    {4, 9, 12, 9, 4},
                                    {5, 12, 15, 12, 5},
                                    {4, 9, 12, 9, 4},
                                    {2, 4, 5, 4, 2}};
    ScalarField2 smoothedField = ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin);
    for(int i = 0; i<nVertPerDim[0]; i++) {
        for(int j = 0; j<nVertPerDim[1]; j++) {
            double smoothedValue = 0;
            //compute the smooth value using the weights from the gaussian kernel
            //we are padding with 0
            for (int i_kernel = 0; i_kernel < 5; i_kernel++) {
                for (int j_kernel = 0; j_kernel < 5; j_kernel++) {
                    if ((i-2+i_kernel >= 0) && (i-2+i_kernel < nVertPerDim[0]) && (j-2+j_kernel >= 0) && (j-2+j_kernel < nVertPerDim[1])) { //we only compute the weighted value inside the grid (we pad with zero outside the grid)
                        smoothedValue += gaussian_kernel[i_kernel][j_kernel]*field.getValueAtVertex({i-2+i_kernel, j-2+j_kernel}) / 115.0;
                    }   
                }   
            }
            smoothedField.setValueAtVertex({i, j}, smoothedValue);
        }
    }
    return smoothedField;
}

}  // namespace inviwo
