/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>

namespace inviwo {

// TODO: Implement a single integration step here
dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position) {
//     Access the vector field with vectorField.interpolate(...)
    dvec2 velocity = vectorField.interpolate(position);
    return velocity;        
}

dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position, const double step, const bool normalizeVelocity, const vec2 vector0) {
    //we compute v1, v2, v3 and v4 for 4th order Runge-Kutta
    dvec2 v1;
    if ((vector0[0] == 0)&&(vector0[1] == 0)) //if vector0 != (0,0), we use it instead of the vector field at the first position
        v1 = vectorField.interpolate(position);
    else
        v1 = vec2(vector0[0], vector0[1]);
    if (normalizeVelocity) v1 = glm::normalize(v1);
    dvec2 v2 = vectorField.interpolate(position + step / 2.0 * v1);
    if (normalizeVelocity) v2 = glm::normalize(v2);
    dvec2 v3 = vectorField.interpolate(position + step/2.0 * v2);
    if (normalizeVelocity) v3 = glm::normalize(v3);
    dvec2 v4 = vectorField.interpolate(position + step * v3);
    if (normalizeVelocity) v4 = glm::normalize(v4);
    //compute velocity
    return v1 / 6.0 + v2 / 3.0 + v3 / 3.0 + v4 / 6.0;
}


void Integrator::drawPoint(const dvec2& p, const vec4& color, IndexBufferRAM* indexBuffer,
                           std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(p[0], p[1], 0), vec3(0, 0, 1), vec3(p[0], p[1], 0), color});
}

// Alias for draw point
void Integrator::drawNextPointInPolyline(const dvec2& p, const vec4& color,
                                         IndexBufferRAM* indexBuffer,
                                         std::vector<BasicMesh::Vertex>& vertices) {
    Integrator::drawPoint(p, color, indexBuffer, vertices);
}

void Integrator::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                                 IndexBufferRAM* indexBuffer,
                                 std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

std::vector<dvec2> Integrator::computeEquidistantStreamline(const dvec2& startPoint, const VectorField2& vectorField, const double stepSize, const int kernelSize) {
    std::vector<dvec2> streamline;
    streamline.push_back(startPoint);
    //get the bounding box of the vector field
    dvec2 BBoxMin_ = vectorField.getBBoxMin();
    dvec2 BBoxMax_ = vectorField.getBBoxMax();

    //compute streamline in forward direction
    dvec2 currentPoint = vec2(startPoint[0], startPoint[1]);
    double step = stepSize; //forward direction: step = +stepSize
    for (int i = 0; i < kernelSize / 2.0; i++) {
        dvec2 normalizedVelocity = glm::normalize(RK4(vectorField, currentPoint, step));
        dvec2 nextPoint = currentPoint + step * normalizedVelocity;

        if ((nextPoint[0] <= BBoxMax_[0]) && (nextPoint[0] >= BBoxMin_[0])&&((nextPoint[1] <= BBoxMax_[1]))&&
            (nextPoint[1] >= BBoxMin_[1])&&(glm::length(normalizedVelocity) != 0)) {
            currentPoint = nextPoint;
            streamline.push_back(currentPoint);
        }
        else {
            break;
        }
    }

    //compute streamline in backward direction
    currentPoint = vec2(startPoint[0], startPoint[1]);
    step = - stepSize; //backward direction: step = -stepSize
    for (int i = 0; i < kernelSize / 2.0; i++) {
        dvec2 normalizedVelocity = glm::normalize(RK4(vectorField, currentPoint, step));
        dvec2 nextPoint = currentPoint + step * normalizedVelocity;

        if ((nextPoint[0] <= BBoxMax_[0]) && (nextPoint[0] >= BBoxMin_[0])&&((nextPoint[1] <= BBoxMax_[1]))&&
            (nextPoint[1] >= BBoxMin_[1])&&(glm::length(normalizedVelocity) != 0)) {
            currentPoint = nextPoint;
            streamline.push_back(currentPoint);
        }
        else {
            break;
        }
    }

    return streamline;
}

std::deque<dvec2> Integrator::computeEquidistantMaxStreamline(const dvec2& startPoint, const VectorField2& vectorField, const double stepSize, const double minVelocity) { //same as below, but get the max length of the streamline
    double eps = 0.000000001;
    std::deque<dvec2> streamline;
    //get the bounding box of the vector field
    dvec2 BBoxMin_ = vectorField.getBBoxMin();
    dvec2 BBoxMax_ = vectorField.getBBoxMax();

    dvec2 currentPoint = vec2(startPoint[0], startPoint[1]);    
    //compute streamline in backward direction
    bool uncompleteStreamline = true;
    double step = - stepSize; //backward direction: step = -stepSize
    while (uncompleteStreamline) {
        dvec2 velocity = RK4(vectorField, currentPoint, step);
        dvec2 normalizedVelocity = glm::normalize(velocity);
        dvec2 nextPoint = currentPoint + step * normalizedVelocity;
        uncompleteStreamline = (nextPoint[0] <= BBoxMax_[0]) && (nextPoint[0] >= BBoxMin_[0] )&& ((nextPoint[1] <= BBoxMax_[1])) && (nextPoint[1] >= BBoxMin_[1]) 
                                && (glm::length(normalizedVelocity) != 0) && (glm::length(nextPoint - startPoint) >= eps) && (glm::length(velocity) >= minVelocity);
        if (uncompleteStreamline) {
            currentPoint = nextPoint; 
            streamline.push_front(currentPoint);
        }
    }

    streamline.push_back(startPoint);

    //compute streamline in forward direction
    uncompleteStreamline = true;
    currentPoint = vec2(startPoint[0], startPoint[1]);
    step = stepSize; //forward direction: step = +stepSize
    while (uncompleteStreamline) {
        dvec2 normalizedVelocity = glm::normalize(RK4(vectorField, currentPoint, step));
        dvec2 nextPoint = currentPoint + step * normalizedVelocity;
        uncompleteStreamline = (nextPoint[0] <= BBoxMax_[0]) && (nextPoint[0] >= BBoxMin_[0]) && ((nextPoint[1] <= BBoxMax_[1])) && (nextPoint[1] >= BBoxMin_[1]) 
                                && (glm::length(normalizedVelocity) != 0) && (glm::length(nextPoint - startPoint) >= eps);
        if (uncompleteStreamline) {
            currentPoint = nextPoint;
            streamline.push_back(currentPoint);
        }
    }


    return streamline;

}

}  // namespace inviwo
