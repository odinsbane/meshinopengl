//
//  Ball.h
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#ifndef __ParallelBalls__Rod__
#define __ParallelBalls__Rod__

#include <iostream>
#include <vector>
#include <mutex>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <memory>
#include <mutex>
#include <glm/geometric.hpp> //glm::dot

namespace Line3D{
    double distance(glm::dvec3 &position, glm::dvec3 &direction, double length, const glm::dvec3 &point);
    const glm::dvec3 origin(0,0,0);

}

class Box3D{
public:
    double bx, by, bz, cx, cy, cz;
    bool contains(glm::dvec3  &point);
    bool contains(Box3D &other);
    bool intersects(Box3D &other);


};
class Rod{
    private:
        std::vector<std::unique_ptr<glm::dvec4>> forces;
        Box3D bounds;

    public:
        double length, diameter, stiffness, friction;
        Rod( double l, double r );
        glm::dvec3 direction;
        glm::dvec3 position;
        std::mutex mutex;

        void applyForce(glm::dvec4* force);
        void clearForces();
        double closestApproach(Rod &other);
        glm::dvec2 intersections(Rod &other);
        double closestApproach(glm::vec3  &point);
        glm::dvec2 getIntersections(glm::dvec3 &point, double radius);
        glm::dvec3 getPoint(double s);
        double collide(Rod &other);
        double prepareForces();
        double update(double dt);
        void updateBounds();
        Box3D& getBounds();


};


#endif
