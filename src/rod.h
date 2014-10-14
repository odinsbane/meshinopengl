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
#include <memory>
#include <mutex>
#include <glm/geometric.hpp> //glm::dot
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <array>
#include <math.h>
#include "Constants.h"

#define GLM_FORCE_RADIANS
#include <glm/gtx/rotate_vector.hpp>
namespace Line3D{
    double distance(glm::dvec3 &position, glm::dvec3 &direction, double length, const glm::dvec3 &point);
    const glm::dvec3 origin(0,0,0);
    double closestApproachPosition(glm::dvec3 &center, glm::dvec3 &direction, double length, const glm::dvec3 &point);
    std::vector<double> sphereBounds(glm::dvec3 &center, glm::dvec3 &direction, double length, glm::dvec3 &point, double radius);
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
        Rod(){}
        double length, diameter, stiffness, friction;
        double alpha_longitudinal = 0.5;
        double alpha_perpendicular = 1.0;
        double alpha_rotational = 1.0;
        Rod( double l, double r );
        glm::dvec3 direction;
        glm::dvec3 position;
        std::mutex mutex;
        glm::dvec3 force, torque;
        void applyForce(glm::dvec4* force);
        void clearForces();
        double closestApproach(Rod &other);
        glm::dvec2 intersections(Rod &other);
        double closestApproach(glm::dvec3  &point);
        std::vector<double> getIntersections(glm::dvec3 &point, double radius);
        glm::dvec3 getPoint(double s);
        double collide(Rod &other);
        double prepareForces();
        double update(double dt);
        void updateBounds();
        Box3D& getBounds();


};

class ProxyRod : public Rod{
    Rod* rod;
public:
    ProxyRod(Rod* r) : Rod(r->length, r->diameter*0.5){
        rod=r;
    }
    void applyForce(glm::dvec4* force){
        rod->applyForce(force);
    }
};

class ActinFilament : public Rod{
    std::vector<void*> bound;
public:
    ActinFilament(double l, double r) : Rod(l,r){}
    bool isBound(ActinFilament* f);
};

class MyosinMotor : public Rod{
    std::array<ActinFilament*, 2> bound;
    public:
        double F0, alpha_s, K_m, tau_B;
        MyosinMotor(double l, double r) : Rod(l,r){bound[0] = 0;bound[1] = 0;}
        ActinFilament* getBound(int head);
        void bind(ActinFilament* f, int head);
        bool isBound(int head){return bound[head]!=0;}
        const static int FRONT=0;
        const static int BACK = 1;

    double bind_length;
};


#endif
