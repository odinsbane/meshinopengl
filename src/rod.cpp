#include <memory>
#include <math.h>
#include "rod.h"

double Line3D::distance(glm::dvec3 &center, glm::dvec3 &direction, double length, const glm::dvec3 &point) {
    glm::dvec3 r = point - center;
    double proj = glm::dot(direction,r);
    double half = length/2.0;
    if(proj<-half){
        //bottom region of end point space
        glm::dvec3 to_bottom(
                center[0] - half*direction[0] - point[0],
                center[1] - half*direction[1] - point[1],
                center[2] - half*direction[2] - point[2]
        );
        return glm::length(to_bottom);
    } else if(proj>half){
        //top region
        glm::dvec3 to_top(
                    center[0] + half*direction[0] - point[0],
                    center[1] + half*direction[1] - point[1],
                    center[2] + half*direction[2] - point[2]
        );
        return glm::length(to_top);
    } else{

        glm::dvec3 normal(r[0] - proj*direction[0], r[1] - proj*direction[1], r[2] - proj*direction[2]);
        return glm::length(normal);
    }
}


bool checkAxis(double a_low, double a_high, double b_low, double b_high){
    return !(b_high<a_low || b_low>a_high);
}

bool Box3D::contains(glm::dvec3 &point){

        return point[0]>=bx && point[0]<=cx
                && point[1]>=by && point[1]<=cy
                && point[2]>=bz && point[2]<=cz;

}

bool Box3D::contains(Box3D &other){
    return other.bx>=bx && other.cx <= cx
            && other.by>=by && other.cy <= cy
            && other.bz>=bz && other.cz <= cz;
};

bool Box3D::intersects(Box3D &other){
    return    checkAxis(bx, cx, other.bx, other.cx)
            && checkAxis(by, cy, other.by, other.cy)
            && checkAxis(bz, cz, other.bz, other.cz);
};


Rod::Rod(double l, double r){
    diameter = r*2;
    length = l;
}

void Rod::applyForce(glm::dvec4 *force) {
    std::lock_guard<std::mutex> lock(mutex);
    forces.push_back(std::unique_ptr<glm::dvec4>(force));
}

void Rod::clearForces(){

    forces.clear();

}

double Rod::closestApproach(Rod &other){
    double dot = glm::dot(direction, other.direction);
    double l_parallel = fabs(dot) * other.length;

    //double[] parallel = new double[]{ dot*direction[0], dot*direction[1], dot*direction[2]};
    glm::dvec3 r = other.position - position;
        double z_norm = glm::dot(direction, r);

        double h = 0.5 * length;

        double o = z_norm + l_parallel * 0.5;
        double p = z_norm - l_parallel * 0.5;

        if (p > h) {
            //no stalk space, pure cap.
            glm::dvec3 pot(position[0] + direction[0] * h, position[1] + direction[1] * h, position[2] + h * direction[2]);
            return Line3D::distance(
                    other.position,
                    other.direction,
                    other.length,
                    pot
            );
        } else if (o > h && p > -h) {

            //top cap space and some stalk space.
            double l_c = (o - h) * other.length / l_parallel;

            double delta = other.length * 0.5 - l_c * 0.5;

            delta = dot < 0 ? -delta : delta;

            glm::dvec3 new_center(
                    other.position[0] + other.direction[0] * delta,
                    other.position[1] + other.direction[1] * delta,
                    other.position[2] + other.direction[2] * delta
            );
            glm::dvec3 pot(position[0] + direction[0] * h, position[1] + direction[1] * h, position[2] + h * direction[2]);
                double cap_distance = Line3D::distance(
                        new_center,
                        other.direction,
                        l_c,
                        pot
                );


                //next step... find region left in stalk.

                //delta is in the opposite direction of the previous delta.
                delta = dot < 0 ? l_c / 2 : -l_c / 2;

            glm::dvec3 stalk_center(
                    other.position[0] + other.direction[0] * delta,
                    other.position[1] + other.direction[1] * delta,
                    other.position[2] + other.direction[2] * delta
            );

            glm::dvec3 r_stalk = stalk_center - position;

                z_norm = glm::dot(direction, r_stalk);

            glm::dvec3 r_perp(
                    r_stalk[0] - z_norm * direction[0],
                    r_stalk[1] - z_norm * direction[1],
                    r_stalk[2] - z_norm * direction[2]
            );

            glm::dvec3 t_perp(
                    other.direction[0] - dot * direction[0],
                    other.direction[1] - dot * direction[1],
                    other.direction[2] - dot * direction[2]
            );

                double m = glm::length(t_perp);
                double stalk_length = other.length - l_c;
                double l_perp = m * stalk_length;

                double stalk_distance;

                if (m == 0) {
                    stalk_distance = glm::length(r_perp);
                } else {
                    //normalize.
                    t_perp = glm::normalize(t_perp);

                    stalk_distance = Line3D::distance(r_perp, t_perp, l_perp, Line3D::origin);
                }


                return stalk_distance < cap_distance ? stalk_distance : cap_distance;

            } else if (o > h && p < -h) {

                //top cap and lengh.
                double l_top = (o - h) * other.length / l_parallel;

                //bottom cap space.
                double l_bottom = -(p + h) * other.length / l_parallel;

                double l_stalk = other.length - l_top - l_bottom;

                double top_delta = other.length * 0.5 - l_top * 0.5;
                double bottom_delta = l_bottom * 0.5 - other.length * 0.5;
                double stalk_delta = (l_bottom - l_top) * 0.5;

                //move in opposite directions if dot is negative.
                top_delta = dot < 0 ? -top_delta : top_delta;
                bottom_delta = dot < 0 ? -bottom_delta : bottom_delta;
                stalk_delta = dot < 0 ? -stalk_delta : stalk_delta;

            glm::dvec3 top_center(
                        other.position[0] + other.direction[0] * top_delta,
                        other.position[1] + other.direction[1] * top_delta,
                        other.position[2] + other.direction[2] * top_delta
            );

            glm::dvec3 pot(position[0] + direction[0] * h, position[1] + direction[1] * h, position[2] + h * direction[2]);
                    double top_distance = Line3D::distance(
                            top_center,
                            other.direction,
                            l_top,
                            pot
                    );


                    //next step... find region left in stalk.

                    //delta is in the opposite direction of the cap region delta.

            glm::dvec3 bottom_center(
                        other.position[0] + other.direction[0] * bottom_delta,
                        other.position[1] + other.direction[1] * bottom_delta,
                        other.position[2] + other.direction[2] * bottom_delta
            );

            glm::dvec3 pot2(position[0] - direction[0] * h, position[1] - direction[1] * h, position[2] - h * direction[2]);
                    double bottom_distance = Line3D::distance(
                            bottom_center,
                            other.direction,
                            l_bottom,
                            pot2
                    );

            glm::dvec3 stalk_center(
                        other.position[0] + other.direction[0] * stalk_delta,
                        other.position[1] + other.direction[1] * stalk_delta,
                        other.position[2] + other.direction[2] * stalk_delta
            );

            glm::dvec3 r_stalk = stalk_center - position;

                    z_norm = glm::dot(direction, r_stalk);

            glm::dvec3 r_perp(
                        r_stalk[0] - z_norm * direction[0],
                        r_stalk[1] - z_norm * direction[1],
                        r_stalk[2] - z_norm * direction[2]
            );


                glm::dvec3 t_perp(
                        other.direction[0] - dot * direction[0],
                        other.direction[1] - dot * direction[1],
                        other.direction[2] - dot * direction[2]
                );

                    double m = glm::length(t_perp);
                    double stalk_distance;
                    double l_perp = m * (l_stalk);
                    if (m == 0) {
                        stalk_distance = glm::length(r_perp);
                    } else {
                        //normalize.
                        t_perp[0] = t_perp[0] / m;
                        t_perp[1] = t_perp[1] / m;
                        t_perp[2] = t_perp[2] / m;

                        stalk_distance = Line3D::distance(r_perp, t_perp, l_perp, Line3D::origin);
                    }
                    double caps = top_distance < bottom_distance ? top_distance : bottom_distance;
                    return caps < stalk_distance ? caps : stalk_distance;

                } else if (o <= h && p > -h) {
                    //only stalk region
            glm::dvec3 r_perp(
                            r[0] - z_norm * direction[0],
                            r[1] - z_norm * direction[1],
                            r[2] - z_norm * direction[2]
            );
            glm::dvec3 t_perp(
                            other.direction[0] - dot * direction[0],
                            other.direction[1] - dot * direction[1],
                            other.direction[2] - dot * direction[2]
            );
                        double m = glm::length(t_perp);
                        double l_perp = m * other.length;

                        if (m == 0) {
                            return glm::length(r_perp);
                        }
                        //normalize.
                        t_perp[0] = t_perp[0] / m;
                        t_perp[1] = t_perp[1] / m;
                        t_perp[2] = t_perp[2] / m;

                        return Line3D::distance(r_perp, t_perp, l_perp, Line3D::origin);

                    } else if (o > -h) {
                        //bottom cap and some stalk.
                        double l_c = -(p + h) * other.length / l_parallel;

                        double delta = other.length * 0.5 - l_c * 0.5;

                        //move in opposite direction as the first case.
                        delta = dot > 0 ? -delta : delta;

            glm::dvec3 new_center(
                                other.position[0] + other.direction[0] * delta,
                                other.position[1] + other.direction[1] * delta,
                                other.position[2] + other.direction[2] * delta
            );

            glm::dvec3 pot(position[0] - direction[0] * h, position[1] - direction[1] * h, position[2] - h * direction[2]);

                            double cap_distance = Line3D::distance(
                                    new_center,
                                    other.direction,
                                    l_c,
                                    pot
                            );


                            //next step... find region left in stalk.

                            //delta is in the opposite direction of the cap region delta.
                            delta = dot > 0 ? l_c / 2 : -l_c / 2;

            glm::dvec3 stalk_center(
                                other.position[0] + other.direction[0] * delta,
                                other.position[1] + other.direction[1] * delta,
                                other.position[2] + other.direction[2] * delta
            );

            glm::dvec3 r_stalk = stalk_center - position;

                            z_norm = glm::dot(direction, r_stalk);

            glm::dvec3 r_perp(
                                r_stalk[0] - z_norm * direction[0],
                                r_stalk[1] - z_norm * direction[1],
                                r_stalk[2] - z_norm * direction[2]
            );


            glm::dvec3 t_perp(
                                other.direction[0] - dot * direction[0],
                                other.direction[1] - dot * direction[1],
                                other.direction[2] - dot * direction[2]
            );

                            double m = glm::length(t_perp);
                            double stalk_distance;
                            double l_perp = m * (other.length - l_c);
                            if (m == 0) {
                                stalk_distance = glm::length(r_perp);
                            } else {
                                //normalize.
                                t_perp[0] = t_perp[0] / m;
                                t_perp[1] = t_perp[1] / m;
                                t_perp[2] = t_perp[2] / m;

                                stalk_distance = Line3D::distance(r_perp, t_perp, l_perp, Line3D::origin);
                            }


                            return stalk_distance < cap_distance ? stalk_distance : cap_distance;


                        } else {
            glm::dvec3 pot(position[0] - direction[0] * h, position[1] - direction[1] * h, position[2] - h * direction[2]);
                            ////no stalk space, pure cap.
                            return Line3D::distance(
                                    other.position,
                                    other.direction,
                                    other.length,
                                    pot
                            );
                        }
}



void Rod::updateBounds(){
    double x1 = position[0] - length*0.5*direction[0];
    double x2 = position[0] + length*0.5*direction[0];
    if(x1>x2){
        bounds.cx = x1 + diameter;
        bounds.bx = x2 - diameter;
    } else{
        bounds.cx = x2 + diameter;
        bounds.bx = x1 - diameter;
    }

    double y1 = position[1] - length*0.5*direction[1];
    double y2 = position[1] + length*0.5*direction[1];
    if(y1>y2){
        bounds.cy = y1 + diameter;
        bounds.by = y2 - diameter;
    } else{
        bounds.cy = y2 + diameter;
        bounds.by = y1 - diameter;
    }

    double z1 = position[2] - length*0.5*direction[2];
    double z2 = position[2] + length*0.5*direction[2];
    if(z1>z2){
        bounds.cz = z1 + diameter;
        bounds.bz = z2 - diameter;
    } else{
        bounds.cz = z2 + diameter;
        bounds.bz = z1 - diameter;
    }
}

Box3D& Rod::getBounds(){
    return bounds;
}