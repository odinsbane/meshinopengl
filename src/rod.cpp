#include <QtCore/qlist.h>
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

double Line3D::closestApproachPosition(glm::dvec3 &center, glm::dvec3 &direction, double length, const glm::dvec3 &point){
    //find which zone point lies in.
    glm::dvec3 r = point - center;
    double proj = glm::dot(direction,r);
    double half = length/2.0;
    if(proj<-half){
        return -half;
    } else if(proj>half){
        return half;
    } else{
        return proj;
    }
}

std::vector<double> Line3D::sphereBounds(glm::dvec3 &center, glm::dvec3 &direction, double length, glm::dvec3 &point, double radius){
    //find which zone point lies in.
    std::vector<double> points;
    glm::dvec3 r = point - center;

    double proj = glm::dot(direction,r);
    double half = length/2.0;

    glm::dvec3 r_perp( r[0] - proj*direction[0], r[1] - proj*direction[1], r[2] - proj*direction[2]);
    double perp = glm::length(r_perp);

    if(perp>radius){
        return points;
    }
    double captured = sqrt(radius*radius - perp*perp);

    double forward = proj + captured;
    double backward = proj - captured;

    if(forward>=-half && forward<=half){
        points.push_back(forward);
    }

    if(backward>=-half && backward<=half){
        points.push_back(backward);
    }

    return points;
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


Rod::Rod(double l, double d){
    diameter = d;
    length = l;
}

void Rod::applyForce(glm::dvec4 *force) {
    std::lock_guard<std::mutex> lock(mutex);
    forces.push_back(std::unique_ptr<glm::dvec4>(force));
}

void Rod::clearForces(){
    std::lock_guard<std::mutex> lock(mutex);
    forces.clear();
    force[0] = 0;
    force[1] = 0;
    force[2] = 0;
    torque[0] = 0;
    torque[1] = 0;
    torque[2] = 0;
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

double Rod::closestApproach(glm::dvec3 &point) {
    return Line3D::distance(position, direction, length, point);
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

/**
* Checks to other rod to see if there is a collision. If there is a collision then a force is applied to both
* rods.
*
* @param other
* return separation distance or -1 if bounding box does not cont
*/

double Rod::collide(Rod &other){
    if(!bounds.intersects(other.bounds)){
        return -1;
    }
    double minimum = 0.5*(other.diameter + diameter);
    double separation = closestApproach(other);
    if(separation<minimum){
        //System.out.println("touching" + separation);
        //touching!
        glm::dvec2 sections = intersections(other);
        glm::dvec3 a = getPoint(sections[0]);
        glm::dvec3 b = other.getPoint(sections[1]);
        glm::dvec3 ab = b - a;
        double mag = glm::length(ab);
        double interference = minimum - separation;
        if(mag==0){
            //TODO
            printf("don't do this!");
        } else{

            double factor = interference*Constants::REPULSION/mag;
            applyForce(new glm::dvec4(
                        -ab[0]*factor,
                        -ab[1]*factor,
                        -ab[2]*factor,
                        sections[0])
            );

            other.applyForce(new glm::dvec4(
                        ab[0]*factor,
                        ab[1]*factor,
                        ab[2]*factor,
                        sections[1]));

        }


    }
    return separation;

}

/**
* Finds the positions of the closest approach for filaments that can be intersecting. returns the coordinates
* as position along the length of the filament.
*
* @param other
* @return
*/
glm::dvec2 Rod::intersections(Rod &other){
    double dot = glm::dot(direction, other.direction);
    double l_parallel = glm::abs(dot)*other.length;

    //double[] parallel = new double[]{ dot*direction[0], dot*direction[1], dot*direction[2]};
    glm::dvec3 r = other.position - position;
    double z_parallel = glm::dot(direction, r);

    double h = 0.5*length;

    double o = z_parallel + l_parallel*0.5;
    double p = z_parallel - l_parallel*0.5;

    glm::dvec2 ret;
    if(p>h){

        //no stalk space, pure cap.
        glm::dvec3 ps(position[0] + direction[0]*h, position[1] + direction[1]*h, position[2] + h*direction[2]);
        ret = glm::dvec2(
                h,
                Line3D::closestApproachPosition(
                    other.position,
                    other.direction,
                    other.length,
                    ps
                )
        );

    } else if(o>h && p>-h){

        //top cap space and some stalk space.
        double l_c = (o - h)*other.length/l_parallel;

        double delta = other.length*0.5 - l_c*0.5;

        delta = dot<0?-delta:delta;

        glm::dvec3 new_center(
            other.position[0] + other.direction[0]*delta,
            other.position[1] + other.direction[1]*delta,
            other.position[2] + other.direction[2]*delta
        );
        glm::dvec3 ps(position[0] + direction[0]*h, position[1] + direction[1]*h, position[2] + h*direction[2]);
        double cap_distance = Line3D::distance(
            new_center,
            other.direction,
            l_c,
            ps
        );

        //next step... find region left in stalk.

        //delta is in the opposite direction of the previous delta.
        delta = dot<0?l_c/2:-l_c/2;

        glm::dvec3 stalk_center(
            other.position[0] + other.direction[0]*delta,
            other.position[1] + other.direction[1]*delta,
            other.position[2] + other.direction[2]*delta
        );

        glm::dvec3 r_stalk = stalk_center - position;

        z_parallel = glm::dot(direction, r_stalk);

        glm::dvec3 r_perp(
            r_stalk[0] - z_parallel*direction[0],
            r_stalk[1] - z_parallel*direction[1],
            r_stalk[2] - z_parallel*direction[2]
        );

        glm::dvec3 t_perp(
            other.direction[0] - dot*direction[0],
            other.direction[1] - dot*direction[1],
            other.direction[2] - dot*direction[2]
        );

        double m = glm::length(t_perp);
        double stalk_length = other.length - l_c;
        double l_perp = m*stalk_length;

        double stalk_distance;

        if(m == 0){
            stalk_distance = glm::length(r_perp);
        } else {
            //normalize.
            t_perp = glm::normalize(t_perp);
            stalk_distance = Line3D::distance(r_perp, t_perp, l_perp, Line3D::origin);
        }

        if(stalk_distance<=cap_distance){
            //if m==0 then use s_stalk is zero, because it can be anywhere along the stalk.
            double s_stalk = 0;
            if(m!=0) {
                double stalk_perp = Line3D::closestApproachPosition(
                    r_perp,
                    t_perp,
                    l_perp,
                    Line3D::origin);

                //scale the stalk perpendicular s to the actual stalk s .
                s_stalk = stalk_perp / l_perp * stalk_length;
            }

            //displace to the other filament coordinate.
            double other_s = s_stalk + delta;

            double stalk_parallel = h - p;

            //moves to the center of the stalk projection.
            double my_s = (h + p)*0.5;
            if(m!=0) {
                if(dot<0){
                    my_s -= s_stalk*stalk_parallel/stalk_length;
                } else {
                    my_s += s_stalk * stalk_parallel / stalk_length;
                }
            }
            ret = glm::dvec2(my_s, other_s);

        } else{
            glm::dvec3 ps = getPoint(h);
            ret = glm::dvec2(
                h,
                Line3D::closestApproachPosition(
                    other.position,
                    other.direction,
                    other.length,
                    ps
                )
            );

        }

    } else if(o>h && p<-h){
    //All three regions are occupied.

        //top cap length.
        double l_top = (o-h)*other.length/l_parallel;

        //bottom cap space.
        double l_bottom = -(p + h)*other.length/l_parallel;

        double l_stalk = other.length - l_top - l_bottom;

        double top_delta = other.length*0.5 - l_top*0.5;
        double bottom_delta = l_bottom*0.5 - other.length*0.5;
        double stalk_delta = (l_bottom - l_top)*0.5;

        //move in opposite directions if dot is negative.
        top_delta = dot<0?-top_delta:top_delta;
        bottom_delta = dot<0?-bottom_delta:bottom_delta;
        stalk_delta = dot<0?-stalk_delta:stalk_delta;

        glm::dvec3 top_center(
            other.position[0] + other.direction[0]*top_delta,
            other.position[1] + other.direction[1]*top_delta,
            other.position[2] + other.direction[2]*top_delta
        );


        glm::dvec3 ps(position[0] + direction[0]*h, position[1] + direction[1]*h, position[2] + h*direction[2]);
        double top_distance = Line3D::distance(
            top_center,
            other.direction,
            l_top,
            ps
        );



        //next step... find region left in stalk.

        //delta is in the opposite direction of the cap region delta.

        glm::dvec3 bottom_center(
                other.position[0] + other.direction[0]*bottom_delta,
                other.position[1] + other.direction[1]*bottom_delta,
                other.position[2] + other.direction[2]*bottom_delta
        );

        glm::dvec3 bs(position[0] - direction[0]*h, position[1] - direction[1]*h, position[2] - h*direction[2]);
        double bottom_distance = Line3D::distance(
            bottom_center,
            other.direction,
            l_bottom,
            bs
        );

        glm::dvec3 stalk_center(
            other.position[0] + other.direction[0]*stalk_delta,
            other.position[1] + other.direction[1]*stalk_delta,
            other.position[2] + other.direction[2]*stalk_delta
        );

        glm::dvec3 r_stalk = stalk_center - position;

        z_parallel = glm::dot(direction, r_stalk);

        glm::dvec3 r_perp(
            r_stalk[0] - z_parallel*direction[0],
            r_stalk[1] - z_parallel*direction[1],
            r_stalk[2] - z_parallel*direction[2]
        );


        glm::dvec3 t_perp(
            other.direction[0] - dot*direction[0],
            other.direction[1] - dot*direction[1],
            other.direction[2] - dot*direction[2]
        );

        double m = glm::length(t_perp);
        double stalk_distance;
        double l_perp = m*(l_stalk);
        if(m == 0){
            stalk_distance = glm::length(r_perp);
        } else {
            //normalize.
            t_perp[0] = t_perp[0] / m;
            t_perp[1] = t_perp[1] / m;
            t_perp[2] = t_perp[2] / m;

            stalk_distance = Line3D::distance(r_perp, t_perp, l_perp, Line3D::origin);
        }
        double caps = top_distance<bottom_distance?top_distance:bottom_distance;
        if(stalk_distance<=caps){

            double s_stalk = 0;
            //leave it zero if it is just a point in perpendicular space.
            if(m!=0) {

                //find the s of the line in perpendicular space
                double stalk_perp = Line3D::closestApproachPosition(
                    r_perp,
                    t_perp,
                    l_perp,
                    Line3D::origin
                );

                //scale the stalk perpendicular to .
                s_stalk = stalk_perp / l_perp * l_stalk;
            }


            double other_s = s_stalk + stalk_delta;
            double stalk_parallel = 2*h;

            double my_s = 0.0;
            if(m!=0) {
                if(dot>0) {
                    my_s += s_stalk * stalk_parallel / l_stalk;
                } else{
                    my_s -= s_stalk * stalk_parallel / l_stalk;
                }
            }
            return glm::dvec2(my_s, other_s);


        }else if(top_distance<bottom_distance){
            glm::dvec3 ps(position[0] + direction[0]*h, position[1] + direction[1]*h, position[2] + h*direction[2]);
            return glm::dvec2(
                h,
                Line3D::closestApproachPosition(
                    other.position,
                    other.direction,
                    other.length,
                    ps
                )
            );

        } else{
            glm::dvec3 ps(position[0] - direction[0]*h, position[1] - direction[1]*h, position[2] - h*direction[2]);
            return glm::dvec2(
                -h,
                Line3D::closestApproachPosition(
                        other.position,
                        other.direction,
                        other.length,
                        ps
                )
            );
        }
    } else if(o<=h && p>-h){
            //only stalk region
        glm::dvec3 r_perp(
            r[0] - z_parallel*direction[0],
            r[1] - z_parallel*direction[1],
            r[2] - z_parallel*direction[2]
        );
        glm::dvec3 t_perp(
            other.direction[0] - dot*direction[0],
            other.direction[1] - dot*direction[1],
            other.direction[2] - dot*direction[2]
        );

        double m = glm::length(t_perp);
        double l_perp = m*other.length;


        //if m==0 then use s_stalk is zero, because it can be anywhere along the stalk.
        double s_stalk = 0;
        if(m!=0) {
            t_perp = glm::normalize(t_perp);

            double stalk_perp = Line3D::closestApproachPosition(
                r_perp,
                t_perp,
                l_perp,
                Line3D::origin
            );

            //scale the stalk perpendicular s to the actual stalk s .
            s_stalk = stalk_perp / l_perp * other.length;
        }

        //displace to the other filament coordinate.
        double other_s = s_stalk;

        double stalk_parallel = o - p;

        //moves to the center of the stalk projection.
        double my_s = (o + p)*0.5;
        if(m!=0) {
            if(dot>0) {
                my_s += s_stalk * stalk_parallel / other.length;
            } else{
                my_s -= s_stalk * stalk_parallel/other.length;
            }
        }
        ret = glm::dvec2(my_s, other_s);

    } else if(o>-h){
        //bottom cap and some stalk.
        double l_c = -(p + h)*other.length/l_parallel;

        double delta = other.length*0.5 - l_c*0.5;

        //move in opposite direction as the first case.
        delta = dot>0?-delta:delta;

        glm::dvec3 new_center(
            other.position[0] + other.direction[0]*delta,
            other.position[1] + other.direction[1]*delta,
            other.position[2] + other.direction[2]*delta
        );
        glm::dvec3 ps(position[0] - direction[0]*h, position[1] - direction[1]*h, position[2] - h*direction[2]);
        double cap_distance = Line3D::distance(
            new_center,
            other.direction,
            l_c,
            ps
        );



        //next step... find region left in stalk.

        //delta is in the opposite direction of the cap region delta.
        delta = dot>0?l_c/2:-l_c/2;

        glm::dvec3 stalk_center(
            other.position[0] + other.direction[0]*delta,
            other.position[1] + other.direction[1]*delta,
            other.position[2] + other.direction[2]*delta
        );

        glm::dvec3 r_stalk = stalk_center - position;

        z_parallel = glm::dot(direction, r_stalk);

        glm::dvec3 r_perp(
            r_stalk[0] - z_parallel*direction[0],
            r_stalk[1] - z_parallel*direction[1],
            r_stalk[2] - z_parallel*direction[2]
        );


        glm::dvec3 t_perp(
            other.direction[0] - dot*direction[0],
            other.direction[1] - dot*direction[1],
            other.direction[2] - dot*direction[2]
        );

        double m = glm::length(t_perp);
        double stalk_length = other.length - l_c;
        double l_perp = m*stalk_length;
        double stalk_distance;
        if(m == 0){
            stalk_distance = glm::length(r_perp);
        } else {
            //normalize.
            t_perp[0] = t_perp[0] / m;
            t_perp[1] = t_perp[1] / m;
            t_perp[2] = t_perp[2] / m;

            stalk_distance = Line3D::distance(r_perp, t_perp, l_perp, Line3D::origin);
        }

        if(stalk_distance<=cap_distance){
            //if m==0 then use s_stalk is zero, because it can be anywhere along the stalk.
            double s_stalk = 0;
            if(m!=0) {
                double stalk_perp = Line3D::closestApproachPosition(
                    r_perp,
                    t_perp,
                    l_perp,
                    Line3D::origin
                );

                //scale the stalk perpendicular s to the actual stalk s .
                s_stalk = stalk_perp / l_perp * stalk_length;
            }

            //displace to the other filament coordinate.
            double other_s = s_stalk + delta;

            double stalk_parallel = o + h;

            //moves to the center of the stalk projection.
            double my_s = (-h + o)*0.5;
            if(m!=0) {
                if(dot<0){
                    my_s += -s_stalk * stalk_parallel / stalk_length;
                } else {
                    my_s += s_stalk * stalk_parallel / stalk_length;
                }
            }
            ret = glm::dvec2(my_s, other_s);
        } else{
            glm::dvec3 ps(position[0] - direction[0]*h, position[1] - direction[1]*h, position[2] - h*direction[2]);
            ret = glm::dvec2(
                -h,
                Line3D::closestApproachPosition(
                        other.position,
                        other.direction,
                        other.length,
                        ps
                )
            );
        }



    } else{
        ////no stalk space, pure cap.
        glm::dvec3 ps(position[0] - direction[0]*h, position[1] - direction[1]*h, position[2] - h*direction[2]);
        ret = glm::dvec2(
            -h,
            Line3D::closestApproachPosition(
                other.position,
                other.direction,
                other.length,
                ps
            )
        );
    }

    return ret;

}


glm::dvec3 Rod::getPoint(double s) {
    return glm::dvec3(position[0] + direction[0]*s, position[1] + direction[1]*s, position[2] + direction[2]*s);
}
glm::dvec3 truncatedCrossProduct(glm::dvec3 &a, glm::dvec4 &b){
    return glm::dvec3(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}
double Rod::prepareForces(){
    std::lock_guard<std::mutex> lock(mutex);

    for(int i = 0; i<forces.size(); i++){
        glm::dvec4 f = *forces[i].get();
        glm::dvec3 t = truncatedCrossProduct(direction, f);
        force[0] += f[0];
        force[1] += f[1];
        force[2] += f[2];
        torque[0] += t[0]*f[3];
        torque[1] += t[1]*f[3];
        torque[2] += t[2]*f[3];
    }
    printf("cum: %15.15e,%15.15e\n", glm::length(force), glm::length(torque));
    //printf("%e,%e,%e\t%e,%e,%e\n", force[0], force[1], force[2], torque[0], torque[1], torque[2]);

    return glm::length(force) + glm::length(torque);
}

double Rod::update(double dt){
    double force_long = glm::dot(force, direction);
    printf("before %15.15e\t", position[0]);
    position[0] = position[0] + dt*(force_long*direction[0]/alpha_longitudinal + (force[0] - force_long*direction[0])/alpha_perpendicular);
    position[1] = position[1] + dt*(force_long*direction[1]/alpha_longitudinal + (force[1] - force_long*direction[1])/alpha_perpendicular);
    position[2] = position[2] + dt*(force_long*direction[2]/alpha_longitudinal + (force[2] - force_long*direction[2])/alpha_perpendicular);
    printf("after: %15.15e\n", position[0]);
    double T = glm::length(torque);

    if(T>0) {
        double omega = T / alpha_rotational;
        double theta = dt * omega;
        torque[0] = torque[0] / T;
        torque[1] = torque[1] / T;
        torque[2] = torque[2] / T;

        direction = glm::rotate(direction, theta, torque);
    }


    updateBounds();
    return T;
}

std::vector<double> Rod::getIntersections(glm::dvec3 &point, double radius) {
    return Line3D::sphereBounds(position, direction, length, point, radius);
}

void MyosinMotor::bind(ActinFilament* f, int head){
    bound[head] = f;
}

ActinFilament *MyosinMotor::getBound(int head) {
    return bound[head];
}

bool ActinFilament::isBound(ActinFilament *filament) {

    for(void* o: bound){
        if(filament==o){
            return true;
        }
    }
    return false;
}

void ActinFilament::bind(ActinFilament *f){
    bound.push_back(f);
}

void ActinFilament::unbind(ActinFilament* f){
    for(auto g = bound.begin(); g!=bound.end(); g++){
        if(*g==f){
            bound.erase(g);
        }
        return;
    }
}