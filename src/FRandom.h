#include <random>
#include <glm/vec3.hpp>
#ifndef _FRANDOM_
#define _FRANDOM_
    class FRandom{
        private:
            std::mt19937_64 engine;
            double max = 1.0*std::mt19937_64::max();
        public:
            FRandom(long seed){
                engine.seed(seed);
            }

            double nextDouble(){
                return engine()/max;
            }
            glm::dvec3 randomDirection(){
                glm::dvec3 v(nextDouble() - 0.5, nextDouble() - 0.5, nextDouble() - 0.5);
                while(v.length()==0){
                    v[0]=nextDouble();
                    v[1]=nextDouble();
                    v[2]=nextDouble();
                }
                return glm::normalize(v);
            }

    };

#endif