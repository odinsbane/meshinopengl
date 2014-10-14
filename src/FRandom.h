#include <random>
#include <glm/vec3.hpp>  //vec3
#include <glm/geometric.hpp> //normalize
#ifndef _FRANDOM_
#define _FRANDOM_
    class FRandom{
        private:
            std::mt19937_64 engine;
            double imax = 1.0/(1.0*std::mt19937_64::max());
        public:
            FRandom(long seed){
                engine.seed(seed);
            }
            FRandom(){
                engine.seed(std::random_device{}());
            }

            double nextDouble(){
                return engine()*imax;
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
            int nextInt(int range){
                return int(nextDouble()*range);
            }
    };

#endif