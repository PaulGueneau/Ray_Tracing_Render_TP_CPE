
/** TP 5ETI IMI- CPE Lyon - 2015-2016 */



#include "lib/common/error_handling.hpp"

#include "image/image.hpp"
#include "ray_tracing/scene/scene_parameter.hpp"
#include "ray_tracing/scene/ray.hpp"
#include "ray_tracing/primitives/intersection_data.hpp"
#include "ray_tracing/primitives/sphere.hpp"
#include "ray_tracing/primitives/plane.hpp"

#include "ray_tracing/render_engine/render_engine.hpp"

#include <iostream>
#include <fstream>
#include <cmath>



using namespace cpe;

int main(int argc,char *argv[])
{

    std::cout<<"***********************"<<std::endl;
    std::cout<<"run "<<argv[0]<<" with "<<argc-1<<" parameters ... \n"<<std::endl;


    try
    {
        //create empty black image
        int const N_pixels = 300;
        image im(N_pixels);
        im.fill({0,0,0});

        //create a 3D scene for ray tracing
        //***********************************//

        scene_parameter scene;

        //setup camera (center,direction,up , distance to the screen,length of the screen)
        scene.set_camera(camera({0,0,-2},{0,0,1},{0,1,0} , 2.0f,1.2f));

        shading_parameter shading;
        
        //add some 3d sphere (center,radius)
        scene.push_back(new sphere( {-0.5f,1.0f,2.0f} , 0.9f ),
                        material(color(0.5f,0.5f,1.0f),shading,0.5f,0));

        scene.push_back(new sphere( {1.0f,0.8f,1.0f} , 0.2f),
                        material(color(1.0f,0.0f,0.0f),shading,0.9f,1.39336f));

        scene.push_back(new sphere( {1.7f,-0.2f,5.0f} , 0.3f),
                        material(color(0.3f,1.0f,0.0f),shading,0.9f,0.0f));
        scene.push_back(new plane( {0.0f,1.0f,0.0f} , {0.0f,-1.0f,0.0f}),
                        material(color(0.8f,0.8f,0.8f),shading,0.9f,1.2443f));
        //add some 3d plane (position,normal)


        //add a light
        scene.push_back(light( {15.0f,-10.0f,-10.0f} ));
        

        //Render the scene
        //***********************************//
        render(im,scene);
        im.save("image.ppm");

        {   // TEST INTERSECTION  
        int question = 8;
        auto p1 = plane({1.0f,0.0f,0.0f} , {-1.0f,0.0f,0.0f});
        auto r1 = ray({0.,0.,0.},{1.,0.,0.});
        auto i1 = intersection_data({1.,0.,0.},{-1.,0.,0.},{1.});
        auto r2 = ray({0.,0.,0.},{-1.,0.,0.});
        auto s1 = sphere({0.0f,0.0f,0.0f} , 1.0f);
        auto s2 = sphere({2.0f,0.0f,0.0f} , 1.0f);

        intersection_data intersect;

        if(question >= 12 )
        {
            if (s1.intersect(r1, intersect) && intersect != i1)
                throw cpe::exception_cpe("TEST FAILED : sphere intersection",EXCEPTION_PARAMETERS_CPE);
            if(s2.intersect(r1, intersect) && intersect != i1)
                throw cpe::exception_cpe("TEST FAILED : sphere intersection",EXCEPTION_PARAMETERS_CPE);
            if(s2.intersect(r2, intersect))
                 throw cpe::exception_cpe("TEST FAILED : sphere intersection",EXCEPTION_PARAMETERS_CPE);
        }
        if(question >= 8 )
        {
            if (p1.intersect(r1, intersect) && intersect != i1)
                throw cpe::exception_cpe("TEST FAILED : plane intersection",EXCEPTION_PARAMETERS_CPE);
            if(p1.intersect(r2, intersect))
                throw cpe::exception_cpe("TEST FAILED : plane intersection",EXCEPTION_PARAMETERS_CPE);
        }
    }   // FIN TEST INTERSECTION    

    }
    catch(cpe::exception_cpe const& e)
    {


        std::cout<<"**********************************"<<std::endl;
        std::cout<<"**********************************"<<std::endl;
        std::cout<<"Exception found"<<std::endl;
        std::cout<<e.info()<<std::endl;
        std::cout<<"\nSaving backtrace in file backtrace_log.txt"<<std::endl;
        std::cout<<"**********************************"<<std::endl;
        std::cout<<"**********************************"<<std::endl;


        std::ofstream back_file("backtrace_log.txt");
        back_file<<e.info()<<std::endl;
        back_file<<e.info_backtrace();
        back_file.close();

        exit(1);
    }


    std::cout<<"Exit Main"<<std::endl;

    return 0;
}
