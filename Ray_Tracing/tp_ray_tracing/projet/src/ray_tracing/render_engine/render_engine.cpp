/*
**    TP CPE Lyon
**    Copyright (C) 2015 Damien Rohmer
**
**    This program is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "render_engine.hpp"

#include "image/image.hpp"
#include "ray_tracing/scene/scene_parameter.hpp"
#include "ray_tracing/scene/ray.hpp"
#include "ray_tracing/primitives/intersection_data.hpp"
#include "ray_tracing/scene/anti_aliasing_table.hpp"

#include <cmath>

namespace cpe
{


void render(image& im,scene_parameter const& scene)
{
    // **************************************************************************** //
    //
    // Current Algorithm:
    //
    // For all pixels of the image
    //    Generate ray from the camera toward in the direction of the current pixel
    //    Compute associated color (ray tracing algorithm)
    //    Set the color to the current pixel
    //
    // **************************************************************************** //

    camera const& cam = scene.get_camera();

    int const Nx = im.Nx();
    int const Ny = im.Ny();
    int n_Sample = 5;
    float sigma = 0.8f;
    anti_aliasing_table aai = anti_aliasing_table(n_Sample,1.0f,sigma);

    // loop over all the pixels of the image
    for(int kx=0 ; kx<Nx ; ++kx)
    {
        float const u = static_cast<float>(kx)/(Nx-1);
        for(int ky=0 ; ky<Ny ; ++ky)
        {
            float const v = static_cast<float>(ky)/(Ny-1);
            //float const w = aai.weight(dx,dy);
            // generate ray and compute color
            //ray const r = ray_generator(cam,u,v);

            color mean_color_aa; 
            for(int dx=0; dx<n_Sample;++dx){
                for(int dy=0 ; dy<n_Sample ; ++dy){
                    float const du = aai.displacement(dx)/(Nx-1); //Nx is the size in pixel in x direction
                    float const dv = aai.displacement(dy)/(Ny-1); //Ny is the size in pixel in y direction
                    ray const r_aa = ray_generator(cam,u+du,v+dv);
                    //std::cout<< mean_color_aa << std::endl;
                    mean_color_aa += aai.weight(dx,dy) * ray_trace(r_aa,scene); //F is a function of (u,v)
                }
            }

            mean_color_aa = mean_color_aa;
            
            im({kx,ky}) = mean_color_aa;
        }
    }

}


ray ray_generator(camera const& cam,float const u,float const v)
{
    // position of the sample on the screen in 3D
    vec3 const p_screen = screen_position(cam,u,v);

    // vector "camera center" to "screen position"
    vec3 const d = p_screen-cam.center();

    // compute the ray
    ray const r(cam.center(),d);

    return r;
}

float signe(float a) {
    if (a > 0) return 1;
    if (a < 0) return -1;
    return 0;
}

color ray_trace(ray const& r,scene_parameter const& scene)
{
    // ********************************************* //
    // ********************************************* //
    //
    // TO DO: Calculer la couleur associee a un rayon dans la scene 3D
    //
    // ********************************************* //
    // ********************************************* //
    //Le code suivant affecte la couleur de base de la premiere intersection trouvee
    //Ce code ne gere pas les reflections.

    int n_reflexion = 3;
    float n_air = 1.0003;
    color sum_color;
    intersection_data intersection; //current intersection
    int intersected_primitive = 0;  //current index of intersected primitive
    ray rayon = r;

    for (int k = 0;k < n_reflexion; k++) {

        bool const is_intersection = compute_intersection(rayon,scene,intersection,intersected_primitive);
        if(is_intersection){  //if the ray intersects a primitive
            material const& mat = scene.get_material(intersected_primitive);
            ray refraction = r;
            /*
            if (k == 0) { // refraction une seule fois

                if (rayon.is_refracted(n_air,mat.refraction(), r.u(),intersection.normal)){
                    float alpha = norm(cross(r.u(),intersection.normal))/tan(asin((n_air/mat.refraction())*norm(cross(r.u(),intersection.normal)))) * signe(dot(refraction.u(),intersection.normal))- dot(refraction.u(),intersection.normal);
                    vec3 direction = rayon.u() + alpha*intersection.normal;
                    refraction(intersection.position, direction);
                    refraction.offset();
                }
                intersection_data intersection_refracted;
                int intersected_primitive_refracted = 0; 
                bool const is_intersection_reflected = compute_intersection(refraction,scene,intersection_refracted,intersected_primitive_refracted);
                sum_color += compute_illumination(mat,intersection_refracted,scene);
            }
            */


            
            ray reflection(intersection.position, rayon.u()-2*dot(rayon.u(),intersection.normal) * intersection.normal );
            reflection.offset();
            intersection_data intersection_reflected;
            int intersected_primitive_reflected = 0;  //current index of intersected reflected primitive

            bool const is_intersection_reflected = compute_intersection(reflection,scene,intersection_reflected,intersected_primitive_reflected);
            sum_color += (mat.reflection()/(1+k))*compute_illumination(mat,intersection,scene);
            rayon = reflection;
            intersection = intersection_reflected;
            intersected_primitive = intersected_primitive_reflected;

        } 
    }

    return sum_color;

}


bool compute_intersection(ray const& r,
                          scene_parameter const& scene,
                          intersection_data& intersection,
                          int& index_intersected_primitive)
{
    // ********************************************* //
    // ********************************************* //
    //
    // TO DO: Calculer l'intersection valide la plus proche
    //        Doit retourner vrai/faux si une intersection est trouvee ou non
    //        Doit mettre a jour intersection et index_intersected_primitive avec l'intersection la plus proche
    //
    // ********************************************* //
    // ********************************************* //

    //Le code suivant renvoie la premiere intersection valide que l'on trouve dans l'ordre de stockage du vecteur
    //Ce code est a modifier afin de trouver la 1ere intersection dans l'espace 3D le long du rayon.

    int const N_primitive = scene.size_primitive();
    bool found_intersection = false;

    float intersection_var_relatif = 100000.0f;

    intersection_data intersection_var;

    for(int k = 0;k<N_primitive;k++)
    {
        primitive_basic const & primitive = scene.get_primitive(k);
        bool is_intersection = primitive.intersect(r,intersection_var);
        if(is_intersection)
        {    found_intersection = true;
            if (intersection_var.relative < intersection_var_relatif) {
                index_intersected_primitive = k;
                intersection = intersection_var;
                intersection_var_relatif = intersection_var.relative;
            }
        }
    }
    return found_intersection;
}


bool is_in_shadow(vec3 const& p,vec3 const& p_light, scene_parameter const& scene)
{ 
    ray rayon_shadow = ray(p,normalized(p_light-p));
    rayon_shadow.offset();

    //rayon_shadow.offset(1e10-5);
    int index_intersected_primitive;
    intersection_data intersection_var;
    if (compute_intersection(rayon_shadow, scene, intersection_var,index_intersected_primitive)) {
        if (norm(p_light-p) > intersection_var.relative) {
            return true;
        }
    }
    

    
    // ********************************************* //
    //
    // TO DO: Calculer si le point p est dans l'ombre de la lumiere situee en p_light
    //
    // ********************************************* //
    return false;

}



color compute_illumination(material const& mat,intersection_data const& intersection,scene_parameter const& scene)
{
    // ********************************************* //
    //
    // TO DO: Mettre en place le calcul d'illumination correct
    //
    //   Pour toutes les lumieres
    //     Si point dans l'ombre
    //       Ajouter a couleur courante une faible part de la couleur du materiau
    //          ex. c += mat.color_object()/10.0f;
    //     Sinon
    //       Calculer illumination au point courant en fonction de la position
    //          de la lumiere et de la camera
    //       Ajouter a couleur courante la contribution suivante:
    //           puissance de la lumiere (L.power)
    //         X couleur de la lumiere (L.c)
    //         X valeur de l'illumination
    //
    // ********************************************* //

    color c;

    vec3 const& p0 = intersection.position;

    int const N_light = scene.size_light();
    for(int k=0; k<N_light ; ++k)
    {
        light const& L = scene.get_light(k);
        bool const shadow = is_in_shadow(p0,L.p,scene);
        if(shadow==false)
            c +=  L.power*L.c*compute_shading(mat.shading(),mat.color_object(),p0,L.p,intersection.normal,scene.get_camera().center());
           // mat.color_object();
        else
        {
            c += mat.color_object()/10.0f;
        }
        
    }

    return c;


}





}


