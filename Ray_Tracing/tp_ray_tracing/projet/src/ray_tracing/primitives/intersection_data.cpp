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

#include <cmath>
#include "intersection_data.hpp"

namespace cpe
{

intersection_data::intersection_data()
    :position(vec3()),normal(vec3(0,0,1)),relative(-1.0f)
{}

intersection_data::intersection_data(vec3 const& position_param,
                      vec3 const& normal_param,
                      float relative_param)
    :position(position_param),normal(normal_param),relative(relative_param)
{}

void intersection_data::set(vec3 const& position_param,
         vec3 const& normal_param,
         float relative_param)
{
    position = position_param;
    normal   = normal_param;
    relative = relative_param;
}


bool operator==(intersection_data const& intersection1, intersection_data const& intersection2)
{
    return norm(intersection1.position-intersection2.position) < 1e-5
    && norm(intersection1.normal-intersection2.normal) < 1e-5
    && std::fabs(intersection1.relative - intersection2.relative) < 1e-5;
}
bool operator!=(intersection_data const& intersection1, intersection_data const& intersection2)
{
    return !(intersection1 == intersection2);
}

std::ostream& operator<<(std::ostream& stream,intersection_data const& intersection)
{
    stream<<intersection.position<<" "<<intersection.normal<<" "<<intersection.relative;
    return stream;
}

}
