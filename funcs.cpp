#include "geom.hpp"
#include <cmath>

plane_t find_plane(const triangle& tr)
{
    Geometry::vector s1 = sub(tr.vertices[1], tr.vertices[0]);
    Geometry::vector s2 = sub(tr.vertices[2], tr.vertices[0]);
    plane_t pl;
    pl.normal = cross(s1, s2);  
    pl.a = pl.normal.x;
    pl.b = pl.normal.y;
    pl.c = pl.normal.z;
    pl.d = (pl.a*tr.vertices[0].x + pl.b*tr.vertices[0].y + pl.c*tr.vertices[0].z) * (-1);
    return pl;
}

bool planes_intersect(const plane_t& pl1, const plane_t& pl2, line& line) //GCT page 529
{
    Geometry::vector d = cross(pl1.normal, pl2.normal);
    if (length(d) <= epsilon)
    {
    return false;
    }
    line.direction = d;

    double s1 = -pl1.d;
    double s2 = -pl2.d;
    double n1n2dot = dot(pl1.normal, pl2.normal);
    double n1normsqr = dot(pl1.normal, pl1.normal);
    double n2normsqr = dot(pl2.normal, pl2.normal);
    double a = (s2 * n1n2dot - s1 * n2normsqr) / (n1n2dot * n1n2dot - n1normsqr * n2normsqr);
    double b = (s1 * n1n2dot - s2 * n1normsqr) / (n1n2dot * n1n2dot - n1normsqr * n2normsqr);
    line.point = add(multiply(pl1.normal, a), multiply(pl2.normal, b));
    return true;
}

Geometry::vector cross(const Geometry::vector& a, const Geometry::vector& b) //векторное произведение, тоже вектор
{
    Geometry::vector n;
    n.x = a.y * b.z - a.z * b.y;
    n.y = a.z * b.x - a.x * b.z;
    n.z = a.x * b.y - a.y * b.x;
    return n;
}

double dot(const Geometry::vector& a, const Geometry::vector& b)
{
    double ans = a.x*b.x + a.y*b.y + a.z*b.z;
    return ans;
}

double length(const Geometry::vector& a) 
{
    double length = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    return length;
}

Geometry::vector multiply(const Geometry::vector& a, double value) 
{
    Geometry::vector res;
    res.x = a.x * value;
    res.y = a.y * value;
    res.z = a.z * value;
    return res;
}

Geometry::vector add(const Geometry::vector& a, const Geometry::vector& b) //вектор
{
    Geometry::vector res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    return res;
} 

Geometry::vector sub(const Geometry::vector& a, const Geometry::vector& b)
{
    Geometry::vector res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    return res;
}

State parallel(const Geometry::vector& v1, const Geometry::vector& v2)
{
    if(std::abs(length (cross(v1,v2)) ) < epsilon) {return parall;}
    else {return non_parall;}
}

State parallel(const plane_t& pl1, const plane_t& pl2)
{
   if (parallel(pl1.normal, pl2.normal) != parall) {return intersecting;}
   else
    {
    if (std::abs(pl1.d-pl2.d) < epsilon) {return coincident;}
    else {return parall;}
    }
}

double dist(const Geometry::vector& p, const plane_t& plane) //GCT page 376
{
    double l_normal = length(plane.normal);
    auto n = multiply( plane.normal, 1/l_normal);
    double r = dot(n, p) + plane.d / l_normal;
    return r;
} 

bool intersection_3D(const triangle& tr1, const triangle& tr2) 
{
    std::vector<Geometry::vector> projections1;
    plane_t plane1 = find_plane(tr1);
    plane_t plane2 = find_plane(tr2);
    line L; //линия пересечения
    planes_intersect(plane1, plane2, L); // L = P + td

    std::vector<double> interval1 = get_interval_3D(tr1, plane2, L);
    std::vector<double> interval2 = get_interval_3D(tr2, plane1, L);

    return overlap(interval1, interval2);
    //если dist0=dist1=dist2=0, то треугольники компланарны 
}

std::vector<double> get_interval_3D(const triangle& tr, const plane_t& pl, const line& L)
{
    double dist0 = dist(tr.vertices[0], pl);
    double dist1 = dist(tr.vertices[1], pl);
    double dist2 = dist(tr.vertices[2], pl);
    Geometry::vector together1;
    Geometry::vector together2;
    Geometry::vector another;
    double dist_another = 0;
    double dist_together1 = 0;
    double dist_together2 = 0;
    std::vector<double> empty;

    if ((dist0 * dist1) > 0)
    { if((dist2*dist1) > 0) {return empty;}
        
    else    {another = tr.vertices[2];
            together1 = tr.vertices[0];
            together2 = tr.vertices[1];
            dist_another = dist2;
            dist_together1 = dist0;
             dist_together2 = dist1;}
              
    }
    else 
    { if ((dist0 * dist2) > 0)
        {another = tr.vertices[1];
        together1 = tr.vertices[0];
        together2 = tr.vertices[2];
        dist_another = dist1;
        dist_together1 = dist0;
         dist_together2 = dist2;}
     else 
        {another = tr.vertices[0];
        together1 = tr.vertices[1];
        together2 = tr.vertices[2];
        dist_another = dist0;
        dist_together1 = dist1;
        dist_together2 = dist2;}
    }

    double p_together1 = dot(L.direction, sub(together1, L.point));
    double p_together2 = dot(L.direction, sub(together2, L.point));
    double p_another = dot(L.direction, sub(another, L.point));

    double t1 = p_together1 + (p_another - p_together1) * 
                dist_together1 / (dist_together1 - dist_another);
    double t2 = p_together2 + (p_another - p_together2) * 
                dist_together2 / (dist_together2 - dist_another);

    std::vector<double> interval;
    double min = 0; double max = 0;
    if (t2 > t1) {min = t1; max = t2;}
    else {min = t1; max = t2;}
    interval.push_back(min);
    interval.push_back(max);
    return interval;
} 

bool overlap(const std::vector<double> interval1, const std::vector<double> interval2)
{
    if ((interval2[0] > interval1[1]) || (interval1[0] > interval2[1]) ) {return false;}
    else {return true;}
}

bool same_side(const triangle& tr, const plane_t& pl)
{
    double dist0 = dist(tr.vertices[0], pl);
    double dist1 = dist(tr.vertices[1], pl);
    double dist2 = dist(tr.vertices[2], pl);

    if ( (dist0>0 && dist1>0 && dist2>0) ||
         (dist0<0 && dist1<0 && dist2<0) ) {return true;}
    else {return false;}
}

std::vector<double> get_interval_2D(const triangle& tr, const Geometry::vector D)
{
    double min = dot(D, tr.vertices[0]);
    double max = min;
    for (int i = 1; i < 3; ++i)
    {
    double value = dot(D, tr.vertices[i]);
    if (value < min) {min = value;}
    else if (value > max) {max = value;}
    }
    std::vector<double> interval;
    interval.push_back(min);
    interval.push_back(max);
    return interval;
}

Geometry::vector Perp(Geometry::vector a)
{
    Geometry::vector b;
    if (abs(a.x) > epsilon)
    {
        b.y = 1;
        b.z = 1;
        b.x = (-1)*(a.y * b.y + a.z * b.z) / (a.x);
    }
    else if (abs(a.y) > epsilon)
    {
        b.x = 1;
        b.z = 1;
        b.y = (-1)*(a.x * b.x + a.z * b.z) / (a.y);
    }
    else if (abs(a.z) > epsilon)
    {
        b.x = 1;
        b.y = 1;
        b.z = (-1)*(a.x * b.x + a.y * b.y) / (a.z); //
    }
    return b;
}

bool intersection_2D(const triangle& tr1, const triangle& tr2)
{
    for (int i0=0, i1=2; i0 < 3; i1=i0, i0++)
    {
    auto D = Perp(sub(tr1.vertices[i0], tr1.vertices[i1]));
    auto interval1 = get_interval_2D(tr1, D);
    auto interval2 = get_interval_2D(tr2, D);
    if (overlap(interval1, interval2) == false) {return false;}
    }
    for (int i0=0, i1=2; i0 < 3; i1=i0, i0++)
    {
    auto D = Perp(sub(tr2.vertices[i0], tr2.vertices[i1]));
    auto interval1 = get_interval_2D(tr1, D);
    auto interval2 = get_interval_2D(tr2, D);
    if (overlap(interval1, interval2) == false) {return false;}
    }

    return true;
}





