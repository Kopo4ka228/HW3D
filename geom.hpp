//Вход: На вход поступает 0<N<1 000 000 наборов точек, представяляющих трехмерные треугольники. 
//Задача: Распечатать номера треугольников, которые пересекаются с каким-либо

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <optional>

const double epsilon = 0.0001; //погрешность для вычислений

enum State {parall, coplanar, coincident, intersecting, non_parall}; 

namespace Geometry
{
    struct vector
    {
        double x, y, z;
        vector() {x = NAN, y = NAN, z = NAN;}
        vector(double val1, double val2, double val3) {x = val1, y = val2, z = val3;}
        void print() const {std::cout << "[" << x << "," << y << "," << z << "]";}
    };
}

struct point
{
    double x = NAN, y = NAN, z = NAN; //координаты точки

    void print() {std::cout << "x = " << x << "y = " << y;}
    bool equal(const point& ref)
    {
        return (std::abs(x-ref.x) < epsilon) && (std::abs(y-ref.y) < epsilon)
            && (std::abs(z-ref.z) < epsilon);
    }

};

struct plane_t
{
    double a = NAN, b = NAN, c = NAN, d = NAN; // ax+by+cz+d=0 - уравнение плоскости
    Geometry::vector normal; //нормальный вектор
    void print() const
    {
    const char* signb = (b > 0.0) ? "+" : "";
    const char* signc = (c > 0.0) ? "+" : "";
    const char* signd = (d > 0.0) ? "+" : "";
    std::cout << a << "x" << signb << b << "y" << signc << c 
              << "z" << signd << d << "d" << " = 0" << std::endl;
    }
};

struct triangle
{
    std::vector<Geometry::vector> vertices;
    plane_t own_plane;
    std::vector<Geometry::vector> Edges;
    
    void find_plane(const triangle& tr); //вычислить уравнение плоскости
    bool intersect(); //узнать, пересекается ли с другим треугольником
    triangle(Geometry::vector a1, Geometry::vector a2, Geometry::vector a3)
    {
    vertices.push_back(a1);
    vertices.push_back(a2);
    vertices.push_back(a3);
    }
    void print() const
    {
    std::cout << "{";
    vertices[0].print(); std::cout << ", ";
    vertices[1].print(); std::cout << ", ";
    vertices[2].print();
    std::cout << "}";
    }
};

struct line
{
    Geometry::vector direction; //направляющий вектор
    Geometry::vector point; //точка прямой
    void print() const
    {
    std::cout << "vector: ";
    direction.print();
    std::cout << std::endl << "point: ";
    point.print();
    std::cout << std::endl;
    }
};


plane_t find_plane(const triangle&); // +

double dot(const Geometry::vector&, const Geometry::vector&); //+
double length(const Geometry::vector&); //+
double dist(const Geometry::vector&, const plane_t&); //+

Geometry::vector multiply(const Geometry::vector& a, double value); //+
Geometry::vector add(const Geometry::vector&, const Geometry::vector&); //++
Geometry::vector sub(const Geometry::vector&, const Geometry::vector&); //+
Geometry::vector cross(const Geometry::vector&, const Geometry::vector&); //+
Geometry::vector Perp(Geometry::vector); //+

bool overlap(const std::vector<double>, const std::vector<double>); //+
bool same_side(const triangle&, const plane_t&);
bool intersection_2D(const triangle&, const triangle&); //+
State parallel(const plane_t& pl1, const plane_t& pl2); //+
State parallel(const Geometry::vector&, const Geometry::vector&); //+
bool planes_intersect(const plane_t&, const plane_t&, line&); //+

std::vector<double> get_interval_2D(const triangle&, const Geometry::vector); //+
std::vector<double> get_interval_3D(const triangle&, const plane_t&, const line&);
