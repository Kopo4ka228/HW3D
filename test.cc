#include "geom.hpp"
#include <fstream>
#include <ctime>
#include <set>

int main()
{
    std::ifstream fin("input.txt");
    std::ofstream out("output.txt");
    if (!fin.is_open()) {std::cout << "file wrong";}
    int N = 0;
    fin >> N;
    std::vector<triangle> arr;
    unsigned int start_time1 =  clock();
    std::set<int> nums_tr; //множество всех номеров треугольников
    std::set<int> nums_inter; //множество пересекающихся
    for (int j = 0; j < N; j++) {
        nums_tr.insert(j);
        std::vector<Geometry::vector> Points;
        for (int i=0; i < 3; i++) {
            double x,y,z;
            fin >> x >> y >> z;
            Geometry::vector P(x,y,z); 
            Points.push_back(P);
        }
        triangle tr(Points[0], Points[1], Points[2]);
        tr.own_plane = find_plane(tr);
        arr.push_back(tr);
    }
    
    unsigned int end_time1 = clock();
    std::cout << end_time1 << "ms" << std::endl;
    
    unsigned int start_time2 = clock();
    for (int j = 0; j < N-1; j++) {
    if (nums_tr.count(j) == 0) continue;
    auto tr0 = arr[j];
    
    for (int i = j+1; i < N; i++)
   {
        auto tr1 = arr[i]; 

        if (same_side(tr1, tr0.own_plane)) {continue;}
        
        State r = parallel(tr1.own_plane, tr0.own_plane);
        switch(r)
        {
            case parall: break;
            case coincident: { 
                if(intersection_2D(tr0,tr1)) 
                {nums_tr.erase(i); nums_inter.insert(i);break;}
                else {break;}
             }
            case intersecting: { 
                if (same_side(tr0, tr1.own_plane)) {break;}
                else
                {
                    line L;
                    planes_intersect(tr0.own_plane, tr1.own_plane, L);
                    auto interval0 = get_interval_3D(tr0, tr1.own_plane, L);
                    auto interval1 = get_interval_3D(tr1, tr0.own_plane, L);
                    if (overlap(interval0,interval1)) {nums_tr.erase(i); nums_inter.insert(i);break;}
                    else {break;}
                }
                    
            }
        } 
    } 
    }

    unsigned int end_time2 =  clock();
    std::cout << end_time2 << "ms" <<std::endl;

    for (auto i : nums_inter) out << i << std::endl;
    out.close();
    fin.close();
}