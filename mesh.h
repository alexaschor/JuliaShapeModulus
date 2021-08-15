#ifndef MESH_H
#define MESH_H
#include <fstream>
#include <string>

#include "SETTINGS.h"

class Mesh {
public:
    std::vector<VEC3F> vertices;
    std::vector<VEC3F> normals;
    std::vector<uint> indices;

    void writeOBJ(std::string filename) {
        std::ofstream out;
        out.open(filename);
        if (out.is_open() == false)
            return;
        out << "g " << "Obj" << std::endl;
        for (size_t i = 0; i < vertices.size(); i++)
            out << "v " << vertices.at(i).x() << " " << vertices.at(i).y() << " " << vertices.at(i).z() << '\n';
        for (size_t i = 0; i < vertices.size(); i++)
            out << "vn " << normals.at(i).x() << " " << normals.at(i).y() << " " << normals.at(i).z() << '\n';
        for (size_t i = 0; i < indices.size(); i += 3)
        {
            out << "f " << indices.at(i) + 1 << "//" << indices.at(i) + 1
                << " " << indices.at(i + 1) + 1 << "//" << indices.at(i + 1) + 1
                << " " << indices.at(i + 2) + 1 << "//" << indices.at(i + 2) + 1
                << '\n';
        }
        out.close();

    }
};

#endif
