#ifndef MESH_H
#define MESH_H
#include <fstream>
#include <iostream>
#include <cstdio>
#include <string>

#include "SETTINGS.h"
#include "triangle.h"

using namespace std;

class Mesh {
public:
    std::vector<VEC3F> vertices;
    std::vector<VEC3F> normals;
    std::vector<uint> indices;

    int numFaces() {
        assert(indices.size() % 3 == 0);
        return indices.size() / 3;
    }

    //readOBJ function adapted from http://www.opengl-tutorial.org
    void readOBJ(std::string filename) {
        FILE * file = fopen(filename.c_str(), "r");
        if (file == NULL){
            printf("Could not open OBJ file %s for reading.\n", filename.c_str());
            exit(1);
        }
        while (true){
            char lineHeader[128];
            // read the first word of the line
            int res = fscanf(file, "%s", lineHeader);
            if (res == EOF)
                break;

            if (strcmp(lineHeader, "v") == 0) {
                VEC3F vertex;
                fscanf(file, "%lf %lf %lf\n", &vertex.x(), &vertex.y(), &vertex.z());
                vertices.push_back(vertex);
            } else if (strcmp(lineHeader, "f") == 0){
                size_t a_i, b_i, c_i;
                int matches = fscanf(file, "%zd %zd %zd\n", &a_i, &b_i, &c_i);
                if (matches != 3) {
                    printf("Encountered malformed face data when reading OBJ %s. Make sure that UV and vertex normals are not included in the file.\n", filename.c_str());
                    exit(1);
                }

                indices.push_back(a_i-1);
                indices.push_back(b_i-1);
                indices.push_back(c_i-1);
            }
        }

        printf("Read %lu vertices and %lu faces from %s\n", vertices.size(), indices.size() / 3, filename.c_str());
    }

    void writeOBJ(std::string filename) {
        std::ofstream out;
        out.open(filename);
        if (out.is_open() == false)
            return;
        out << "g " << "Obj" << std::endl;
        for (size_t i = 0; i < vertices.size(); i++)
            out << "v " << vertices.at(i).x() << " " << vertices.at(i).y() << " " << vertices.at(i).z() << '\n';
        for (size_t i = 0; i < normals.size(); i++)
            out << "vn " << normals.at(i).x() << " " << normals.at(i).y() << " " << normals.at(i).z() << '\n';
        for (size_t i = 0; i < indices.size(); i += 3) {
            out << "f " << indices.at(i) + 1 << "//" << indices.at(i) + 1
                << " " << indices.at(i + 1) + 1 << "//" << indices.at(i + 1) + 1
                << " " << indices.at(i + 2) + 1 << "//" << indices.at(i + 2) + 1
                << '\n';
        }
        out.close();

        std::cout << "Wrote " << vertices.size() << " vertices and " << indices.size() / 3 << " faces to " << filename << std::endl;
    }

    Triangle triangle(int idx) {
        idx *= 3;
        Triangle out(&(vertices[indices[idx]]), &(vertices[indices[idx] + 1]), &(vertices[indices[idx] + 2]));
        return out;
    }

};

#endif
