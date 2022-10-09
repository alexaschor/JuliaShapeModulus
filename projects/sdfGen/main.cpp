//SDFGen - A simple grid-based signed distance field (level set) generator for triangle meshes.
//Written by Christopher Batty (christopherbatty@yahoo.com, www.cs.columbia.edu/~batty)
//...primarily using code from Robert Bridson's website (www.cs.ubc.ca/~rbridson)
//This code is public domain. Feel free to mess with it, let me know if you like it.

#include "SETTINGS.h"
#include "makelevelset3.h"
#include "field.h"
#include "projects/sdfGen/vec.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

using namespace std;


int main(int argc, char* argv[]) {

    if(argc < 2) {
        cout << "USAGE: " << endl;
        cout << "To generate an SDF from a mesh with automatically generated bounds:" << endl;
        cout << " " << argv[0] << " <*.obj input> <resolution> <*.f3d output> <padding cells>\n";
        cout << "To get the bounds for a mesh sequence:" << endl;
        cout << " " << argv[0] << " BOUNDS <obj 1> <obj 2> ... <obj N>\n";
        cout << "To generate an SDF from a mesh with specified bounds:" << endl;
        cout << " " << argv[0] << " <*.obj input> <resolution> <*.f3d output> <min X> <min Y> <min Z> <max X> <max Y> <max Z>\n";
        exit(-1);
    }

    if (string(argv[1]) == "BOUNDS") {
        //start with a massive inside out bound box.
        Vec3f min_box(numeric_limits<float>::max(),numeric_limits<float>::max(),numeric_limits<float>::max()),
              max_box(-numeric_limits<float>::max(),-numeric_limits<float>::max(),-numeric_limits<float>::max());

        for (int i = 2; i < argc; ++i) {
            string filename(argv[i]);

            if(filename.size() < 5 || filename.substr(filename.size()-4) != string(".obj")) {
                cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
                exit(-1);
            }

            cout << "Reading data.\n";

            ifstream infile(argv[i]);
            if(!infile) {
                cerr << "Failed to open " << argv[i] << " Terminating.\n";
                exit(-1);
            }

            int ignored_lines = 0;
            int total_lines = 0;
            string line;
            while(!infile.eof()) {
                getline(infile, line);

                //.obj files sometimes contain vertex normals indicated by "vn"
                if(line.substr(0,1) == string("v") && line.substr(0,2) != string("vn")){
                    stringstream data(line);
                    char c;
                    Vec3f point;
                    data >> c >> point[0] >> point[1] >> point[2];
                    update_minmax(point, min_box, max_box);
                }
                else if( line.substr(0,2) == string("vn") ){
                    cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
                    exit(-2);
                }
                else if(line.substr(0,1) == string("f")) { }
                else {
                    ++ignored_lines;
                }
                ++total_lines;
            }
            infile.close();

            if(ignored_lines > 0)
                cout << "Warning: " << ignored_lines << " of " << total_lines << " lines were ignored since they did not contain faces or vertices.\n";

        }

        PRINTV3(min_box);
        PRINTV3(max_box);

        exit(0);
    }

    string filename(argv[1]);
    if(filename.size() < 5 || filename.substr(filename.size()-4) != string(".obj")) {
        cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
        exit(-1);
    }

    int res = atoi(argv[2]);

    int padding = 1;

    if (argc == 5) {
        padding = atoi(argv[4]);
    }

    //start with a massive inside out bound box.
    Vec3f min_box(numeric_limits<float>::max(),numeric_limits<float>::max(),numeric_limits<float>::max()),
          max_box(-numeric_limits<float>::max(),-numeric_limits<float>::max(),-numeric_limits<float>::max());

    cout << "Reading data.\n";

    ifstream infile(argv[1]);
    if(!infile) {
        cerr << "Failed to open. Terminating.\n";
        exit(-1);
    }

    int ignored_lines = 0;
    string line;
    vector<Vec3f> vertList;
    vector<Vec3ui> faceList;
    while(!infile.eof()) {
        getline(infile, line);

        //.obj files sometimes contain vertex normals indicated by "vn"
        if(line.substr(0,1) == string("v") && line.substr(0,2) != string("vn")){
            stringstream data(line);
            char c;
            Vec3f point;
            data >> c >> point[0] >> point[1] >> point[2];
            vertList.push_back(point);
            update_minmax(point, min_box, max_box);
        }
        else if(line.substr(0,1) == string("f")) {
            stringstream data(line);
            char c;
            int v0,v1,v2;
            data >> c >> v0 >> v1 >> v2;
            faceList.push_back(Vec3ui(v0-1,v1-1,v2-1));
        }
        else if( line.substr(0,2) == string("vn") ){
            cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
            exit(-2);
        }
        else {
            ++ignored_lines;
        }
    }
    infile.close();

    if(ignored_lines > 0)
        cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";

    cout << "Read in " << vertList.size() << " vertices and " << faceList.size() << " faces." << endl;


    PRINT("Bounds from mesh:");
    PRINTV3(min_box);
    PRINTV3(max_box);

    if (argc >= 10) {
        float minX = atof(argv[4]);
        float minY = atof(argv[5]);
        float minZ = atof(argv[6]);

        float maxX = atof(argv[7]);
        float maxY = atof(argv[8]);
        float maxZ = atof(argv[9]);

        if ((minX > min_box[0]) || (minY > min_box[1]) || (minZ > min_box[2]) ||
            (maxX < max_box[0]) || (maxY < max_box[1]) || (maxZ < max_box[1]))
        {
            PRINT("WARNING: User-specified bounds will truncate mesh!");
        }

        min_box = Vec3f(minX, minY, minZ);
        max_box = Vec3f(maxX, maxY, maxZ);

        PRINT("NOTE: Overriding with user-specified bounds");
        PRINTV3(min_box);
        PRINTV3(max_box);
    }


    // find the center
    Vec3f center = (min_box + max_box) * 0.5;
    Vec3f lengths = max_box - min_box;
    double maxLength = lengths[0] > lengths[1] ? lengths[0] : lengths[1];
    maxLength = maxLength > lengths[2] ? maxLength : lengths[2];
    lengths = Vec3f(maxLength, maxLength, maxLength);

    min_box = center - lengths * 0.5;
    max_box = center + lengths * 0.5;

    float dx = maxLength / res;
    cout << " dx: " << dx << endl;

    //Add padding around the box.
    Vec3f unit(1,1,1);
    min_box -= padding*dx*unit;
    max_box += padding*dx*unit;

    // now re-do dx to include the padding as well
    lengths = max_box - min_box;
    dx = lengths[0] / res;
    Vec3ui sizes = Vec3ui((max_box - min_box)/dx);

    cout << "Computing signed distance field.\n";
    SDFArray3F phi_grid;
    make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], phi_grid);

    string outname(argv[3]);

    //Very hackily strip off file suffix.
    cout << "Writing results to: " << outname << "\n";

    // create a FIELD_3D object to save
    int xRes = phi_grid.ni;
    int yRes = phi_grid.nj;
    int zRes = phi_grid.nk;
    Vec3f myLengths(lengths[0], lengths[1], lengths[2]);
    Vec3f myCenter(center[0], center[1], center[2]);
    cout << " Writing field3d file with center: " << myCenter << " lengths: " << myLengths << endl;

    // FIXME ugly conversion, really should standardize vector
    // types

    VEC3F minV3F(min_box[0], min_box[1], min_box[2]);
    VEC3F maxV3F(max_box[0], max_box[1], max_box[2]);

    ArrayGrid3D field(xRes, yRes, zRes);
    field.setMapBox(AABB(minV3F, maxV3F));

    for (int z = 0; z < zRes; z++)
        for (int y = 0; y < yRes; y++)
            for (int x = 0; x < xRes; x++)
                field.at(x,y,z) = phi_grid(x,y,z);

    field.writeF3D(outname, true);

    cout << "Processing complete.\n";

    return 0;
}
