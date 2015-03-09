//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <cassert>
using namespace std;

int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

struct Sphere
{
    mat4 transform;
    mat4 m_inverse;
    float r, g, b;

    //add lighting later
    Sphere(float x, float y, float z, 
           float sclx, float scly, float sclz,
           float rr, float gg, float bb) : r(rr), g(gg), b(bb) { //add lighting later
        transform = Translate(x,y,z);
        transform *= Scale(sclx, scly, sclz);

        assert(InvertMatrix(transform, m_inverse));
    }
};


// TODO: add structs for spheres, lights and anything else you may need.

vector<vec4> g_colors;
vector<Sphere> g_spheres;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

//ambient light
float amb_r, amb_g, amb_b;
//background color
float back_r, back_g, back_b;
//output filename
char file[21]; //20 character string with 1 null byte
// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;  
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    if (vs[0] == "") return; //ignore newlines
    const string input[11] = {"NEAR", "LEFT", "RIGHT", "BOTTOM", "TOP", "RES","SPHERE", "LIGHT", "BACK", "AMBIENT", "OUTPUT"};
                            //  0       1       2         3       4      5       6      7        8        9          10             
    map<string, int> Mapagrosenhr;
    for (int i = 0; i < 11; i++)
        Mapagrosenhr[input[i]] = i; //map the above strings to int
    try { //in case file has lines with incorrect formatting not in the map
        switch(Mapagrosenhr.at(vs[0])) { 
            case 0: // NEAR
                g_near = toFloat(vs[1]);
                break;
            case 1: // LEFT
                g_left = toFloat(vs[1]);
                break;
            case 2: // RIGHT
                g_right = toFloat(vs[1]);
                break;
            case 3: // BOTTOM
                g_bottom = toFloat(vs[1]);
                break;
            case 4: // TOP
                g_top = toFloat(vs[1]);
                break;
            case 5: // RES: given in template
                g_width  = (int)toFloat(vs[1]);
                g_height = (int)toFloat(vs[2]);
                g_colors.resize(g_width * g_height);
                break;
            case 6: // SPHERE
                //TODO: add a sphere
                assert(g_spheres.size() + 1 <= 5); //assert size limit
                g_spheres.push_back( Sphere(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]), 
                                            toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]),
                                            toFloat(vs[8]), toFloat(vs[9]), toFloat(vs[10])) );
                
                break;
            case 7: // LIGHT
                //TODO: add a light source
                break;
            case 8: // BACK
                back_r = toFloat(vs[1]),
                back_g = toFloat(vs[2]),
                back_b = toFloat(vs[3]);
                break;
            case 9: // AMBIENT
                //TODO: add ambient light
                break;
            case 10: // OUTPUT
                //TODO: add output file
                // cout << vs[1] << vs[1].size()<< endl;
                if (vs[1].size() > 20) {
                    cerr << "Error: filename must be 20 characters or less with no spaces\n";
                    exit(1);
                }
                for(int i = 0; i < vs[1].size(); i++)
                    file[i] = vs[1][i];
                file[vs[1].size()] = '\0';
                // cout << file << endl;
                break;
            default: //if you get here you really broke something
                cout << "How the heck did you get here?\n";
                break;
        }
    }
    catch (const out_of_range& oor){
        // if (vs[0] == "\n" || vs[0] == "\r\n") return;
        cout << "Error " << oor.what() << ": incorrect input file format with " << vs[0] << endl;
        exit(1);
    }
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.


// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.

    return vec4(back_r, back_g, back_b, 1.0f);
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    vec4 dir;
    float x = g_left   + (float)(ix / g_width)  * (g_right - g_left),
          y = g_bottom + (float)(iy / g_height) * (g_top - g_bottom),
          z = -g_near;
    dir = vec4(x,y,z, 0.0f);
    dir = normalize(dir);
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
    
    // TODO: change file name based on input file name.
    savePPM(g_width, g_height, file, buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
    cout << "Spheres: " << g_spheres.size() << endl;
	return 0;
}
