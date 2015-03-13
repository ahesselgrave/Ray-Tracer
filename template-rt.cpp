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
    float Ka, Kd, Ks, Kr, n;

    //add lighting later
    Sphere(float x, float y, float z, 
           float sclx, float scly, float sclz,
           float rr, float gg, float bb,
           float a, float d, float s,
           float r, float nn) : r(rr), g(gg), b(bb), Ka(a), Kd(d), Ks(s), Kr(r), n(nn) {

        transform  = Translate(x,y,z);
        transform *= Scale(sclx, scly, sclz);

        assert(InvertMatrix(transform, m_inverse));
    }
};

struct Light
{
    float Ir, Ig, Ib; //intensities
    vec4 pos;

    Light(float x, float y, float z,
          float Ir, float Ig, float Ib) : Ir(Ir), Ig(Ig), Ib(Ib) {
        pos = vec4(x, y, z, 1.0f);
    }
};

vector<vec4> g_colors;
vector<Sphere> g_spheres;
vector<Light> g_lights;

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
//zero vector
const vec4 zeroVec = vec4(0.0f, 0.0f, 0.0f, 0.0f);
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
                                            toFloat(vs[8]), toFloat(vs[9]), toFloat(vs[10]),
                                            toFloat(vs[11]), toFloat(vs[12]), toFloat(vs[13]),
                                            toFloat(vs[14]), toFloat(vs[15])) );
                break;
            case 7: // LIGHT
                assert(g_lights.size() + 1 <= 5); //assert size limit
                g_lights.push_back( Light(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]),
                                          toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7])) );
                break;
            case 8: // BACK
                back_r = toFloat(vs[1]),
                back_g = toFloat(vs[2]),
                back_b = toFloat(vs[3]);
                break;
            case 9: // AMBIENT
                amb_r = toFloat(vs[1]),
                amb_g = toFloat(vs[2]),
                amb_b = toFloat(vs[3]);
                break;
            case 10: // OUTPUT
                if (vs[1].size() > 20) {
                    cerr << "Error: filename must be 20 characters or less with no spaces\n";
                    exit(1);
                }
                for(int i = 0; i < vs[1].size(); i++)
                    file[i] = vs[1][i];
                file[vs[1].size()] = '\0';
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

float posQuad(float A, float B, float C) {
    return -B/A + std::sqrt(B*B-A*C) / A;
}

float negQuad(float A, float B, float C) {
    return -B/A - std::sqrt(B*B-A*C) / A;
}


// TODO: add your ray-sphere intersection routine here.
vec4 lightContribution(const Ray& ray, const vec4& N, 
                       const vec4& hitpoint, const int& sphereNum);

vec4 rayIntersectsSphere(const Ray& ray, int recursionLevel) {
    vec4 S = ray.origin, c = ray.dir;
    vector<float> t_hVec; //store all t_h values
    for(int i = 0; i < g_spheres.size(); i++) {
        //initialize vec3 untransformed spheres
        vec4 S_prime4 = g_spheres[i].m_inverse * S;
        vec3 S_prime3(S_prime4.x, S_prime4.y, S_prime4.z);

        vec4 c_prime4 = g_spheres[i].m_inverse * c;
        vec3 c_prime3(c_prime4.x, c_prime4.y, c_prime4.z);


        float A = dot(c_prime3, c_prime3);
        float B = dot(S_prime3, c_prime3);
        float C = dot(S_prime3, S_prime3)- 1;

        float discriminant = B*B - A*C;
        if (discriminant < 0) continue;
        else {
            float pos_t = posQuad(A,B,C);
            float neg_t = negQuad(A,B,C);

            float t_h = pos_t < neg_t ? pos_t : neg_t; //assign to lowest t
            if (t_h <= 1) t_hVec.push_back(0); //we will ignore all 0 values
            else t_hVec.push_back(t_h);
            {
                //color starts out as sphere color * ambient
                vec4 color = vec4(g_spheres[i].r * g_spheres[i].Ka, 
                                  g_spheres[i].g * g_spheres[i].Ka, 
                                  g_spheres[i].b * g_spheres[i].Ka, 1.0f);
                //add color from light sources
                //first get normal vector
                vec4 hitpoint = S + t_h * c;
                vec4 unitHitpoint = S_prime4 + t_h * c_prime4;
                vec4 normal = vec4(unitHitpoint.x, unitHitpoint.y, unitHitpoint.z, 0.0f);
                vec4 trans_normal = transpose(g_spheres[i].m_inverse) * normal;
                trans_normal = normalize(trans_normal);

                color += lightContribution(ray, trans_normal, hitpoint, i);
                return color;
            }
        }
    }
    //if it intersects with no spheres, return background color
    return vec4(back_r, back_g, back_b, 1.0f);
}

vec4 lightContribution(const Ray& ray, const vec4& N, const vec4& hitpoint, const int& sphereNum) {
    //assume N is normalized
    vec4 color = zeroVec; //initialize color contribution as empty
    for (int i = 0; i < g_lights.size(); i++) {
        Light light = g_lights[i]; //easier access

        vec4 L = light.pos - hitpoint;
        L = normalize(L);

        vec4 R = ((2 * N) * dot(N,L)) - L;
        R = normalize(R);

        //diffuse
        float cosDifAngle = dot(N,L);
        // if (cosDifAngle )
        vec4 diffuse = vec4(light.Ir*g_spheres[sphereNum].Kd*cosDifAngle,
                            light.Ig*g_spheres[sphereNum].Kd*cosDifAngle,
                            light.Ib*g_spheres[sphereNum].Kd*cosDifAngle, 1.0f);
        color += diffuse;

        //specular
        vec4 v = ray.dir;
        float cosnSpecAngle = dot(R,v);
        cosnSpecAngle = pow(cosnSpecAngle, g_spheres[i].n);
        vec4 specular = vec4(light.Ir*g_spheres[sphereNum].Ks*cosnSpecAngle,
                             light.Ig*g_spheres[sphereNum].Ks*cosnSpecAngle,
                             light.Ib*g_spheres[sphereNum].Ks*cosnSpecAngle, 1.0f);
        color += specular;
    }
    //ambient
    vec4 ambient = vec4(amb_r, amb_g, amb_b, 1.0f);
    ambient *= g_spheres[sphereNum].Ka;
    color += ambient;
    color.w = 1.0f;
    return color;
}

// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
    return rayIntersectsSphere(ray, 1);
    return vec4(back_r, back_g, back_b, 1.0f);
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    vec4 dir;
    float x = g_left   + ((float)ix / (float)g_width)  * (g_right - g_left),
          y = g_bottom + ((float)iy / (float)g_height) * (g_top - g_bottom),
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
    for (int i = 0; i < g_colors.size(); i++){
        for (int j = 0; j < 3; j++){
            if (g_colors[i][j] < 0.0f) g_colors[i][j] = 0.0f;
            else if (g_colors[i][j] > 1.0f) g_colors[i][j] = 1.0f;
            else continue;
        }
    }
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
    
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
	return 0;
}
