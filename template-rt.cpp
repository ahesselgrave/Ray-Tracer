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
#include <algorithm>
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
const vec4 black = vec4(0.0f, 0.0f, 0.0f, 1.0f);
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
                exit(1);
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
inline float posQuad(float A, float B, float D) {
    return -B/A + sqrt(D)/A;
}

inline float negQuad(float A, float B, float D){
    return -B/A - sqrt(D)/A;
}

vec4 lightContribution(const vec4& P, const vec4& nn, int k, const Ray& ray) {
    //find L, r, v
    vec4 baseColor = vec4(g_spheres[k].r, g_spheres[k].g, g_spheres[k].b, 1.0f);
    float Ka = g_spheres[k].Ka;   
    float Kd = g_spheres[k].Kd;
    float Ks = g_spheres[k].Ks;
    vec4 color = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    vec3 n = vec3(nn.x, nn.y, nn.z);
    n = normalize(n);

    for (int i = 0; i < g_lights.size(); i++){
        //check for shadows first
        vec4 S = P;
        vec4 c = g_lights[i].pos - P;
        bool shadowExists = false;
        for (int j = 0; j < g_spheres.size(); j++){ //for every sphere
            vec4 S_prime4 = g_spheres[j].m_inverse * S;
            vec3 S_prime3 = vec3(S_prime4.x, S_prime4.y, S_prime4.z);

            vec4 c_prime4 = g_spheres[j].m_inverse * c;
            vec3 c_prime3 = vec3(c_prime4.x, c_prime4.y, c_prime4.z);

            //find th
            float A = dot(c_prime3, c_prime3),     //|c|^2
                  B = dot(S_prime3, c_prime3),     //S . c
                  C = dot(S_prime3, S_prime3) - 1; //|S|^2 - 1
            float D = B*B - A*C; //discriminant
            if (D < 0) continue; //move along, no intersection here
            float th1 = posQuad(A, B, D),
                  th2 = negQuad(A, B, D);
            if ((th1 > 0.0001 && th1 < 1) || (th2 > 0.001 && th2 < 1)) {
                shadowExists = true;
                break;
            }
        }
        if (!shadowExists) {
            vec4 LL = g_lights[i].pos - P;
            vec3 L = vec3(LL.x, LL.y, LL.z);
            L = normalize(L);
    
            vec3 r = (2*dot(n,L) * n) - L;
            r = normalize(r);
         
            vec4 vv = ray.origin - P;
            vec3 v = vec3(vv.x, vv.y, vv.z);
            v = normalize(v);
    
            float Ir = g_lights[i].Ir,
                  Ig = g_lights[i].Ig,
                  Ib = g_lights[i].Ib;
    
            //diffuse
            float cosDifAngle = dot(n,L);
            vec4 diffuse = black, specular = black;
            if (cosDifAngle < 0)
                diffuse = vec4(0.0f, 0.0f, 0.0f, 1.0f);
            else {
                diffuse = vec4(Kd * Ir * cosDifAngle * baseColor.x,
                               Kd * Ig * cosDifAngle * baseColor.y,
                               Kd * Ib * cosDifAngle * baseColor.z, 1.0f);
                //specular
                if (dot(r,v) < 0)
                    specular = vec4(0.0f, 0.0f, 0.0f, 1.0f);
                else {
                    double cosnSpecAngle = pow(dot(r,v), g_spheres[k].n);
                    specular = vec4(Ks * Ir * cosnSpecAngle,
                                    Ks * Ig * cosnSpecAngle,
                                    Ks * Ib * cosnSpecAngle, 1.0f);
                }
            }   
            color += (diffuse + specular);
        }
    }
    return color;
}

// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray, int depth)
{
    if (depth > 3) return black; //base case

    vec4 S = ray.origin, c = ray.dir;
    map<float, int> sphereMap;
    for (int i = 0; i < g_spheres.size(); i++){ //for every sphere
        vec4 S_prime4 = g_spheres[i].m_inverse * S;
        vec3 S_prime3 = vec3(S_prime4.x, S_prime4.y, S_prime4.z);

        vec4 c_prime4 = g_spheres[i].m_inverse * c;
        vec3 c_prime3 = vec3(c_prime4.x, c_prime4.y, c_prime4.z);

        //find th
        float A = dot(c_prime3, c_prime3),     //|c|^2
              B = dot(S_prime3, c_prime3),     //S . c
              C = dot(S_prime3, S_prime3) - 1; //|S|^2 - 1
        float D = B*B - A*C; //discriminant
        if (D < 0) continue; //move along, no intersection here
        float th1 = posQuad(A, B, D),
              th2 = negQuad(A, B, D);
        float th  = th1 < th2 ? th1 : th2;
        if (th > 1)
            sphereMap[th] = i;
    }
    //maps are default sorted, so first element should be first sphere
    if (sphereMap.size() == 0) { //no intersection
        if(depth == 0)
            return vec4(back_r, back_g, back_b, 1.0f);
        else 
        {
            cerr << "Depth " <<depth<<" and no hits!" << endl;
            return black;
        }
    }
    map<float, int>::iterator iter = sphereMap.begin();
    float th = iter->first;
    int i = iter->second;

    //insert lighting here
    vec4 P = S + (c * th); //sphere hitpoint
    vec4 S_prime4 = g_spheres[i].m_inverse * S;
    vec4 c_prime4 = g_spheres[i].m_inverse * c;

    vec4 N = S_prime4 + (c_prime4 * th); //unit sphere hitpoint
    N.w = 0.0f; //now normal vec
    N = transpose(g_spheres[i].m_inverse)*N; //M^-t * N 
    N = normalize(N);

    //ambient
    vec4 ambient(amb_r * g_spheres[i].Ka * g_spheres[i].r,
                 amb_g * g_spheres[i].Ka * g_spheres[i].g,
                 amb_b * g_spheres[i].Ka * g_spheres[i].b, 1.0f);

    vec4 color = ambient;
    color += lightContribution(P, N, i, ray);
    //reflection
    vec4 v = -2*dot(N,c)*N + c;
    v = normalize(v);

    Ray reflection;
    reflection.origin = P;
    reflection.dir = v;

    vec4 reflectColor = g_spheres[i].Kr * trace(reflection, depth + 1);
    color += reflectColor;
    color.w = 1.0f;
    return color;
}

vec4 getDir(int ix, int iy)
{
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
    vec4 color = trace(ray, 0);
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
