#include <iostream>
#include <fstream>
#include <sstream>
#include<bits/stdc++.h>
#include <math.h>
#include <limits.h>
#include <vector>
#include <algorithm>
#include "glm/glm.hpp"

class FovType {
    public:
    float h;
    float v;
};

class ImsizeType{
    public:
    int width;
    int height;
};

class ViewingwindowType{
    public:
    float height;
    float width;
    glm::vec3 ul;
    glm::vec3 ur;
    glm::vec3 ll;
    glm::vec3 lr;
};

class ColorType{
    public:
    float r;
    float g;
    float b;
};

class TextureType{
    public:
    int width;
    int height;
    glm::vec3 **texturemap;
};

class MaterialType{
    public:
    glm::vec3 od;
    glm::vec3 os;
    float ka;
    float kd;
    float ks;
    int n;
    bool justrgb = true;
};

class RayType{
    public:
    glm::vec3 position;
    glm::vec3 direction;
};

class SphereType {
    public:
    float radius;
    glm::vec3 position;
    MaterialType m;
    bool justrgb;
    TextureType texture;
    bool textureStatus = false;
};

class LightType{
    public:
    glm::vec3 position;
    int type;
    glm::vec3 color;
};

class FaceType{
    public:
    bool SS = false;
    glm::vec3 position;
    MaterialType m;
    int index1;
    int index2;
    int index3;
    TextureType texture;
    bool textureStatus = false;

};

class TriangleType{
    public:
    glm::vec3 p0;
    glm::vec3 p1;
    glm::vec3 p2;
    glm::vec3 e1;
    glm::vec3 e2;
    float alpha;
    float beta;
    float gamma;
};

class PlaneType{
    public:
    bool SS;
    float A;
    float B;
    float C;
    float D;
    MaterialType m;
    TriangleType triangle;
    glm::vec3 n0;
    glm::vec3 n1;
    glm::vec3 n2;
    glm::vec2 vt0;
    glm::vec2 vt1;
    glm::vec2 vt2;
    TextureType texture;
    bool textureStatus = false;
};

std::vector<LightType> lights;  //global variables
std::vector<SphereType> objects;
std::vector<glm::vec3> points;
std::vector<FaceType> faces;
std::vector<PlaneType> planes;
std::vector<glm::vec3> shading_Points;
std::vector<TextureType> textures;
std::vector<glm::vec2> texture_Points;

FaceType Plane(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, PlaneType planeToPush){
    glm::vec3 e1 = p1 - p0;
    glm::vec3 e2 = p2 - p0;
    glm::vec3 n = glm::cross(e1, e2);
    glm::normalize(n);
    FaceType face;
    face.position = n; 
    return face;
}

void Populate_planes(std::vector<glm::vec3> points, std::vector<FaceType> faces){
    for(int i = 0; i < faces.size(); i++){ 
        glm::vec3 p0 = points.at(faces.at(i).position.x);
        glm::vec3 p1 = points.at(faces.at(i).position.y);
        glm::vec3 p2 = points.at(faces.at(i).position.z);
        PlaneType planeToPush;
        planeToPush.triangle.p0 = p0;
        planeToPush.triangle.p1 = p1;
        planeToPush.triangle.p2 = p2;
        planeToPush.triangle.e1 = p1 - p0;
        planeToPush.triangle.e2 = p2 - p0;
        FaceType temp = Plane(p0, p1, p2, planeToPush);
        planeToPush.A = temp.position.x;
        planeToPush.B = temp.position.y;
        planeToPush.C = temp.position.z;
        planeToPush.D = -1 * (planeToPush.A * p0.x + planeToPush.B * p0.y + planeToPush.C * p0.z);
        planeToPush.m = faces.at(i).m;
        if(faces.at(i).SS == true){
            planeToPush.SS = true;
            planeToPush.n0 = shading_Points.at(faces.at(i).index1);
            planeToPush.n1 = shading_Points.at(faces.at(i).index2);
            planeToPush.n2 = shading_Points.at(faces.at(i).index3);
        }
        else{
            planeToPush.SS = false;
        }
        planeToPush.textureStatus = faces.at(i).textureStatus;
        if(planeToPush.textureStatus == true){
            planeToPush.texture = faces.at(i).texture;
            planeToPush.vt0 = texture_Points.at(faces.at(i).index1);
            planeToPush.vt1 = texture_Points.at(faces.at(i).index2);
            planeToPush.vt2 = texture_Points.at(faces.at(i).index3);
        }
        planes.push_back(planeToPush);
    }
}

int Shadow_status(glm::vec3 intersectionPoint, glm::vec3 L, int o_num, float light_distance, bool shapetype){
    if(shapetype == 1){
        o_num = -1;
    }
    int current_obj = -1; 
    int iter = -1;
    for(auto sphere : objects){
        iter += 1; 
        if(iter != o_num){
            float xc = intersectionPoint.x - sphere.position.x;
            float yc = intersectionPoint.y - sphere.position.y;
            float zc = intersectionPoint.z - sphere.position.z;
            float A = 1.0;
            float B = 2 * (L.x * (xc) + L.y * (yc) + L.z * (zc));
            float C = (xc) * (xc) + (yc) * (yc) + (zc) * (zc) - sphere.radius * sphere.radius;
            float disc = (B * B) - (4 * C);
            if(disc > 0){
                float t1 = (-B + sqrt(disc)) / (2 * A);
                float t2 = (-B - sqrt(disc)) / (2 * A);
                if(t1 > 0 && t1 < light_distance){
                    return 0;
                }
                if(t2 > 0 && t2 < light_distance){
                    return 0;
                }
            }
            else if(disc == 0){
                float t1 = (-B / 2 * A);
                if(t1 > 0 && t1 < light_distance){
                    return 0;
                }
            }
        }
    }
    return 1;
}

int Shadow_Status_triangle(glm::vec3 intersectionPoint, glm::vec3 L, int o_num, float light_distance, bool shapetype){
    if(o_num < planes.size()){
       PlaneType plane = planes.at(o_num); 
    }
    else return 1;
    if(shapetype == 0){
        o_num = -1;
    }
    for(int i = 0; i < planes.size(); i++){
        if(i != o_num){
            float denom = planes.at(i).A * L.x + planes.at(i).B *  L.y + planes.at(i).C * L.z;
            if(denom == 0){
                std::cout << "denom is zero\n";
                return 0;
            }
            float numer = -1.0 * (planes.at(i).A * L.x + planes.at(i).B * L.y + planes.at(i).C * L.z + planes.at(i).D);
            float t3 = numer/denom;
            if(t3 > 0 && t3 < light_distance){
                return 0;
            }
        }
    } 
    return 1;
}

int Full_Shadow_Status(glm::vec3 intersectionPoint, glm::vec3 L, int o_num, float light_distance, bool shapetype){ //0 for sphere, 1 for triangle;
    int result0 = Shadow_Status_triangle(intersectionPoint, L, o_num, light_distance, shapetype);
    int result1 = Shadow_status(intersectionPoint, L, o_num, light_distance, shapetype);
    if(result0 == 0 || result1 == 0){
        return 0;
    }
    else return 1;
}

ColorType Shade_ray(SphereType sphere){
    ColorType colorToReturn;
    colorToReturn.r = sphere.m.od.x;
    colorToReturn.g = sphere.m.od.y;
    colorToReturn.b = sphere.m.od.z;
    return colorToReturn;
}

ColorType Shade_ray2(SphereType sphere, glm::vec3 intersectionPoint, RayType ray, int o_num){   //for use when there is lighting in the scene
    ColorType colorToReturn;
    glm::vec3 vectorColor;
    glm::vec3 N;
    glm::vec3 H;
    glm::vec3 L;
    glm::vec3 sum;
    sum.x = 0;
    sum.y = 0;
    sum.z = 0;
    glm::vec3 newVec;
    float light_distance;
    glm::vec3 v = glm::normalize(intersectionPoint - ray.position);
    N = (intersectionPoint - sphere.position)/ sphere.radius;
    if(sphere.textureStatus == true){
        float phi_radians = acos((intersectionPoint.z - sphere.position.z)/sphere.radius);
        float theta_radians = atan2((intersectionPoint.y - sphere.position.y),(intersectionPoint.x - sphere.position.x));
        float v = phi_radians / M_PI;
        float u;
        if(theta_radians < 0){
            u = (theta_radians + M_PI)/(2*M_PI);
        }
        else{
        u = theta_radians /(2* M_PI);   
        }
        sphere.m.od = (sphere.texture.texturemap[(int)(v * (sphere.texture.height -1))][(int)(u * (sphere.texture.width -1))])/(float)(255.0);
    }
    for(int i = 0; i < lights.size(); i++){
        if(lights.at(i).type == 0){ //directional
            L = glm::normalize(float(-1) * lights.at(i).position);
            light_distance = FLT_MAX; //acts as infinity 
        }
        else if(lights.at(i).type == 1){
            L = glm::normalize(lights.at(i).position - intersectionPoint);
            light_distance = sqrt(pow(lights.at(i).position.x - intersectionPoint.x, 2) + pow(lights.at(i).position.y - intersectionPoint.y, 2) + pow(lights.at(i).position.z - intersectionPoint.z, 2));
        }
        H = glm::normalize(L + v);
        newVec = lights.at(i).color * ((sphere.m.kd) * sphere.m.od * (std::max(float(0), glm::dot(N, L))) + sphere.m.ks * sphere.m.os * float(pow(std::max(float(0), glm::dot(N, H)), sphere.m.n))); //calculate value in sum
        sum += newVec * float(Full_Shadow_Status(intersectionPoint, L, o_num, light_distance, false));
    }
    vectorColor.x = sphere.m.ka * sphere.m.od.x;
    vectorColor.y = sphere.m.ka * sphere.m.od.y;
    vectorColor.z = sphere.m.ka * sphere.m.od.z;
    vectorColor += sum;
    colorToReturn.r = std::min(vectorColor.x, 1.0f);
    colorToReturn.g = std::min(vectorColor.y, 1.0f);
    colorToReturn.b = std::min(vectorColor.z, 1.0f);
    return colorToReturn;
}

ColorType Shade_Ray_Triangle2(PlaneType plane, glm::vec3 intersectionPoint, RayType ray, int o_num, float min_d){
    ColorType colorToReturn;
    glm::vec3 vectorColor;
    glm::vec3 N;
    glm::vec3 H;
    glm::vec3 L;
    glm::vec3 sum;
    sum.x = 0;
    sum.y = 0;
    sum.z = 0;
    glm::vec3 newVec;
    float light_distance;
    glm::vec3 v = glm::normalize(ray.position - intersectionPoint);
    if(plane.SS == true){
        N = glm::normalize(float(plane.triangle.alpha) * plane.n0 + float(plane.triangle.beta) * plane.n1 + float(plane.triangle.gamma) * plane.n2);
    }
    else{
        N = glm::normalize(glm::cross(plane.triangle.e1, plane.triangle.e2));
    }
    if(plane.textureStatus == true){
        float u = (plane.triangle.alpha * plane.vt0.x) + (plane.triangle.beta * plane.vt1.x) + (plane.triangle.gamma * plane.vt2.x);
        float v = (plane.triangle.alpha * plane.vt0.y) + (plane.triangle.beta * plane.vt1.y) + (plane.triangle.gamma * plane.vt2.y);
        std::ofstream newfile;
        // newfile.open("created.ppm");
        // for(int i = 0; i < plane.texture.height - 1; i++){
        //     for(int j = 0; j < plane.texture.width - 1; j++){
        //         newfile << plane.texture.texturemap[i][j].x << " " << plane.texture.texturemap[i][j].y << " " << plane.texture.texturemap[i][j].z << std::endl;
        //     }
        // }
        // newfile.close();
        int v2 = (int)(v * (plane.texture.height -1));
        int u2 = (int)(u * (plane.texture.width -1));
        std::cout << u2 << " " << v2 << std::endl;
        plane.m.od = (plane.texture.texturemap[(int)(v * (plane.texture.height -1))][(int)(u * (plane.texture.width -1))])/(float)(255.0);
    }
    for(int i = 0; i < lights.size(); i++){
        if(lights.at(i).type == 0){ //directional
            L = glm::normalize(float(-1) * lights.at(i).position);
            light_distance = FLT_MAX; //acts as infinity
        }
        else if(lights.at(i).type == 1){
            L = glm::normalize(lights.at(i).position - intersectionPoint);
            light_distance = sqrt(pow(lights.at(i).position.x - intersectionPoint.x, 2) + pow(lights.at(i).position.y - intersectionPoint.y, 2) + pow(lights.at(i).position.z - intersectionPoint.z, 2));
        }
        H = glm::normalize(L + v);
        newVec = lights.at(i).color * ((plane.m.kd) * plane.m.od * (std::max(float(0), glm::dot(N, L))) + plane.m.ks * plane.m.os * float(pow(std::max(float(0), glm::dot(N, H)), plane.m.n))); //calculate value in sum
        sum += newVec * float(Full_Shadow_Status(intersectionPoint, L, o_num, light_distance, true)); //intersectionpt + N
    }

    vectorColor.x = plane.m.ka * plane.m.od.x;
    vectorColor.y = plane.m.ka * plane.m.od.y;
    vectorColor.z = plane.m.ka * plane.m.od.z;
    vectorColor += sum;
    colorToReturn.r = std::min(vectorColor.x, 1.0f);
    colorToReturn.g = std::min(vectorColor.y, 1.0f);
    colorToReturn.b = std::min(vectorColor.z, 1.0f);
    return colorToReturn;
}

ColorType Trace_Ray(RayType ray, std::vector<SphereType> objects, ColorType background){
    float min_d = FLT_MAX;
    int current_obj = -1; 
    int iter = -1;
    glm::vec3 intersectionPoint;
    int object_type;
    for(auto sphere : objects){
        iter += 1; 
        float xc = ray.position.x - sphere.position.x;
        float yc = ray.position.y - sphere.position.y;
        float zc = ray.position.z - sphere.position.z;
        float A = 1.0;
        float B = 2 * (ray.direction.x * (xc) + ray.direction.y * (yc) + ray.direction.z * (zc));
        float C = (xc) * (xc) + (yc) * (yc) + (zc) * (zc) - sphere.radius * sphere.radius;
        float disc = (B * B) - (4 * C);
        if(disc > 0){
            float t1 = (-B + sqrt(disc)) / (2 * A);
            float t2 = (-B - sqrt(disc)) / (2 * A);
            if(t1 > 0 && t1 < min_d){
                min_d = t1;
                current_obj = iter;
                object_type = 1; //sphere
            }
            if(t2 > 0 && t2 < min_d){
                min_d = t2;
                current_obj = iter;
                object_type = 1; //sphere
            }
        }
        else if(disc == 0){
            float t1 = (-B / 2 * A);
            if(t1 > 0 && t1 < min_d){
                min_d = t1;
                current_obj = iter;
                object_type = 1; //sphere
            }
        }
    }
    for(int i = 0; i < planes.size(); i++){
        float denom = planes.at(i).A * ray.direction.x + planes.at(i).B *  ray.direction.y + planes.at(i).C * ray.direction.z;
        if(denom == 0){
            std::cout << "denom is zero\n";
            return background;
        }
        float numer = -1.0 * (planes.at(i).A * ray.position.x + planes.at(i).B * ray.position.y + planes.at(i).C * ray.position.z + planes.at(i).D);
        float t3 = numer/denom;
        if(t3 > 0 && t3 < min_d){
            glm::vec3 p;
            p.x = ray.position.x + t3 * ray.direction.x;
            p.y = ray.position.y + t3 * ray.direction.y;
            p.z = ray.position.z + t3 * ray.direction.z;
            glm::vec3 ep = p - planes.at(i).triangle.p0;
            float d11 = glm::dot(planes.at(i).triangle.e1, planes.at(i).triangle.e1);
            float d22 = glm::dot(planes.at(i).triangle.e2, planes.at(i).triangle.e2);
            float d12 = glm::dot(planes.at(i).triangle.e1, planes.at(i).triangle.e2);
            float d1p = glm::dot(planes.at(i).triangle.e1, ep);
            float d2p = glm::dot(planes.at(i).triangle.e2, ep);
            float det = (d11 * d22 - d12 * d12);
            if(det == 0){
                std::cout << "det is zero\n";
            }
            float beta = (d22 * d1p - d12 * d2p)/det;
            float gamma = (d11 * d2p - d12 * d1p)/det;
            float alpha = 1 - (beta + gamma);
            planes.at(i).triangle.alpha = alpha;
            planes.at(i).triangle.beta = beta;
            planes.at(i).triangle.gamma = gamma;
            if(beta > 0 && beta < 1 && gamma > 0 && gamma < 1 && (beta + gamma) > 0 && (beta + gamma) < 1){
                min_d = t3;
                object_type = 2;
                current_obj = i;
            }
        }
    }
    if(current_obj != -1 && object_type == 1){  //sphere
        intersectionPoint = min_d * ray.direction + ray.position;
        if(objects.at(current_obj).m.justrgb == false){
            return Shade_ray2(objects.at(current_obj), intersectionPoint, ray, current_obj);
        }
        else return Shade_ray(objects.at(current_obj));
    }
    else if(current_obj != -1 && object_type == 2){  //triangle
        intersectionPoint = min_d * ray.direction + ray.position;
        return Shade_Ray_Triangle2(planes.at(current_obj), intersectionPoint, ray, current_obj, min_d);
    }
    return background;
}

bool is_number(const std::string& s){ //determines if a string is a number.
    for(char const &ch : s){
        if(std::isdigit(ch) == 0 && ch != '-' && ch != '.'){
            return false;
        }
    }
    return true;
}

int raycast(std::string filename){
    bool initial_texture_status = false;
    std::ifstream file;
    file.open(filename);
    if(file.fail()){
        std::cout << " fail\n";
        return 0;
    }
    if(filename.substr(filename.size() - 4) != ".txt"){
        std::cout << "not a text file!\n";
        return 0;
    }
    glm::vec3 dummy;    //make points start at index 1
    dummy.x = 0;
    dummy.y = 0;
    dummy.z = 0;
    glm::vec2 dummy2;
    dummy2.x = 0;
    dummy2.y = 0;
    points.push_back(dummy);
    shading_Points.push_back(dummy);
    texture_Points.push_back(dummy2);
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::vector<MaterialType> object_materials;
    glm::vec3 eye;
    glm::vec3 vdir;
    glm::vec3 up;
    FovType fov;
    ImsizeType imsize;
    ColorType bkgcolor;
    ViewingwindowType window;
    RayType ray;
    bool skip = false;
    std::string subs;
    while(buffer){  // 1. read information from text file
        if(skip == false){
            buffer >> subs;
        }
        skip = false;
        if(subs == "eye"){  // get eye information.
            buffer >> subs;
            if(is_number(subs)){
                eye.x = stof(subs);
            }
            else{
                return 0;
            }
            buffer >> subs;
            if(is_number(subs)){
                eye.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                eye.z = stof(subs);
            }
            else return 0;
        }
        else if(subs == "viewdir"){  // get viewdir information.
            buffer >> subs;
            if(is_number(subs)){
                vdir.x = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                vdir.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                vdir.z = stof(subs);
            } 
            else return 0;
        }
        else if(subs == "updir"){    // get updir information.
            buffer >> subs;
            if(is_number(subs)){
                up.x = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                up.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                up.z = stof(subs);
            }
            else return 0;
        }
        else if(subs == "hfov"){ // get fov information.
            buffer >> subs;
            if(is_number(subs)){
                fov.h = stof(subs);
            }
            else return 0;
        }
        else if(subs == "imsize"){   // get size information
            buffer >> subs;
            if(is_number(subs)){
                imsize.width = stoi(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                imsize.height = stoi(subs);
            }
            else return 0;
        }
        else if(subs == "bkgcolor"){   // get background color information
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                bkgcolor.r = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                bkgcolor.g = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                bkgcolor.b = stof(subs);
            }
            else return 0;
        }
        else if(subs == "mtlcolor"){   // get background color information
            buffer >> subs;
            MaterialType newmat;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.od.x = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.od.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.od.z = stof(subs);
            }
            else return 0;
            buffer >> subs; //reads in sphere
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.os.x = stof(subs);
                newmat.justrgb = false;
            }
            else{
                skip = true;
                goto onlyRGB;   
            }
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.os.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.os.z = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.ka = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.kd = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newmat.ks = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                newmat.n = stoi(subs);
            }
            else return 0;
            onlyRGB:
            object_materials.push_back(newmat);
        }
        else if(subs == "sphere"){
            buffer >> subs;
            SphereType newsphere;
            if(is_number(subs)){
                newsphere.position.x = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                newsphere.position.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                newsphere.position.z = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                newsphere.radius = stof(subs);
            }
            else return 0;

            newsphere.m = object_materials.back();
            if(initial_texture_status == true){
                newsphere.textureStatus = true;
                newsphere.texture = textures.back();
            }
            objects.push_back(newsphere);
        }
        else if(subs == "light"){
            buffer >> subs;
            LightType newlight;
            if(is_number(subs)){
                newlight.position.x = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                newlight.position.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                newlight.position.z = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stoi(subs) == 1 || stoi(subs) == 0){
                newlight.type = stoi(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newlight.color.x = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newlight.color.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs) && stof(subs) <= 1 && stof(subs) >= 0){
                newlight.color.z = stof(subs);
            }
            else return 0;
            lights.push_back(newlight);
        }
        else if(subs == "v"){
            buffer >> subs;
            glm::vec3 vecToAdd;
            if(is_number(subs)){
                vecToAdd.x = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                vecToAdd.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                vecToAdd.z = stof(subs);
            }
            else return 0;
            points.push_back(vecToAdd);
        } 
        else if(subs == "f"){
            buffer >> subs;
            FaceType faceToAdd;
            std::string delimiter;
            std::string firstNum;
            std::string secondNum;
            if(subs.find("//") != std::string::npos){  //smooth shaded
                faceToAdd.SS = true;
                delimiter = "//";
            }
            else if(subs.find("/") != std::string::npos){    //flat shaded and textured
                faceToAdd.textureStatus = true;
                faceToAdd.texture = textures.back();
                delimiter = "/";
            }
            if(is_number(subs)){    //flat shading no texture
                faceToAdd.position.x = stoi(subs);
            }
            else if(faceToAdd.SS == true){  //smooth shaded
                firstNum = subs.substr(0, subs.find(delimiter));
                faceToAdd.position.x = stoi(firstNum);
                secondNum = subs.substr(subs.find(delimiter) + 2, subs.length() - 1);
                faceToAdd.index1 = stoi(secondNum);
            }
            else if(faceToAdd.textureStatus == true){  //flat shaded with texture
                firstNum = subs.substr(0, subs.find(delimiter));
                faceToAdd.position.x = stoi(firstNum);
                secondNum = subs.substr(subs.find(delimiter) + 1, subs.length() - 1);
                faceToAdd.index1 = stoi(secondNum);
            }
            buffer >> subs;
            if(is_number(subs)){
                faceToAdd.position.y = stoi(subs);
            }
            else if(faceToAdd.SS == true){
                firstNum = subs.substr(0, subs.find(delimiter));
                faceToAdd.position.y = stoi(firstNum);
                secondNum = subs.substr(subs.find(delimiter) + 2, subs.length() - 1);
                faceToAdd.index2 = stoi(secondNum);
            }
            else if(faceToAdd.textureStatus == true){
                firstNum = subs.substr(0, subs.find(delimiter));
                faceToAdd.position.y = stoi(firstNum);
                secondNum = subs.substr(subs.find(delimiter) + 1, subs.length() - 1);
                faceToAdd.index2 = stoi(secondNum);
            }
            buffer >> subs;
            if(is_number(subs)){
                faceToAdd.position.z = stoi(subs);
            }
            else if(faceToAdd.SS == true){
                firstNum = subs.substr(0, subs.find(delimiter));
                faceToAdd.position.z = stoi(firstNum);
                secondNum = subs.substr(subs.find(delimiter) + 2, subs.length() - 1);
                faceToAdd.index3 = stoi(secondNum);
            }
            else if(faceToAdd.textureStatus == true){
                firstNum = subs.substr(0, subs.find(delimiter));
                faceToAdd.position.z = stoi(firstNum);
                secondNum = subs.substr(subs.find(delimiter) + 1, subs.length() - 1);
                faceToAdd.index3 = stoi(secondNum);
            }
            faceToAdd.m = object_materials.back();
            faces.push_back(faceToAdd);
        }
        else if(subs == "vn"){
            buffer >> subs;
            glm::vec3 shadePt;
            if(is_number(subs)){
                shadePt.x = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                shadePt.y = stof(subs);
            }
            else return 0;
            buffer >> subs;
            if(is_number(subs)){
                shadePt.z = stof(subs);
            }
            else return 0;  
            shading_Points.push_back(glm::normalize(shadePt));
        }
        else if(subs == "texture"){ 
            initial_texture_status = true;
            TextureType texturetopush; 
            std::stringstream tbuffer;
            std::string tsubs;  
            buffer >> subs;
            std::ifstream file;
            file.open(subs);
            if(file.fail()){
                std::cout << "file fail\n";
                return 0;
            }
            if(subs.substr(subs.size() - 4) != ".ppm"){
                std::cout << "not a ppm file!\n";
                return 0;
            }
            tbuffer << file.rdbuf();  
            tbuffer >> tsubs;
            if(tsubs != "P3"){
                std::cout << "not a p3 file\n";
                return -1;
            }
            tbuffer >> tsubs;
            if(is_number(tsubs)){
                texturetopush.width = stof(tsubs);
            }
            else return 0;
            tbuffer >> tsubs; 
            if(is_number(tsubs)){
                texturetopush.height = stof(tsubs);
            }
            else return 0;
            tbuffer >> tsubs;
            tbuffer >> tsubs;
            glm::vec3** image = new glm::vec3* [texturetopush.height];       //initialize
            for(int j = 0; j < texturetopush.height -1; j++){
                image[j] = new glm::vec3[texturetopush.width];
            }
            glm::vec3 basic;
            for(int i = 0; i < texturetopush.height -1; i++){
                for(int j = 0; j < texturetopush.width - 1; j++){
                    basic.x = stof(tsubs);
                    tbuffer >> tsubs;
                    basic.y = stof(tsubs);
                    tbuffer >> tsubs;
                    basic.z = stof(tsubs);
                    tbuffer >> tsubs;
                    image[i][j] = basic;
                }
            }
            texturetopush.texturemap = image;
            textures.push_back(texturetopush);
        }
        else if(subs == "vt"){
            glm::vec2 texture_Point_ToAdd;
            buffer >> subs;
            if(is_number(subs)){
                texture_Point_ToAdd.x = stoi(subs);
            }
            else return -1;
            buffer >> subs;
            if(is_number(subs)){
                texture_Point_ToAdd.y = stoi(subs);
            }
            texture_Points.push_back(texture_Point_ToAdd);
        }
    }
    Populate_planes(points, faces);
    ColorType arr[imsize.width * imsize.height];  // 2. Define an array of sufficent size to store the color values
    float aspect_ratio = (float)imsize.width/(float)imsize.height;  // define viewing window
    window.width = 2 * 5 * (tan(.5 * (fov.h * M_PI / 180.0)));
    window.height = window.width / aspect_ratio;
    glm::vec3 u = glm::cross(vdir, up); // calculate u and v
    glm::vec3 v = glm::cross(u, vdir);
    u = glm::normalize(u);
    v = glm::normalize(v);
    if(u.x == -0){     //cross managed to compute -0, this is to fix incase.
        u.x = 0;
    }
    if(u.y == -0){
        u.y = 0;
    }
    if(u.z == -0){
        u.z = 0;
    }
    if(v.x == -0){
        v.x = 0;
    }
    if(v.y == -0){
        u.y = 0;
    }
    if(v.z == -0){
        v.z = 0;
    }
    vdir = glm::normalize(vdir);
    window.ul = eye + vdir * (float)5 - (u * (float)(window.width/2)) + (v * (float)(window.height/2));;    // calculate four courners of viewing window
    window.ur = eye + vdir * (float)5 + (u * (float)(window.width/2)) + (v * (float)(window.height/2));;
    window.ll = eye + vdir * (float)5 - (u * (float)(window.width/2)) - (v * (float)(window.height/2));;
    window.lr = eye + vdir * (float)5 + (u * (float)(window.width/2)) - (v * (float)(window.height/2));;
    glm::vec3 deltaH;
    glm::vec3 deltaV;
    deltaH = (window.ur - window.ul) / (float)(imsize.width -1);
    deltaV = (window.ll - window.ul) / (float)(imsize.height -1);
    glm::vec3 point;
    std::string newfile = filename; // A string storing the output file's name is created and assigned with the previous file name, .txt is removed and then out.txt is appended.
    newfile.pop_back();
    newfile.pop_back();
    newfile.pop_back();
    newfile.pop_back();
    newfile.append(".ppm");
    std::ofstream outfile(newfile); // create new file with our created name.
    outfile << "P3" << "\n" << imsize.width << " " << imsize.height << "\n" << "255\n"; //write header
    ray.position = eye;
    for(int j = 0; j < imsize.height; j++){  // 4. for each pixel in the output image
        for(int i = 0; i < imsize.width; i++){
            point = window.ul + (float)i * deltaH + (float)j * deltaV;
            ray.direction = glm::normalize((point - eye));
            arr[i] = Trace_Ray(ray, objects, bkgcolor);
            outfile << int(arr[i].r * 255) << " " << int(arr[i].g * 255) << " " << int(arr[i].b * 255) << "\n";
        }
    }
    outfile.close();
    return 1;
}

int main(int argc, char **argv){
    if(argc == 2){
    raycast(argv[1]);
    }
    else{
        std::cout << "One file needed\n";
        return 0;
    }
}
