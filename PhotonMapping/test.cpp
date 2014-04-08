#include </Users/terence/dev/cprocessing/cprocessing.hpp>
#include <vector>

using namespace std;
using namespace cprocessing;

// ----- Scene Description -----
int szImg = 512;                  //Image Size
int nrTypes = 2;                  //2 Object Types (Sphere = 0, Plane = 1)
vector<int> nrObjects = {2,5};          //2 Spheres, 5 Planes
float gAmbient = 0.1;             //Ambient Lighting
vector<float> gOrigin = {0.0,0.0,0.0};  //World Origin for Convenient Re-Use Below (Constant)
vector<float> Light = {0.0,1.2,3.75};   //Point Light-Source Position
vector<vector<float>> spheres = {{1.0,0.0,4.0,0.5},{-0.6,-1.0,4.5,0.5}};         //Sphere Center & Radius
vector<vector<float>> planes  = {{0, 1.5},{1, -1.5},{0, -1.5},{1, 1.5},{2,5.0}}; //Plane Axis & Distance-to-Origin

// ----- Photon Mapping -----
int nrPhotons = 1000;             //Number of Photons Emitted
int nrBounces = 3;                //Number of Times Each Photon Bounces
bool lightPhotons = true;      //Enable Photon Lighting?
float sqRadius = 0.7;             //Photon Integration Area (Squared for Efficiency)
float exposure = 50.0;            //Number of Photons Integrated at Brightest Pixel
vector<vector<int> > numPhotons = {{0,0},{0,0,0,0,0}};              //Photon Count for Each Scene Object

vector<vector<vector<vector<vector<float> > > > > photons(2,
    vector<vector<vector<vector<float> > > >(5, vector<vector<vector<float> > >(5000,
        vector<vector<float> >(3, vector<float>(3))))); //Allocated Memory for Per-Object Photon Info

// ----- Raytracing Globals -----
bool gIntersect = false;       //For Latest Raytracing Call... Was Anything Intersected by the Ray?
int gType;                        //... Type of the Intersected Object (Sphere or Plane)
int gIndex;                       //... Index of the Intersected Object (Which Sphere/Plane Was It?)
float gSqDist, gDist = -1.0;      //... Distance from Ray Origin to Intersection
vector<float> gPoint = {0.0, 0.0, 0.0}; //... Point At Which the Ray Intersected the Object


bool empty = true;
int pRow, pCol, pIteration, pMax;     //Pixel Rendering Order
bool odd(int x) {return x % 2 != 0;}


float min(float a, float b) {
    return a <= b ? a : b;
}

float dot3(vector<float> a, vector<float> b){     //Dot Product 3-Vectors
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

vector<float> mul3c(vector<float> a, float c){    //Multiply 3-Vector with Scalar
  vector<float> result = {c*a[0], c*a[1], c*a[2]};
  return result;
}

vector<float> normalize3(vector<float> v){        //Normalize 3-Vector

    float L = sqrt(dot3(v,v));
    return mul3c(v, 1.0/L);
}

vector<float> sub3(vector<float> a, vector<float> b){   //Subtract 3-Vectors
    vector<float> result = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return result;
}

vector<float> add3(vector<float> a, vector<float> b){   //Add 3-Vectors
    vector<float> result = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    return result;
}

float random(float Min, float Max)
{
    return ((float(rand()) / float(RAND_MAX)) * (Max - Min)) + Min;
}

vector<float> rand3(float s){               //Random 3-Vector
    vector<float> rand = {random(-s,s),random(-s,s),random(-s,s)};
    return rand;
}

bool gatedSqDist3(vector<float> a, vector<float> b, float sqradius){ //Gated Squared Distance
    float c = a[0] - b[0];          //Efficient When Determining if Thousands of Points
    float d = c*c;                  //Are Within a Radius of a Point (and Most Are Not!)
    if (d > sqradius) return false; //Gate 1 - If this dimension alone is larger than
    c = a[1] - b[1];                //         the search radius, no need to continue
    d += c*c;
    if (d > sqradius) return false; //Gate 2
    c = a[2] - b[2];
    d += c*c;
    if (d > sqradius) return false; //Gate 3
    gSqDist = d;      return true ; //Store Squared Distance Itself in Global State
}

//---------------------------------------------------------------------------------------
//Ray-Geometry Intersections  -----------------------------------------------------------
//---------------------------------------------------------------------------------------

void checkDistance(float lDist, int p, int i){
    if (lDist < gDist && lDist > 0.0){ //Closest Intersection So Far in Forward Direction of Ray?
        gType = p; gIndex = i; gDist = lDist; gIntersect = true;} //Save Intersection in Global State
}

void raySphere(int idx, vector<float> r, vector<float> o) //Ray-Sphere Intersection: r=Ray Direction, o=Ray Origin
{
    vector<float> s = sub3(spheres[idx],o);  //s=Sphere Center Translated into Coordinate Frame of Ray Origin
    float radius = spheres[idx][3];    //radius=Sphere Radius

    //Intersection of Sphere and Line     =       Quadratic Function of Distance
    float A = dot3(r,r);                       // Remember This From High School? :
    float B = -2.0 * dot3(s,r);                //    A x^2 +     B x +               C  = 0
    float C = dot3(s,s) - sq(radius);          // (r'r)x^2 - (2s'r)x + (s's - radius^2) = 0
    float D = B*B - 4*A*C;                     // Precompute Discriminant

    if (D > 0.0){                              //Solution Exists only if sqrt(D) is Real (not Imaginary)
        float sign = (C < -0.00001) ? 1 : -1;    //Ray Originates Inside Sphere If C < 0
        float lDist = (-B + sign*sqrt(D))/(2*A); //Solve Quadratic Equation for Distance to Intersection
        checkDistance(lDist,0,idx);}             //Is This Closest Intersection So Far?
}

void rayPlane(int idx, vector<float> r, vector<float> o){ //Ray-Plane Intersection
    int axis = (int) planes[idx][0];            //Determine Orientation of Axis-Aligned Plane
    if (r[axis] != 0.0){                        //Parallel Ray -> No Intersection
        float lDist = (planes[idx][1] - o[axis]) / r[axis]; //Solve Linear Equation (rx = p-o)
        checkDistance(lDist,1,idx);}
}

void rayObject(int type, int idx, vector<float> r, vector<float> o){
    if (type == 0) raySphere(idx,r,o); else rayPlane(idx,r,o);
}


//---------------------------------------------------------------------------------------
// Lighting -----------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

float lightDiffuse(vector<float> N, vector<float> P){  //Diffuse Lighting at Point P with Surface Normal N
    vector<float> L = normalize3( sub3(Light,P) ); //Light Vector (Point to Light)
    return dot3(N,L);                        //Dot Product = cos (Light-to-Surface-Normal Angle)
}

vector<float> sphereNormal(int idx, vector<float> P){
    return normalize3(sub3(P,spheres[idx])); //Surface Normal (Center to Point)
}

vector<float> planeNormal(int idx, vector<float> P, vector<float> O){
    int axis = (int) planes[idx][0];
    vector<float> N = {0.0,0.0,0.0};
    N[axis] = O[axis] - planes[idx][1];      //Vector From Surface to Light
    return normalize3(N);
}

vector<float> surfaceNormal(int type, int index, vector<float> P, vector<float> Inside){
    if (type == 0) {return sphereNormal(index,P);}
    else           {return planeNormal(index,P,Inside);}
}

float lightObject(int type, int idx, vector<float> P, float lightAmbient){
    float i = lightDiffuse( surfaceNormal(type, idx, P, Light) , P );
    return min(1.0, cprocessing::max(i, lightAmbient));   //Add in Ambient Light by Constraining Min Value
}

//---------------------------------------------------------------------------------------
// Raytracing ---------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

void raytrace(vector<float> ray, vector<float> origin)
{
    gIntersect = false; //No Intersections Along This Ray Yet
    gDist = 999999.9;   //Maximum Distance to Any Object

    for (int t = 0; t < nrTypes; t++)
        for (int i = 0; i < nrObjects[t]; i++)
            rayObject(t,i,ray,origin);
}

vector<float> reflect(vector<float> ray, vector<float> fromPoint){                //Reflect Ray
    vector<float> N = surfaceNormal(gType, gIndex, gPoint, fromPoint);  //Surface Normal
    return normalize3(sub3(ray, mul3c(N,(2 * dot3(ray,N)))));     //Approximation to Reflection
}





//---------------------------------------------------------------------------------------
//Photon Mapping ------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

vector<float> gatherPhotons(vector<float> p, int type, int id){
    vector<float> energy = {0.0,0.0,0.0};
    vector<float> N = surfaceNormal(type, id, p, gOrigin);                   //Surface Normal at Current Point
    for (int i = 0; i < numPhotons[type][id]; i++){                    //Photons Which Hit Current Object
        if (gatedSqDist3(p,photons[type][id][i][0],sqRadius)){           //Is Photon Close to Point?
            float weight = cprocessing::max(0.0f, -dot3(N, photons[type][id][i][1] ));   //Single Photon Diffuse Lighting
            weight *= (1.0 - sqrt(gSqDist)) / exposure;                    //Weight by Photon-Point Distance
            energy = add3(energy, mul3c(photons[type][id][i][2], weight)); //Add Photon's Energy to Total
        }}
    return energy;
}

void storePhoton(int type, int id, vector<float> location, vector<float> direction, vector<float> energy){
    photons[type][id][numPhotons[type][id]][0] = location;  //Location
    photons[type][id][numPhotons[type][id]][1] = direction; //Direction
    photons[type][id][numPhotons[type][id]][2] = energy;    //Attenuated Energy (Color)
    numPhotons[type][id]++;
}

void shadowPhoton(vector<float> ray){                               //Shadow Photons
    vector<float> shadow = {-0.25,-0.25,-0.25};
    vector<float> tPoint = gPoint;
    int tType = gType, tIndex = gIndex;                         //Save State
    vector<float> bumpedPoint = add3(gPoint,mul3c(ray,0.00001));      //Start Just Beyond Last Intersection
    raytrace(ray, bumpedPoint);                                 //Trace to Next Intersection (In Shadow)
    vector<float> shadowPoint = add3( mul3c(ray,gDist), bumpedPoint); //3D Point
    storePhoton(gType, gIndex, shadowPoint, ray, shadow);
    gPoint = tPoint; gType = tType; gIndex = tIndex;            //Restore State
}

vector<float> filterColor(vector<float> rgbIn, float r, float g, float b){ //e.g. White Light Hits Red Wall
    vector<float> rgbOut = {r,g,b};
    for (int c=0; c<3; c++) rgbOut[c] = min(rgbOut[c],rgbIn[c]); //Absorb Some Wavelengths (R,G,B)
    return rgbOut;
}

vector<float> getColor(vector<float> rgbIn, int type, int index){ //Specifies Material Color of Each Object
    if      (type == 1 && index == 0) { return filterColor(rgbIn, 0.0, 1.0, 0.0);}
    else if (type == 1 && index == 2) { return filterColor(rgbIn, 1.0, 0.0, 0.0);}
    else                              { return filterColor(rgbIn, 1.0, 1.0, 1.0);}
}


vector<float> computePixelColor(float x, float y){
    vector<float> rgb = {0.0,0.0,0.0};
    vector<float> ray = {  (float)(x/szImg - 0.5) ,       //Convert Pixels to Image Plane Coordinates
        -(float)(y/szImg - 0.5), 1.0f}; //Focal Length = 1.0
    raytrace(ray, gOrigin);                //Raytrace!!! - Intersected Objects are Stored in Global State

    if (gIntersect){                       //Intersection
        gPoint = mul3c(ray,gDist);           //3D Point of Intersection

        if (gType == 0 && gIndex == 1){      //Mirror Surface on This Specific Object
            ray = reflect(ray,gOrigin);        //Reflect Ray Off the Surface
            raytrace(ray, gPoint);             //Follow the Reflected Ray
            if (gIntersect){ gPoint = add3( mul3c(ray,gDist), gPoint); }} //3D Point of Intersection

        if (lightPhotons){                   //Lighting via Photon Mapping
            rgb = gatherPhotons(gPoint,gType,gIndex);}
        else{                                //Lighting via Standard Illumination Model (Diffuse + Ambient)
            int tType = gType, tIndex = gIndex;//Remember Intersected Object
            float i = gAmbient;                //If in Shadow, Use Ambient Color of Original Object
            raytrace( sub3(gPoint,Light) , Light);  //Raytrace from Light to Object
            if (tType == gType && tIndex == gIndex) //Ray from Light->Object Hits Object First?
                i = lightObject(gType, gIndex, gPoint, gAmbient); //Not In Shadow - Compute Lighting
            rgb[0]=i; rgb[1]=i; rgb[2]=i;
            rgb = getColor(rgb,tType,tIndex);}
    }
    return rgb;
}

void drawPhoton(vector<float> rgb, vector<float> p){           //Photon Visualization
    if (p[2] > 0.0){                       //Only Draw if In Front of Camera
        int x = (szImg/2) + (int)(szImg *  p[0]/p[2]); //Project 3D Points into Scene
        int y = (szImg/2) + (int)(szImg * -p[1]/p[2]); //Don't Draw Outside Image
        if (y <= szImg) {stroke(255.0*rgb[0],255.0*rgb[1],255.0*rgb[2]); point(x,y);}}
}

void emitPhotons(){
//    randomSeed(0);                               //Ensure Same Photons Each Time
    for (int t = 0; t < nrTypes; t++)            //Initialize Photon Count to Zero for Each Object
        for (int i = 0; i < nrObjects[t]; i++)
            numPhotons[t][i] = 0;

    for (int i = 0; i < nrPhotons * 3.0; i++){
        int bounces = 1;
        vector<float> rgb = {1.0,1.0,1.0};               //Initial Photon Color is White
        vector<float> ray = normalize3( rand3(1.0) );    //Randomize Direction of Photon Emission
        vector<float> prevPoint = Light;                 //Emit From Point Light Source

        //Spread Out Light Source, But Don't Allow Photons Outside Room/Inside Sphere
        while (prevPoint[1] >= Light[1]){ prevPoint = add3(Light, mul3c(normalize3(rand3(1.0)), 0.75));}
        if (abs(prevPoint[0]) > 1.5 || abs(prevPoint[1]) > 1.2 ||
            gatedSqDist3(prevPoint,spheres[0],spheres[0][3]*spheres[0][3])) bounces = nrBounces+1;

        raytrace(ray, prevPoint);                          //Trace the Photon's Path

        while (gIntersect && bounces <= nrBounces){        //Intersection With New Object
            gPoint = add3( mul3c(ray,gDist), prevPoint);   //3D Point of Intersection
            rgb = mul3c (getColor(rgb,gType,gIndex), 1.0/sqrt(bounces));
            storePhoton(gType, gIndex, gPoint, ray, rgb);  //Store Photon Info
            drawPhoton(rgb, gPoint);                       //Draw Photon
            shadowPhoton(ray);                             //Shadow Photon
            ray = reflect(ray,prevPoint);                  //Bounce the Photon
            raytrace(ray, gPoint);                         //Trace It to Next Location
            prevPoint = gPoint;
            bounces++;}
    }
}




void drawInterface() {
//    String path = "/Users/terence/dev/PhotonMapping/processing/data/";
//    stroke(221,221,204); fill(221,221,204); rect(0,szImg,szImg,48); //Fill Background with Page Color
//    img1=loadImage(path + "1_32.png"); img2=loadImage(path + "2_32.png"); img3=loadImage(path + "3_32.png"); //Load Images

//    textFont(font); //Display Text
//    if (!view3D) {fill(0); img3.filter(GRAY);} else fill(160);
//    text("Ray Tracing", 64, szImg + 28);
//    if (lightPhotons || view3D) {fill(0); img1.filter(GRAY);} else fill(160);
//    text("Photon Mapping", 368, szImg + 28);
//    if (!lightPhotons || view3D) img2.filter(GRAY);

//    stroke(0); fill(255);  //Draw Buttons with Icons
//    rect(198,519,33,33); image(img1,199,520);
//    rect(240,519,33,33); image(img2,241,520);
//    rect(282,519,33,33); image(img3,283,520);
}

void switchToMode(char i, int x){

    drawInterface();
}

void render(){ //Render Several Lines of Pixels at Once Before Drawing
    int x,y,iterations = 0;
    vector<float> rgb = {0.0,0.0,0.0};

    if (pCol >= pMax) {pRow++; pCol = 0;
        if (pRow >= pMax) {pIteration++; pRow = 0; pMax = int(pow(2,pIteration));}}
    bool pNeedsDrawing = (pIteration == 1 || odd(pRow) || (!odd(pRow) && odd(pCol)));
    x = pCol * (szImg/pMax); y = pRow * (szImg/pMax);
    pCol++;

    if (pNeedsDrawing){
        iterations++;
        rgb = mul3c( computePixelColor(x,y), 255.0);               //All the Magic Happens in Here!
        stroke(rgb[0],rgb[1],rgb[2]); cprocessing::fill(rgb[0],rgb[1],rgb[2]);  //Stroke & Fill
        rect(x,y,(szImg/pMax)-1,(szImg/pMax)-1);}                  //Draw the Possibly Enlarged Pixel

    if (pRow == szImg-1) {empty = false;}
}


void setup(){
    size(szImg,szImg + 48);
    emitPhotons();
}

void draw(){

        if (empty){
            stroke(0); fill(0); rect(0,0,szImg-1,szImg-1); //Black Out Drawing Area
            emitPhotons(); empty = false;
        }

}

