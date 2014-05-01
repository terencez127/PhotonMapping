#include "main.h"
#include <vector>
#include "GLUT/glut.h"
#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;
using namespace pm;

int mainWin, subWin1, subWin2, currentWindow;
int mouseX = 0;
int mouseY = 0;
bool bMousePressed = false;
bool bKeyPressed = false;
unsigned char key = 255;
int keyCode;

int width;
int height;

int frameRate = 60;

color strokeColor (0,0,0);
color fillColor   (255,255,255);

bool emitDone = false;

float min(float a, float b) {
    return a <= b ? a : b;
}

//Dot Product 3-Vectors
float dot3(vector<float> const &a, vector<float> const &b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

//Multiply 3-Vector with Scalar
vector<float> mul3c(vector<float> const &a, float c){
    vector<float> result = {c*a[0], c*a[1], c*a[2]};
    return result;
}

//Normalize 3-Vector
vector<float> normalize3(vector<float> const &v){
    float L = sqrt(dot3(v,v));
    return mul3c(v, 1.0/L);
}

//Subtract 3-Vectors
vector<float> sub3(vector<float> const &a, vector<float> const &b){
    vector<float> result = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return result;
}

//Add 3-Vectors
vector<float> add3(vector<float> const &a, vector<float> const &b){
    vector<float> result = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    return result;
}

float random(float Min, float Max)
{
    return ((float(rand()) / float(RAND_MAX)) * (Max - Min)) + Min;
}

//Random 3-Vector
vector<float> rand3(float s){
    vector<float> rand = {random(-s,s),random(-s,s),random(-s,s)};
    return rand;
}

template<class C>
C constrain(const C& a, const C& minv, const C& maxv) { return min(maxv,max(minv,a)); }


double clamp(double v) {
    return v > 255 ? 255 : v < 0 ? 0 :  v;
}

pm::color::color(double val1, double val2, double val3, double valA) {
    rgba[0] = clamp(val1);
    rgba[1] = clamp(val2);
    rgba[2] = clamp(val3);
    rgba[3] = clamp(valA);
};


color::color(double gray, double alpha){
    unsigned char val = clamp(gray);
    rgba[0] = val;
    rgba[1] = val;
    rgba[2] = val;
    rgba[3] = clamp(alpha);
}

static std::vector<unsigned> sphereIdx;


void quad (double x0, double y0,
           double x1, double y1,
           double x2, double y2,
           double x3, double y3)
{
    GLdouble vertices[] = {
        x0, y0,
        x1, y1,
        x2, y2,
        x3, y3
    };
    // activate and specify pointer to vertex array
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_DOUBLE, 0, vertices);

    if (fillColor.rgba[3] > 0) {
        // See if filled triangle is required
        glColor4ubv (fillColor.rgba);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDrawArrays(GL_QUADS,0,4);
    }
    if (strokeColor.rgba[3] > 0) {
        // See if outline triangle is required
        glColor4ubv (strokeColor.rgba);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawArrays(GL_QUADS,0,4);
    }
    // deactivate vertex arrays after drawing
    glDisableClientState(GL_VERTEX_ARRAY);

}

/// Draws a point.
void drawPoint (double x, double y, double z = 0)
{
    if (strokeColor.rgba[3] > 0) {
        // Draw point using the stroke color
        glColor4ubv (strokeColor.rgba);
        glBegin (GL_POINTS);
        glVertex3d (x,y,z);
        glEnd();
    }
}

void pm::rect (double x, double y, double a, double b)
{
    quad (x, y, x+a, y, x+a, y+b, x, y+b);
}

void background (const color& c) {
    glClearColor (c.rgba[0] * (1.0/255),
                  c.rgba[1] * (1.0/255),
                  c.rgba[2] * (1.0/255),
                  c.rgba[3] * (1.0/255));
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}



int resolution = 512;
int nrTypes = 2;
vector<int> nrObjects = {3,5};
vector<float> gOrigin = {0.0,0.0,0.0};
vector<float> Light = {0.0,1.2,3.75};
vector<vector<float> > spheres = {{1.0,0.0,4.0,0.5}, {-0.6,-1.0,4.5,0.5}, {0.0,1.0,4.5,0.3}};
vector<vector<float>> planes  = {{0, 1.5},{1, -1.5},{0, -1.5},{1, 1.5},{2,5.0}};


int nrPhotons = 400;
int nrBounces = 3;
float sqRadius = 0.7;
float exposure = 30.0;
vector<vector<int> > numPhotons = {{0,0,0},{0,0,0,0,0}};

//  photon map, [type][index][photon][photon_info][]
vector<vector<vector<vector<vector<float> > > > > photons
(2,vector<vector<vector<vector<float> > > >
 (6, vector<vector<vector<float> > >
  (5000,vector<vector<float> >
   (3, vector<float>(3)))));

bool gIntersect = false;
int gType;  //  Type of the Intersected Object (Sphere or Plane)
int gIndex; //  Index of the Intersected Object
float gSqDist, gDist = -1.0;    //  Distance from Ray Origin to Intersection
vector<float> gPoint = {0.0, 0.0, 0.0}; //  The Intersection point

bool empty = true;
bool mode3D = true;
int pRow, pCol, pIteration, pMax;
bool odd(int x) {return x % 2 != 0;}


int prevMouseX = -9999, prevMouseY = -9999, sphereIndex = -1;
float s = 130.0;
bool mouseDragging = false;
void mouseReleased() {
    prevMouseX = -9999;
    prevMouseY = -9999;
    mouseDragging = false;

}


bool gatedSqDist3(vector<float> const &a, vector<float> const &b, float sqradius){
    float c = a[0] - b[0];
    float d = c*c;
    if (d > sqradius) return false;
    c = a[1] - b[1];
    d += c*c;
    if (d > sqradius) return false;
    c = a[2] - b[2];
    d += c*c;
    if (d > sqradius) return false;
    gSqDist = d;      return true ;
}



void checkDistance(float lDist, int p, int i){
    if (lDist < gDist && lDist > 0.0){
        gType = p; gIndex = i; gDist = lDist; gIntersect = true;}
}

void raySphere(int idx, vector<float> const &direction, vector<float> const &origin)
{
    vector<float> s = sub3(spheres[idx],origin);
    float radius = spheres[idx][3];


    float A = dot3(direction,direction);
    float B = -2.0 * dot3(s,direction);
    float C = dot3(s,s) - radius*radius;
    float D = B*B - 4*A*C;

    if (D > 0.0){
        float sign = (C < -0.00001) ? 1 : -1;
        float lDist = (-B + sign*sqrt(D))/(2*A);
        checkDistance(lDist,0,idx);}
}

void rayPlane(int idx, vector<float> const &r, vector<float> const &o){
    int axis = (int) planes[idx][0];
    if (r[axis] != 0.0){
        float lDist = (planes[idx][1] - o[axis]) / r[axis];
        checkDistance(lDist,1,idx);}
}

void rayObject(int type, int idx, vector<float> const &r, vector<float> const &o){
    if (type == 0) raySphere(idx,r,o); else rayPlane(idx,r,o);
}



float lightDiffuse(vector<float> const &N, vector<float> const &P){
    vector<float> L = normalize3( sub3(Light,P) );
    return dot3(N,L);
}

vector<float> sphereNormal(int idx, vector<float>& P){
    return normalize3(sub3(P,spheres[idx]));
}

vector<float> planeNormal(int idx, vector<float>& P, vector<float>& O){
    int axis = (int) planes[idx][0];
    vector<float> N = {0.0,0.0,0.0};
    N[axis] = O[axis] - planes[idx][1];
    return normalize3(N);
}

vector<float> surfaceNormal(int type, int index, vector<float>& P, vector<float>& Inside){
    if (type == 0) {return sphereNormal(index,P);}
    else           {return planeNormal(index,P,Inside);}
}

float lightObject(int type, int idx, vector<float>& P, float lightAmbient){
    float i = lightDiffuse( surfaceNormal(type, idx, P, Light) , P );
    return min(1.0, max(i, lightAmbient));
}


// Raytracing
void raytrace(vector<float> const &ray, vector<float> const &origin)
{
    gIntersect = false;
    gDist = 999999.9;

    for (int t = 0; t < nrTypes; t++)
        for (int i = 0; i < nrObjects[t]; i++)
            rayObject(t,i,ray,origin);
}

vector<float> reflect(vector<float>& ray, vector<float>& fromPoint){
    vector<float> N = surfaceNormal(gType, gIndex, gPoint, fromPoint);
    return normalize3(sub3(ray, mul3c(N,(2 * dot3(ray,N)))));
}


//Photon Mapping
vector<float> gatherPhotons(vector<float>& p, int type, int id){
    vector<float> energy = {0.0,0.0,0.0};
    vector<float> N = surfaceNormal(type, id, p, gOrigin);
    for (int i = 0; i < numPhotons[type][id]; i++){
        if (gatedSqDist3(p,photons[type][id][i][0],sqRadius)){
            float weight = max(0.0f, -dot3(N, photons[type][id][i][1] ));
            weight *= (1.0 - sqrt(gSqDist)) / exposure;
            energy = add3(energy, mul3c(photons[type][id][i][2], weight));
        }}
    return energy;
}

void storePhoton(int type, int id, vector<float>& location, vector<float>& direction, vector<float>& energy){
    photons[type][id][numPhotons[type][id]][0] = location;
    photons[type][id][numPhotons[type][id]][1] = direction;
    photons[type][id][numPhotons[type][id]][2] = energy;
    numPhotons[type][id]++;
}

void shadowPhoton(vector<float>& ray){
    vector<float> shadow = {-0.25,-0.25,-0.25};
    vector<float> tPoint = gPoint;
    int tType = gType, tIndex = gIndex;
    vector<float> bumpedPoint = add3(gPoint,mul3c(ray,0.00001));
    raytrace(ray, bumpedPoint);
    vector<float> shadowPoint = add3( mul3c(ray,gDist), bumpedPoint);
    storePhoton(gType, gIndex, shadowPoint, ray, shadow);
    gPoint = tPoint; gType = tType; gIndex = tIndex;
}


vector<float> computePixelColor(float x, float y){
    vector<float> rgb = {0.0,0.0,0.0};
    vector<float> ray = {  (float)(x/resolution - 0.5) ,
        -(float)(y/resolution - 0.5), 1.0f};
    raytrace(ray, gOrigin);

    if (gIntersect){
        gPoint = mul3c(ray,gDist);

        if (gType == 0 && gIndex != 0){
            ray = reflect(ray,gOrigin);
            raytrace(ray, gPoint);
            if (gIntersect){ gPoint = add3( mul3c(ray,gDist), gPoint); }}

        rgb = gatherPhotons(gPoint,gType,gIndex);
    }
    return rgb;
}

void drawPhoton(vector<float> const &rgb, vector<float> const &p, bool mode){
    if (mode && p[2] > 0.0){

        int x = (resolution/2) + (int)(resolution *  p[0]/p[2]);
        int y = (resolution/2) + (int)(resolution * -p[1]/p[2]);
        if (y <= resolution) {
            strokeColor = color(255.0*rgb[0],255.0*rgb[1],255.0*rgb[2]);
//            glutSetWindow(subWin2);
            drawPoint(x,y);
//            glutSetWindow(subWin2);
//            drawPoint(x,y);
        }

    }
}


void emitPhotons(){


    emitDone = false;
    for (int t = 0; t < nrTypes; t++)
        for (int i = 0; i < nrObjects[t]; i++)
            numPhotons[t][i] = 0;

      for (int i = 0; i < nrPhotons; i++){
        int bounces = 1;
        vector<float> rgb = {1.0,1.0,1.0};
        vector<float> ray = normalize3( rand3(1.0) );
        vector<float> prevPoint = Light;


        while (prevPoint[1] >= Light[1]){ prevPoint = add3(Light, mul3c(normalize3(rand3(1.0)), 0.75));}
        if (abs(prevPoint[0]) > 1.5 || abs(prevPoint[1]) > 1.2 ||
            gatedSqDist3(prevPoint,spheres[0],spheres[0][3]*spheres[0][3])) bounces = nrBounces+1;

        raytrace(ray, prevPoint);    //Trace the Photon's Path

        while (gIntersect && bounces <= nrBounces){
            gPoint = add3( mul3c(ray,gDist), prevPoint);
            rgb = mul3c (rgb, 1.0/sqrt(bounces));
            storePhoton(gType, gIndex, gPoint, ray, rgb);
            drawPhoton(rgb, gPoint, mode3D);
            shadowPhoton(ray);
            ray = reflect(ray,prevPoint); //Bounce the Photon
            raytrace(ray, gPoint);        //Trace It to Next Location
            prevPoint = gPoint;
            bounces++;
        }
    }

    emitDone = true;
}

void display1 () {
    glutSetWindow(subWin1);
    display_scene();
}

void display2 () {
    glutSetWindow(subWin2);
    display_photon();
}

void renderAll() {
//        printf("all ");
    display1();
    display2();
    glutPostRedisplay();
}

void reset() {
    pRow=0; pCol=0; pIteration=1; pMax=2;
    empty=true;
    if (!mode3D) emitPhotons();
    glutPostRedisplay();
}

void render(){
    int x,y,iterations = 0;
    vector<float> rgb = {0.0,0.0,0.0};

    while (iterations < (mouseDragging ? 1024 : max(pMax, 256) )){
//        printf("%d %d %d %d\n", pIteration, pRow, pCol, pMax);
        if (pCol >= pMax) {
            pRow++;
            pCol = 0;
            if (pRow >= pMax) {
                pIteration++;
                pRow = 0;
                pMax = int(pow(2, pIteration));
            }
        }

        bool pNeedsDrawing = (pIteration == 1 || odd(pRow) || (!odd(pRow) && odd(pCol)));
        x = pCol * (resolution/pMax); y = pRow * (resolution/pMax);
        pCol++;

        if (pNeedsDrawing){
            iterations++;
            rgb = mul3c( computePixelColor(x,y), 255.0);
            strokeColor = color(rgb[0],rgb[1],rgb[2]); fillColor = color(rgb[0],rgb[1],rgb[2]);
            rect(x,y,(resolution/pMax)-1,(resolution/pMax)-1);
        }
    }

    if (pRow == resolution-1) {empty = false;}
}

void mousePressed(){
    sphereIndex = 3;
    vector<float> mouse3 = {(mouseX - resolution/2)/s, -(mouseY - resolution/2)/s, 0.5f*(spheres[0][2] + spheres[1][2])};
    if (gatedSqDist3(mouse3,spheres[0],spheres[0][3])) sphereIndex = 0;
    else if (gatedSqDist3(mouse3,spheres[1],spheres[1][3])) sphereIndex = 1;
    else if (gatedSqDist3(mouse3,spheres[2],spheres[2][3])) sphereIndex = 2;

    currentWindow = glutGetWindow();
}

void mouseDragged(){
    if (prevMouseX > -9999 && sphereIndex > -1){
        if (sphereIndex < nrObjects[0]){ //Drag Sphere
            spheres[sphereIndex][0] += (mouseX - prevMouseX)/s;
            spheres[sphereIndex][1] -= (mouseY - prevMouseY)/s;
        } else { //Drag Light
            Light[0] += (mouseX - prevMouseX)/s; Light[0] = constrain(Light[0],-1.4f,1.4f);
            Light[1] -= (mouseY - prevMouseY)/s; Light[1] = constrain(Light[1],-0.4f,1.2f);
        }
        reset();
    }

    prevMouseX = mouseX; prevMouseY = mouseY; mouseDragging = true;
}

static void init () {
    glEnable(GL_DEPTH_TEST);

    glDepthFunc(GL_LEQUAL);

    glPolygonOffset (1., -1.);

    glEnable(GL_RESCALE_NORMAL);

    float ambient [] = {0, 0, 0, 1};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

    glFrontFace(GL_CW);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glColorMaterial(GL_FRONT_AND_BACK, GL_SPECULAR);
}


void display () {
    glutSetWindow(mainWin);
//    glutSwapBuffers();
}


void display_scene () {

    glLoadIdentity();

    gluLookAt (width/4.0, height/2.0, (height/2.0) / tan(M_PI*60.0 / 360.0),
               width/4.0, height/2.0, 0, 0, 1, 0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double cameraZ = height/2.0 / tan(M_PI*60/360);
    gluPerspective(60, width/2.0/height, cameraZ/100.0, cameraZ*10.0);


    glScalef(1,-1,1);

    glMatrixMode(GL_MODELVIEW);



            render();


    glutSwapBuffers() ;
}


void display_photon () {
    glLoadIdentity();

    gluLookAt (width/4.0, height/2.0, (height/2.0) / tan(M_PI*60.0 / 360.0),
               width/4.0, height/2.0, 0, 0, 1, 0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double cameraZ = height/2.0 / tan(M_PI*60/360);
    gluPerspective(60, width/2.0/height, cameraZ/100.0, cameraZ*10.0);

    glScalef(1,-1,1);

    glMatrixMode(GL_MODELVIEW);

   if (empty){
       strokeColor = 0;
       fillColor = 0;
       rect(0,0,resolution-1,resolution-1);
       emitPhotons();
       empty = false;
   }

    int a = 0,b = 0;

        for (int i = 0; i < 3; i++) {
            a += numPhotons[0][i];
            for (int j = 0; j < numPhotons[0][i]; j++) {
                drawPhoton(photons[0][i][j][0], photons[0][i][j][2], true);
            }

        }

        for (int i = 0; i < 5; i++) {
            b += numPhotons[1][i];
            for (int j = 0; j < numPhotons[1][i]; j++) {
                drawPhoton(photons[0][i][j][0], photons[0][i][j][2], true);
            }

        }

//    printf("%d %d\n", a,b);

    glutSwapBuffers() ;
}

static void reshape (int wid, int hgt)
{
    printf("%d, %d\n", wid, hgt);

    width = wid;
    height = hgt;

    background (0);
    glutSetWindow(mainWin);
    init();

    glutSetWindow(subWin1);
    glutPositionWindow(0, 0);
    glutReshapeWindow(wid/2, hgt);
    glViewport(0,0,wid/2,hgt);


    glutSetWindow(subWin2);
    glutPositionWindow(wid/2, 0);
    glutReshapeWindow(wid/2, hgt);
    glViewport(0,0,wid/2,hgt);

}


void refresh (int) {
    renderAll();
    glutTimerFunc (1000/frameRate, refresh, 0);

}


void mousemotion (int x, int y) {
    mouseX = x;
    mouseY = y;
    if (bMousePressed) {
        mouseDragged();
    }
}


static void mouse (int button, int state, int x, int y) {
    mouseX = x;
    mouseY = y;

    bMousePressed = state == GLUT_DOWN;
    if (bMousePressed) {
        mousePressed();
    }
    else {
        mouseReleased();
    }
}


static void keyboard (unsigned char ch, int x, int y) {
    bKeyPressed = true;
    key = ch;
    keyCode = ch;

    if (ch=='1') {
        mode3D = false;
        reset();
    } else if (ch=='2') {

        mode3D = true;
        reset();

    } else if (ch=='w') {
        sqRadius += 0.1;
        printf("%f\n", sqRadius);
        reset();
    } else if (ch=='s') {
        sqRadius -= 0.1;
        if (sqRadius < 0.1) sqRadius = 0.1;
        printf("%f\n", sqRadius);
        reset();
    }
}



static void special (int ch, int x, int y) {
    bKeyPressed = true;
    keyCode = 0;
    switch (ch) {
        case GLUT_KEY_LEFT:
            nrBounces -= 1;
            if (nrBounces < 1) nrBounces = 1;
            printf("Bounce time limit: %d\n", nrBounces);
            reset();
            break;
        case GLUT_KEY_UP:
            nrPhotons += 100;
            printf("Photon amount: %d\n", nrPhotons);
            reset();
            break;
        case GLUT_KEY_RIGHT:
            nrBounces += 1;
            printf("Bounce time limit: %d\n", nrBounces);
            reset();
            break;
        case GLUT_KEY_DOWN:
            nrPhotons -= 100;
            if (nrPhotons < 0) nrPhotons = 0;
            printf("Photon amount: %d\n", nrPhotons);
            reset();
            break;
    }

}

int main (int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
//    glutTimerFunc (1000/frameRate, refresh, 0);


    glutInitWindowPosition(50, 100);
    glutInitWindowSize (resolution*2, resolution);
    width = resolution*2;
    height = resolution;
    mainWin = glutCreateWindow ("Photon Mapping");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(renderAll);
    glutMotionFunc (mousemotion);
    glutMouseFunc (mouse);
    glutKeyboardFunc (keyboard);
    glutSpecialFunc(special);
//
//    glutInitWindowPosition(resolution + 50, 100);
//    glutInitWindowSize (resolution, resolution);
//    subWin2 = glutCreateWindow("Photon Visualization");
//    glutDisplayFunc(display_photon);
//    glutMotionFunc (mousemotion);
//    glutMouseFunc (mouse);
//    glutKeyboardFunc (keyboard);
//    glutSpecialFunc(special);

        subWin1 = glutCreateSubWindow(mainWin, 0, 0, resolution, resolution);
        glutDisplayFunc(display1);
        glutMotionFunc (mousemotion);
        glutMouseFunc (mouse);
        glutKeyboardFunc (keyboard);
        glutSpecialFunc(special);


        subWin2 = glutCreateSubWindow(mainWin, resolution, 0, resolution, resolution);
        glutDisplayFunc(display2);
        glutMotionFunc (mousemotion);
        glutMouseFunc (mouse);
        glutKeyboardFunc (keyboard);
        glutSpecialFunc(special);

    glutSetWindow(subWin1);
    emitPhotons();
    reset();
    glutMainLoop();
}
