#include "main.h"
#include <vector>
#include <GLUT/glut.h>
#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;
using namespace pm;


int mouseX = 0;
int mouseY = 0;
bool bMousePressed = false;
int mouseButton = LEFT;
bool bKeyPressed = false;
unsigned char key = 255;
int keyCode;

int width;
int height;


int frameRate = 60;

int initialized = false;

color strokeColor (0,0,0);
color fillColor   (255,255,255);

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

void pm::stroke (const color& c) {
    strokeColor = c;
}


void pm::fill (const color& c) {
    fillColor = c;
}

static unsigned globColorMode = RGB;

static double maxColor = 255;

static void hsb_to_rgb (double h, double s, double v,
                        double& r, double& g, double& b)
{
    double tmp = h*5.9999;
    int hi = int (tmp);
    double f = tmp-hi;
    double p = v * (1-s);
    double q = v * (1-f*s);
    double t = v * (1-(1-f)*s);
    if (hi==0) {
        r = v; g = t; b = p;
    } else if (hi==1) {
        r = q; g = v; b = p;
    } else if (hi==2) {
        r = p; g = v; b = t;
    } else if (hi == 3) {
        r = p; g = q; b = v;
    } else if (hi == 4) {
        r = t; g = p; b = v;
    } else {
        r = v; g = p; b = q;
    }
}

inline unsigned char clamp(double v) {
    return v > 255 ? 255 : v < 0 ? 0 : (unsigned char) v;
}

pm::color::color(double val1, double val2, double val3, double valA) {
    //SCale the values to a range of 255
    val1 = val1/maxColor;
    val2 = val2/maxColor;
    val3 = val3/maxColor;
    if (valA == MAXCOLOR) valA = 1.0;
    else valA = valA/maxColor;

    if (globColorMode != RGB) {
        hsb_to_rgb(val1, val2, val3, val1, val2, val3);
    }

    rgba[0] = clamp(val1*255);
    rgba[1] = clamp(val2*255);
    rgba[2] = clamp(val3*255);
    rgba[3] = clamp(valA*255);
};


color::color(double gray, double alpha){
    if (alpha == MAXCOLOR) { alpha = maxColor;}
    unsigned char val = clamp(gray/maxColor*255);
    rgba[0] = val;
    rgba[1] = val;
    rgba[2] = val;
    rgba[3] = clamp(alpha/maxColor*255);
}

static unsigned rectMode = CORNER;
static std::vector<PVector> ellipseVtx;

static std::vector<PVector> sphereVtx;
static std::vector<unsigned> sphereIdx;


void pm::quad (double x0, double y0,
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
void pm::point (double x, double y, double z)
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
    // Make changes to arguments to reflect the current rectMode
    switch (::rectMode) {
		case CORNER:
			quad (x, y, x+a, y, x+a, y+b, x, y+b);
			break;
		case CENTER:
			quad (x-a/2, y-b/2, x+a/2, y-b/2, x+a/2, y+b/2, x-a/2, y+b/2);
            break;
        case RADIUS:
			quad (x-a, y-b, x+a, y-b, x+a, y+b, x-a, y+b);
		   	break;
		case CORNERS:
			quad (x, y, a, y, a, b, x, b);
		   	break;
    }
}

void pm::background (const color& c) {
    glClearColor (c.rgba[0] * (1.0/255),
                  c.rgba[1] * (1.0/255),
                  c.rgba[2] * (1.0/255),
                  c.rgba[3] * (1.0/255));
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}



int resolution = 512;
int nrTypes = 2;
vector<int> nrObjects = {3,5};
float gAmbient = 0.1;
vector<float> gOrigin = {0.0,0.0,0.0};
vector<float> Light = {0.0,1.2,3.75};
vector< vector<float> > spheres = {{1.0,0.0,4.0,0.5}, {-0.6,-1.0,4.5,0.5}, {0.0,1.0,4.5,0.3}};
vector<vector<float>> planes  = {{0, 1.5},{1, -1.5},{0, -1.5},{1, 1.5},{2,5.0}};


int nrPhotons = 400;
int nrBounces = 3;
bool lightPhotons = true;
float sqRadius = 0.7;
float exposure = 30.0;
vector<vector<int> > numPhotons = {{0,0},{0,0,0,0,0}};

vector<vector<vector<vector<vector<float> > > > > photons(2,
                                                          vector<vector<vector<vector<float> > > >(5, vector<vector<vector<float> > >(5000,
                                                                                                                                      vector<vector<float> >(3, vector<float>(5)))));

bool gIntersect = false;
int gType;
int gIndex;
float gSqDist, gDist = -1.0;
vector<float> gPoint = {0.0, 0.0, 0.0};

bool empty = true;
bool view3D = true;
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

void raySphere(int idx, vector<float> const &r, vector<float> const &o)
{
    vector<float> s = sub3(spheres[idx],o);
    float radius = spheres[idx][3];


    float A = dot3(r,r);
    float B = -2.0 * dot3(s,r);
    float C = dot3(s,s) - sq(radius);
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
    return min(1.0, pm::max(i, lightAmbient));
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
            float weight = pm::max(0.0f, -dot3(N, photons[type][id][i][1] ));
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

vector<float> filterColor(vector<float>& rgbIn, float r, float g, float b){
    vector<float> rgbOut = {r,g,b};
    for (int c=0; c<3; c++) rgbOut[c] = min(rgbOut[c],rgbIn[c]);
    return rgbOut;
}

vector<float> getColor(vector<float>& rgbIn, int type, int index){
    if      (type == 1 && index == 0) { return filterColor(rgbIn, 0.0, 1.0, 0.0);}
    else if (type == 1 && index == 1) { return filterColor(rgbIn, 1.0, 0.0, 0.0);}
    else if (type == 1 && index == 2) { return filterColor(rgbIn, 0.0, 0.0, 1.0);}
    else
    { return filterColor(rgbIn, 1.0, 1.0, 1.0);}
}


vector<float> computePixelColor(float x, float y){
    vector<float> rgb = {0.0,0.0,0.0};
    vector<float> ray = {  (float)(x/resolution - 0.5) ,
        -(float)(y/resolution - 0.5), 1.0f};
    raytrace(ray, gOrigin);

    if (gIntersect){
        gPoint = mul3c(ray,gDist);

        if (gType == 0 && gIndex == 1){
            ray = reflect(ray,gOrigin);
            raytrace(ray, gPoint);
            if (gIntersect){ gPoint = add3( mul3c(ray,gDist), gPoint); }}

        if (lightPhotons){
            rgb = gatherPhotons(gPoint,gType,gIndex);}
        else{
            int tType = gType, tIndex = gIndex;
            float i = gAmbient;
            raytrace( sub3(gPoint,Light) , Light);
            if (tType == gType && tIndex == gIndex)
                i = lightObject(gType, gIndex, gPoint, gAmbient);
            rgb[0]=i; rgb[1]=i; rgb[2]=i;
            rgb = getColor(rgb,tType,tIndex);}
    }
    return rgb;
}

void drawPhoton(vector<float> const &rgb, vector<float> const &p){
    if (view3D && p[2] > 0.0){
        int x = (resolution/2) + (int)(resolution *  p[0]/p[2]);
        int y = (resolution/2) + (int)(resolution * -p[1]/p[2]);
        if (y <= resolution) {stroke(255.0*rgb[0],255.0*rgb[1],255.0*rgb[2]); point(x,y);}}
}

void emitPhotons(){
    for (int t = 0; t < nrTypes; t++)
        for (int i = 0; i < nrObjects[t]; i++)
            numPhotons[t][i] = 0;

    for (int i = 0; i < (view3D ? nrPhotons * 3.0 : nrPhotons); i++){
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
            rgb = mul3c (getColor(rgb,gType,gIndex), 1.0/sqrt(bounces));
            storePhoton(gType, gIndex, gPoint, ray, rgb);
            drawPhoton(rgb, gPoint);
            shadowPhoton(ray);
            ray = reflect(ray,prevPoint); //Bounce the Photon
            raytrace(ray, gPoint);        //Trace It to Next Location
            prevPoint = gPoint;
            bounces++;}
    }
}

void resetRender() {
    pRow=0; pCol=0; pIteration=1; pMax=2;
    empty=true; if (lightPhotons && !view3D) emitPhotons();
    glutPostRedisplay();
}

void render(){
    int x,y,iterations = 0;
    vector<float> rgb = {0.0,0.0,0.0};

    while (iterations < (mouseDragging ? 1024 : pm::max(pMax, 512) )){
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
            stroke(rgb[0],rgb[1],rgb[2]); pm::fill(rgb[0],rgb[1],rgb[2]);
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
}

void mouseDragged(){
    if (prevMouseX > -9999 && sphereIndex > -1){
        if (sphereIndex < nrObjects[0]){ //Drag Sphere
            spheres[sphereIndex][0] += (mouseX - prevMouseX)/s;
            spheres[sphereIndex][1] -= (mouseY - prevMouseY)/s;}
        else{ //Drag Light
            Light[0] += (mouseX - prevMouseX)/s; Light[0] = constrain(Light[0],-1.4f,1.4f);
            Light[1] -= (mouseY - prevMouseY)/s; Light[1] = constrain(Light[1],-0.4f,1.2f);}
        resetRender();}
    prevMouseX = mouseX; prevMouseY = mouseY; mouseDragging = true;
}

void setup(){
    size(resolution,resolution);
    emitPhotons();
    resetRender();
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
    glLoadIdentity();

    gluLookAt (width/2.0, height/2.0, (height/2.0) / tan(M_PI*60.0 / 360.0),
               width/2.0, height/2.0, 0, 0, 1, 0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double cameraZ = height/2.0 / tan(M_PI*60/360);
    gluPerspective(60, width*1.0/height, cameraZ/100.0, cameraZ*10.0);


    glScalef(1,-1,1);

    glMatrixMode(GL_MODELVIEW);


    glDisable(GL_LIGHTING);

    if (view3D){
        if (empty){
            stroke(0); fill(0); rect(0,0,resolution-1,resolution-1);
            emitPhotons(); empty = false;
//            render();
        }
    }
    else{
        if (empty) render();
    }


    glutSwapBuffers() ;
}

static void reshape (int wid, int hgt)
{
    glViewport(0,0,wid,hgt);

    width = wid;
    height = hgt;

    background (200);

    init();

}


void refresh (int) {
    glutPostRedisplay();
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

    if (button == GLUT_LEFT_BUTTON) {
        mouseButton = LEFT;
    } else if (button == GLUT_RIGHT_BUTTON) {
        mouseButton = RIGHT;
    } else {
        mouseButton = CENTER;
    }
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
        view3D = false;
        lightPhotons = false;
        resetRender();
    } else if (ch=='2') {
        view3D = false;
        lightPhotons = true;
        resetRender();
    } else if (ch=='3') {
        view3D = true;
        resetRender();
    } else if (ch=='w') {
        exposure += 0.1;
        resetRender();
    } else if (ch=='d') {
        exposure -= 0.1;
        resetRender();
    }
}



static void special (int ch, int x, int y) {
    bKeyPressed = true;
    keyCode = 0;
    switch (ch) {
        case GLUT_KEY_LEFT:
            nrBounces -= 1;
            resetRender();
            break;
        case GLUT_KEY_UP:
            nrPhotons += 100;
            resetRender();
            break;
        case GLUT_KEY_RIGHT:
            nrBounces += 1;
            resetRender();
            break;
        case GLUT_KEY_DOWN:
            nrPhotons -= 100;
            resetRender();
            break;
    }

}


void pm::size (unsigned w, unsigned h, const char* name) {
    if (initialized) {
        glutReshapeWindow (w, h);
        glutSetWindowTitle (name);
    } else {
        glutInitWindowSize (w, h);
        width = w;
        height = h;
        glutCreateWindow (name);
        glutReshapeFunc(reshape);
        glutDisplayFunc(display);
        glutMotionFunc (mousemotion);
        glutMouseFunc (mouse);
        glutKeyboardFunc (keyboard);
        glutSpecialFunc(special);
        initialized = true;
    }
    
}

int main (int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutTimerFunc (1000/frameRate, refresh, 0);
    setup();
    glutMainLoop();
}
