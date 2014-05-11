/* JSB: Separated header file from main file of Jensen, 1991 PhotonMap class:
   Implementation in PhotonMap.cpp */ 

#ifndef PHOTONMAP_H
#define PHOTONMAP_H
#include <vector>

using namespace std;

/* This is the photon
 * The power is not compressed so the
 * size is 28 bytes
*/
//**********************
typedef struct Photon {
//**********************
  float pos[3];                  // photon position
  short plane;                   // splitting plane for kd-tree
  unsigned char theta, phi;      // incoming direction
  float power[3];                // photon power (uncompressed)
} Photon;


/* This structure is used only to locate the
 * nearest photons
*/
//**********************
typedef struct NearestPhotons {
//**********************
    int max; 
    int found; 
    int got_heap; 
    float pos[3]; 
    float *dist2; 
    const Photon **index; 
} NearestPhotons;


/* This is the Photon_map class
 */
//****************
class PhotonMap {
//****************

 public:
    PhotonMap(int max_phot);
    ~PhotonMap();

    Photon* getPhotons(void);
    int get_stored_photons(void);
  
    void store(
      const vector<float>& power,        // photon power
      const vector<float>& pos,          // photon position
      const vector<float>& dir );        // photon direction

    void scale_photon_power(
      const float scale);          // 1/(number of emitted photons)

    void balance(void);            // balance the kd-tree (before use!)

    void irradiance_estimate(
     vector<float >& irrad,              // returned irradiance
      const vector<float>& pos,          // surface position
      const vector<float >& normal,       // surface normal at pos
      const float max_dist,        // max distance to look for photons
      const int nphotons ) const;  // number of photons to use

    void locate_photons(
      NearestPhotons *const np,    // np is used to locate the photons
      const int index) const;      // call with index = 1

    void photon_dir(
      float *dir,                  // direction of photon (returned)
      const Photon *p) const;      // the photon

 private:
    void balance_segment(
      Photon **pbal, 
      Photon **porg, 
      const int index,
      const int start, 
      const int end );

    void median_split(
      Photon **p, 
      const int start, 
      const int end,
      const int median, 
      const int axis );
  
    Photon *photons; 

    int stored_photons; 
    int half_stored_photons; 
    int max_photons; 
    int prev_scale; 

    float costheta[256]; 
    float sintheta[256]; 
    float cosphi[256]; 
    float sinphi[256]; 
  
    float bbox_min[3];     // use bbox_min
    float bbox_max[3];     // use bbox_max
};

#endif
