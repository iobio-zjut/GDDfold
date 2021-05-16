
#ifndef INCLUDED_protocols_abinitio_LJAngleRotation_hh
#define INCLUDED_protocols_abinitio_LJAngleRotation_hh

#include<vector>
#include<map>
#include<string>
#include <core/types.hh>
#include <math.h>
#include <numeric/xyzVector.hh>

using namespace std;
using namespace core;


namespace protocols {
namespace abinitio {
  
class LJAngleRotation {
public:
	LJAngleRotation() {};
	
	void init();
	
	void get_parameters(map<string,Real> &parametersMap);
		
	Real score(vector<Real> &rotationAngles);
	
	Real difference(vector<Real> &Angles_1, vector<Real> &Angles_2);
	
	vector<vector<Real> > angle_learning();
	
	Size dimension(){
		return rotationAxis.size();
	};
	
	void rotation_parameters(vector<vector<Real> > &rotation_axis, vector<numeric::xyzVector<Real> > &fixed_points, vector<numeric::xyzVector<Real> > &rotation_points, vector<pair<pair<Size, Size>, Real > > &L_contact_index_conf);
  
private:

	vector<numeric::xyzVector<Real> > rotationPoints;
	vector<numeric::xyzVector<Real> > fixedPoints;
	vector<pair<pair<Size, Size>, Real > > contact_index_conf;
//	vector<pair<pair<Size, Size>, vector<Real> > > restraint_distance_map;
	
	vector<vector<Real> > rotationAxis;
	
private:
	Size NP2;
	Size G2;
	Real F;
	Real CR;
	Real KTl;
	Real KT_reciprocal;
//	Real Driver_Angle;
	Real Max_Disturbance;
	Real Use_best;
	Size N_Candidate;
	Real Greedy_strategy;
//	Size Convergence_G;
	
public:

	static Real TrialScore;
	static bool isSimilar(Real &score){
		return ( fabs(score - TrialScore) < 0.00001 );
	}
};

} //abinitio
} //protocols

#endif

