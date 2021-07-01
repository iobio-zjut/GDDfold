// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/ClassicAbinitio.hh
/// @brief header file for ClassicAbinitio protocol
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka


#ifndef INCLUDED_protocols_abinitio_ClassicAbinitio_hh
#define INCLUDED_protocols_abinitio_ClassicAbinitio_hh

// Unit Headers


#include <protocols/simple_moves/FragmentMover.fwd.hh>

#include <core/types.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/abinitio/Protocol.hh>
#include <protocols/moves/TrialMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <vector>

#include <core/fragment/FragSet.fwd.hh>
#include <utility/vector1.hh>


#include <protocols/simple_moves/BBGaussianMover.hh>

#include <core/pose/Pose.hh>

#include <math.h>
#include<fstream>
#include<map>
//#include "spline.h"
///========================================================================
#include <protocols/abinitio/LJAngleRotation.fwd.hh>
#include <protocols/abinitio/LJAngleRotation.hh>
using namespace std;

extern core::pose::Pose nativePose;
extern core::pose::Pose initPose;
extern ofstream Local_acc;

namespace protocols {
namespace abinitio {

/// Move these forward declarations to ClassicAbinitio.fwd.hh
class ClassicAbinitio;
typedef utility::pointer::shared_ptr< ClassicAbinitio > ClassicAbinitioOP;

//@ brief The Classic Abinitio protocol from rosetta++
/*!
@ detail
general usage:
ClassicAbinitio  abinitio;
abinitio.init( pose );
...
while(nstruct) {
abinitio.apply( pose );
}

call ClassicAbinitio::register_options() before core::init::init to add relevant options to the applications help

, with the following
stages, all of which uses a different ScoreFunction based on the cen_std.wts in minirosetta_database:

- Stage 1: large (usually 9mer) randomly selected fragment insertions, only VDW term turned on.
Uses score0.wts_patch and runs for either a maximum of 2000 cycles or until all moveable phi/psi values
have been changed.

- Stage 2: large randomly selected fragment insertions, more score terms turned on. Uses score1.wts_patch
and runs for 2000 cycles.

- Stage 3: uses large randomly selected fragment insertions, although the size of the fragment insertions
is tunable via the set_apply_large_frags( bool ) method. Alternates between score2.wts_patch and score5.wts_patch,
running tunable numbers of 2000-cycle iterations between the two scoring functions.

- Stage 4: uses small (usually 3mer) fragment insertions with the fragment selection based on the Gunn cost for
finding local fragment moves. Runs for 4000-cycles and uses score3.wts_patch.

The class implements the basic abinito approach as known from rosetta++. We tried to set this up, such that
behaviour of the protocol can be changed in many different ways ( see, e.g., FoldConstraints ). To be able to change the
behaviour of the protocol easily the class-apply function and methods called therein (e.g., prepare_XXX() / do_XXX_cycles() ) should
not directly change moves or trials. A reference to the currently used score-function should be obtained by
mc().score_function() ...

Behaviour can be changed in the following ways:

use non-classic FragmentMover  --> eg. not uniformly sampled fragments, but using some weighting
--> large and small moves doesn't have to be 3mers and 9mers... use other movers...
---> or other fragets for the "convenience constructor"
use custom trial classes --> overload update_moves()

change sampling behaviour:
overload prepare_XXX() methods: these are called before the cycling for a certain stage begins
overload do_stageX_cycles() : the actual loops over trial-moves ...

change scoring functions:
overload set_default_scores()
weight-changes effective for all stages: set_score_weight()

*/


class ClassicAbinitio : public Protocol {
	typedef Protocol Parent;
public:
	enum StageID {
		ALL_STAGES = 0,
		STAGE_1,
		STAGE_2,
		STAGE_3a,
		STAGE_3b,
		STAGE_4,
		STAGE_4rot,
		STAGE_5
	};
	/// @brief This constructor does not work -- Fix it before using it.
	// constructor: supply mover classes for Fragment Moves
	ClassicAbinitio(
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small,
		int  /*dummy otherwise the two constructors are ambiguous */
	);

	/// @brief constructor: supply fragsets for large and small fragment moves
	ClassicAbinitio(
		core::fragment::FragSetCOP fragset_small,
		core::fragment::FragSetCOP fragset_large,
		core::kinematics::MoveMapCOP movemap
	);

	/// @brief Explicit copy constructor since this class contains OPs of other classes
	ClassicAbinitio( ClassicAbinitio const & src );

	/// @brief Explicit destructor since this class contains OPs of other classes
	~ClassicAbinitio() override;

	//@brief setup moves, mc-object, scores
	//@details can't call this from constructor; virtual functions don't operate
	//until construction has completed.

	void init( core::pose::Pose const& pose ) override;

	//@brief ClassicAbinitio has virtual functions... use this to obtain a new instance

	moves::MoverOP clone() const override;

	//@brief run protocol on pose
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	//@brief return FramgentMover for smooth_small fragment insertions (i.e., stage4 moves)
	simple_moves::FragmentMoverOP smooth_move_small();

	//@brief return FragmentMover for small fragment insertions ( i.e., stage3/4 moves )
	simple_moves::FragmentMoverOP brute_move_small();

	//@brief return FragmentMover for large fragment insertions (i.e., stage1/2 moves )
	simple_moves::FragmentMoverOP brute_move_large();

	//@brief change the movemap ( is propagated to mover-objects )
	//@detail overload if your extension stores additional moves as member variables
	virtual void set_movemap ( core::kinematics::MoveMapCOP mm );

	//@brief returns current movemap
	core::kinematics::MoveMapCOP movemap();

	//@brief set new instances of FragmentMovers
	void set_moves (
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small
	);

	//@brief set new monte-carlo object
	void set_mc( moves::MonteCarloOP );

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarlo & mc() {
		return *mc_;
	}

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarlo const & mc() const {
		return *mc_;
	}

	//@brief override cycle setting ( init() -> sets it according to cmd-line options )
	virtual void set_cycles( core::Real increase_cycles = 1.0 );

	//@brief
	Size total_trials() const {
		return total_trials_;
	}

	/// @brief for debugging, one wants to have access to the native pose.
	// void set_native_pose( core::pose::Pose const & pose );

	//@brief set weight - effective for all scoring functions  stage == -1 --> set weight for all stages
	// mod -1 (ignored) ,  in stage3 mod = 1 --> 3a , mod = 2 --> 3b
	virtual void set_score_weight( core::scoring::ScoreType, core::Real setting, StageID stage = ALL_STAGES );

protected:
	//@brief called to notify about changes regarding movers... new movemap / new instances of FragmentMover
	virtual void update_moves();

	//@brief called by init() --- calls all set_default_XXX methods
	virtual void set_defaults( core::pose::Pose const& pose );

	//@brief read out cmd-line options
	virtual void set_default_options();

	//@brief register cmd-line options in option system ( call before core::init::init )
public:
	static void register_options();
protected:
	//@brief construct default monto-carlo object
	virtual void set_default_mc(
		core::pose::Pose const& pose,
		core::scoring::ScoreFunction const& scorefxn
	);

	//@brief assigns default score-functions to the 4 stages--> score0 / score1 / (score2/score5) / score3
	virtual void set_default_scores();

	//@brief currently used score function ( depends on stage )
	core::scoring::ScoreFunction const& current_scorefxn() const;

	//@brief set current scorefunction
	void current_scorefxn( core::scoring::ScoreFunction const& scorefxn );

	//@brief If appropriate for the stage then recover low mc.
	void recover_low( core::pose::Pose& pose, StageID stage );

	//@brief this is called if the scorefxn is replaced with a new one
	virtual void replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real intra_stage_progress );

	//@brief set individual weight of current scorefunction
	//@details NOTE: does not change the predefined weights, this information is lost each time the score is set
	//  weight will be overwritten with default ->at beginning of stages, each iteration of loop in stage3
	// we could change this behaviour by having a pointer to the current score and making sure that nobody can
	// change the score inside the mc-object (only const accessor )
	void set_current_weight( core::scoring::ScoreType type, core::Real setting );

	//@brief run cycles for different scoring_stages, return number of steps used
	virtual bool do_stage1_cycles( core::pose::Pose &pose );
	virtual bool do_stage2_cycles( core::pose::Pose &pose );
	virtual bool do_stage3_cycles( core::pose::Pose &pose );
	virtual bool do_stage4_cycles( core::pose::Pose &pose );
	virtual bool do_stage5_cycles( core::pose::Pose &pose );//vats

	//@brief returns true if pose is < 3.0 A rms to last pose sent to this function
	// bool convergence_check( core::pose::Pose const & pose );

	//@brief returns the Mover that is applied inside the stage1 loop
	virtual moves::TrialMoverOP
	stage1_mover( core::pose::Pose &pose,  moves::TrialMoverOP trials_in );


	//@brief returns the Mover that is applied inside the stage3 double loop
	virtual moves::TrialMoverOP
	stage2_mover( core::pose::Pose &pose,  moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage3 double loop
	virtual moves::TrialMoverOP
	stage3_mover( core::pose::Pose & pose, int lct1, int lct2,  moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage4 loop
	virtual moves::TrialMoverOP
	stage4_mover( core::pose::Pose &pose, int kk, moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage4 loop
	virtual moves::TrialMoverOP
	stage4rot_mover( core::pose::Pose &pose, int kk, moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage5 loop
	virtual moves::TrialMoverOP
	stage5_mover( core::pose::Pose &pose, moves::TrialMoverOP trials_in );//vats

	//@brief called by update_moves() creates the instances of TrialMover with the FragmentMoves
	virtual void set_trials();

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_large();

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_small();

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_smooth();

	// anything you want to have done before the stages ?
	//@brief prepare_stageX is called before do_stageX_cycles... overload to change status/scoring/conformation....
	virtual bool prepare_stage1( core::pose::Pose &pose );
	virtual bool prepare_stage2( core::pose::Pose &pose );
	virtual bool prepare_stage3( core::pose::Pose &pose );
	virtual bool prepare_stage4( core::pose::Pose &pose );
	virtual bool prepare_stage5( core::pose::Pose &pose );//vats

	//@brief called in each iteration of inner loop in stage3 before stage3_cycles_ of trials commence
	virtual bool prepare_loop_in_stage3(
		core::pose::Pose&,
		Size, /* loop_iteration*/
		Size  /* total_iterations */
	);

	//@brief called in each iteration of the loop in stage4 before the stage4_cycles_ of trials commence
	virtual bool prepare_loop_in_stage4(
		core::pose::Pose&,
		Size, /* loop_iteration*/
		Size  /* total_iterations */
	);


public:
	//@brief accessor for boolean flag: just_smooth_cycles
	inline
	bool just_smooth_cycles() const {
		return just_smooth_cycles_;
	}

	//@brief Accessor for number of stage1 cycles
	inline
	Size stage1_cycles() const {
		return stage1_cycles_;
	}

	//@brief Accessor for number of stage2
	inline
	Size stage2_cycles() const {
		return stage2_cycles_;
	}

	//@brief Accessor for number of stage3 cycles
	inline
	Size stage3_cycles() const {
		return stage3_cycles_;
	}

	//@brief Setter for number of stage4 cycles
	inline
	void set_stage4_cycles(Size stage4_cycles_new) {
		stage4_cycles_ = stage4_cycles_new;
	}

	//@brief Accessor for number of stage4 cycles
	inline
	Size stage4_cycles() const {
		return stage4_cycles_;
	}

	//@brief Accessor for number of stage5 cycles  //vats
	inline
	Size stage5_cycles() const {
		return stage5_cycles_;
	}


	//@brief query this flag if you have time-intensive stuff and cut it short
	bool bQuickTest() const {
		return bQuickTest_;
	}

	void set_skip_stage1 ( bool setting ) {
		bSkipStage1_ = setting;
	}

	void set_skip_stage2 ( bool setting ) {
		bSkipStage2_ = setting;
	}

protected:
	void output_debug_structure( core::pose::Pose & pose, std::string prefix ) override;

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarloOP mc_ptr();

protected:
	//@brief  cycles_numbers  -- modified by set_cycles(), set_default_options()
	Size stage1_cycles_; //score0
	Size stage2_cycles_; //score1
	Size stage3_cycles_; //score25
	Size stage4_cycles_; //score3
	Size stage5_cycles_; //score3

	core::Real stage4_cycles_pack_rate_;

public:
		Size insert_po,ind;   //插入的位置
             //选择片段
 
 		vector<vector<Real> > d2;
 		 vector<Real> resnum;
 		 int findElement(vector<Real> v, Real key);

private:
	//@brief scoring functions -- modified by set_default_scores() / set_score_weight()
	core::scoring::ScoreFunctionOP score_stage1_;  //score0
	core::scoring::ScoreFunctionOP score_stage2_;  //score1
	core::scoring::ScoreFunctionOP score_stage3a_; //score2
	core::scoring::ScoreFunctionOP score_stage3b_; //score5
	core::scoring::ScoreFunctionOP score_stage4_;  //score3
	core::scoring::ScoreFunctionOP score_stage4rot_;  //score_cenrot
	core::scoring::ScoreFunctionOP score_stage4rot_sc_;  //score_cenrot
	core::scoring::ScoreFunctionOP score_stage5_; //score3 //vats lets try score3 first

	//@brief flags
	bool apply_large_frags_;   //above contig_cut2
	bool short_insert_region_; //below contig_cut3
	bool just_smooth_cycles_;
	bool bQuickTest_;  // this flag might land in base-class
	bool close_chbrk_; // closing chainbreaks in old style pose_abinitio

	//@brief a temperature
	core::Real temperature_;

	//@brief movemap --> which dofs can be moved during abinitio
	core::kinematics::MoveMapCOP movemap_;

	//@brief a MonteCarlo object -- set_default_mc() , access: mc()
	moves::MonteCarloOP mc_;

	// Large and small Fragments
	simple_moves::FragmentMoverOP brute_move_small_;
	simple_moves::FragmentMoverOP brute_move_large_;
	simple_moves::FragmentMoverOP smooth_move_small_;

	// TrialMovers
	moves::TrialMoverOP trial_large_;
	moves::TrialMoverOP trial_small_;
	moves::TrialMoverOP smooth_trial_small_;

	moves::TrialMoverOP trial_small_pack_;
	moves::TrialMoverOP smooth_trial_small_pack_;

	// Cenrot Sidechain Mover
	minimization_packing::PackRotamersMoverOP pack_rotamers_;

	Size total_trials_;

public:
	bool bSkipStage1_;
	bool bSkipStage2_;
	bool bSkipStage3_;
	bool bSkipStage4_;
	bool bSkipStage5_;
//用户添加变量

	utility::vector1< StageID > recover_low_stages_;

private:
	/// @brief Private, unimplemented assignment operator to prevent assignment of this class.
	/// Copy-constructor copying only.
	ClassicAbinitio const & operator = ( ClassicAbinitio const & src );

///================================= @CGLFold ================================================================
private:
	///@brief save contact
	struct Contact 
	{
	    Size residue1;
	    Size residue2;
	    Real confidence;
	    
	    
	    Contact() {residue1 = 0; residue2 = 0; confidence = 0;}
	    Contact(Size res1, Size res2, Real con):residue1(res1), residue2(res2), confidence(con) {}
	};
	vector<Contact> Contacts;
	Real Contact_score(core::pose::Pose &pose);
	Real Cscore_full;
	void norm(Real &angel);
	void convex_mutation1(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
	Real Ave_Cscore;
	Real Average_Cscore_of_Populatin(vector<Real> &population_Cscore);
	Real Average_Cscore_of_Populatin(Real &old_Cscore, Real &new_Cscore, Real &average_Cscore);
	void convex_mutation2(vector< core::pose::Pose >& mypose,vector<Real>& population_Cscore,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id,Real rate);
	bool Distance_accept(Real target_Dscore, Real trial_Dscore);
	//LT
	void gassu_mutation(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
	void convex_crossover(core::pose::Pose& target, core::pose::Pose& posemutation, int seq_length);
	void Output_last_generation(vector< core::pose::Pose >& population, vector< Real >& population_Dscore, char c);
	Real Ave_Dscore;
	Real Ave_Rmsd;
	bool evaluation(core::pose::Pose& pre_pose,core::pose::Pose& pro_pose,Real rate);
	Real Average_Dscore_of_Populatin(vector<Real> &population_Dscore);
	Real Average_Dscore_of_Populatin(Real &old_Dscore, Real &new_Dscore, Real &average_Dscore);
	///@brief Distance 
	struct Distance 
	{
	    Size res1;
	    Size res2;
	    Real dist;
		Real pro;
	    Real mean;
	    Real var;
	    
	    Distance() {res1 = 0; res2 = 0; dist = 0; pro = 0; mean = 0; var = 0;}
	    Distance(Size re1, Size re2, Real dis,  Real pr, Real mea, Real va):res1(re1), res2(re2), dist(dis), pro(pr), mean(mea), var(va) {}
	};
	
	vector<pair<pair<Size, Size>, vector<Real> >  > Distance_splines;
	
	
	vector<Distance> Distance_list;
	Real Distance_score(core::pose::Pose &pose);
	Real Local_distance_score(core::pose::Pose &pose,vector<pair<pair<Size, Size>, vector<Real> > > &restraint_dist_map);
	Real Local_distance_score1(core::pose::Pose &pose,vector<pair<pair<Size, Size>, vector<Real> > > &restraint_dist_map);
	void choose_frag(core::pose::Pose &pose);
	//void calculate_distance(core::pose::Pose &pose );

	struct SS_unit 
	{
	    Real L;
	    Real H;
	    Real E;
	    
	    SS_unit() {L = 0; H = 0; E = 0;}
	    SS_unit(Size l, Size h, Real e):L(l), H(h), E(e) {}
	};
	vector<SS_unit> SS_map;
	string SS;
	Real SSscore_full;
	Real SS_score(core::pose::Pose &pose);
private:
	void read_distance();
	void read_contact();
	void read_SS();
	void read_parameters();
	void get_parameters(map<string,Real> &parametersMap);
	
	bool boltzmann_accept( Real const &targetEnergy, Real const &trialEnergy );
	bool boltzmann_accept( Real const &targetEnergy, Real const &trialEnergy, Real const &KT );
	void minmover(pose::Pose & pose);
	
	void Output_last_generation(vector<core::pose::Pose> &population, vector<Real> &population_energy, vector< core::pose::Pose >& targetpose , vector< Real >& targetpose_energy);
//	void Output_pdb_Global(vector<core::pose::Pose> &population);
//	void Output_pdb_Local(vector<core::pose::Pose> &population);
	
	void Output_SPICKER_All_Data(core::pose::Pose &pose, Real const &energy);
	void SPICKER_Demand_All(core::pose::Pose &nativePose);
	
	vector<core::pose::Pose> generate_population(Size const &population_size);
	
	
	vector<core::pose::Pose> generate_population_random(Size const &population_size, const Size& seq_length);
	void learn_period(vector<core::pose::Pose> &population, Size const &population_size, const Size& seq_length, const Size& G_length, scoring::ScoreFunction const & scorefxn);
	vector<core::pose::Pose> k_mediods( const Size& k_number, vector<core::pose::Pose> &population, Size const &population_size, vector<Real> &state_tag);
	vector< vector<Real> > compute_distmatrix(vector<core::pose::Pose> &population1, vector<core::pose::Pose> &population2, Size const &population_size1, Size const &population_size2);
	void entropy_switch_strategy(vector<core::pose::Pose> &population, Size const &population_size, const Size& seq_length, const Size& G_global_, const Size& G_local, scoring::ScoreFunction const & scorefxn, vector<Real> &state_tag, const Size& k_number, vector<core::pose::Pose> &mediods_set);
	void cluster(vector<core::pose::Pose> &population, Size const &population_size, vector<Real> &state_tag, const Size& k_number, vector<core::pose::Pose> &mediods_set);
	void minmover(vector<core::pose::Pose> &population, Size const &population_size, const Size& seq_length, const Size& G_local, scoring::ScoreFunction const & scorefxn);
	//void minmover(core::pose::Pose &Ctrial, scoring::ScoreFunction const & scorefxn);
	void mutation_rand_1(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
	void mutation_current2rand_1(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
	void mutation_best_1(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
	
	core::pose::Pose generate_a_new_pose();
	
	vector<Real> calculate_population_energy_2(vector<core::pose::Pose> &population);
	vector<Real> calculate_population_energy(vector<core::pose::Pose> &population,scoring::ScoreFunction const & scorefxn);
	vector<Real> calculate_population_Cscore(vector<core::pose::Pose> &population);
	vector<Real> calculate_population_Dscore(vector<core::pose::Pose> &population);
	vector<Real> calculate_population_Rmsd(vector<core::pose::Pose> &population);
	vector<Real> calculate_population_SSscore(vector<core::pose::Pose> &population);
	
	
	//void rotation_parameters( vector<numeric::xyzVector<Real> > &fixed_points, vector<numeric::xyzVector<Real> > &rotation_points );
  	
	///@note 旋转点坐标集
	vector<numeric::xyzVector<Real> > rotationPoints;
	///@note 固定点坐标集
	vector<numeric::xyzVector<Real> > fixedPoints;
	
	void SS_search(vector<core::pose::Pose> &population, vector<Real> &population_energy, vector<Real> &population_SS);
	void Global_search(vector<core::pose::Pose> &population, vector<Real> &population_energy, vector<Real> &population_Dscore, vector<Real> &population_Rmsd);
	void Local_search(vector<core::pose::Pose> &population, vector<Real> &population_energy, vector<Real> &population_Dscore, vector<Real> &population_Rmsd, LJAngleRotationOP LJAR);
	
	///@brief set parameters
	Size NP_;
	double NP_reciprocal;
	Size G_LP_;
	Size G_global_;
	Size G_local_;
	Size Knumber;
	Size frag_;
	Real KTc_;
	bool Only_Global_;
//	bool Only_Local_;

	map<string,Real> parametersMap;
	
	vector<Real> PopulationRmsd;
	vector<Real> PopulationEnergy;
	
	
	Size Global_Selection_Mode_;
	Size Local_Selection_Mode_;
	
	ofstream out_put;
	ofstream Ave_acc;
	ofstream Ave_entropy;
	ofstream Global_acc;
//	ofstream Local_acc;
	ofstream dis_map;
	
//	Size Convergence_Length_E_;
//	Size Convergence_Length_C_;
//	Size Convergence_bias_E_;
//	Size Convergence_bias_C_;
	
	///@note Spicker paramters for all accept points 
	ofstream SPICKER_All_Data;
	Size Spicker_All_file_num;
	Size Spicker_All_num;
	///@note Spicker paramters for Global phase accept points
	ofstream SPICKER_Global_Data;
	Size Spicker_Global_file_num;
	Size Spicker_Global_num;
	///@note Spicker paramters for Local phase accept points
	ofstream SPICKER_Local_Data;
	Size Spicker_Local_file_num;
	Size Spicker_Local_num;
	
	Size proteinLength;
	int lenth_of_sequence;
	Size similar_times;
	simple_moves::ClassicFragmentMoverOP FragAssem_;
	simple_moves::ClassicFragmentMoverOP bms, bml, sms;
	
public:	
	///@brief 相似性判断函数，作为泛型算法 find_if 的谓词
	static Real TrialEnergy;
	static bool isSimilar(Real &energy){
		return ( fabs(energy - TrialEnergy) < 0.00001 );
	}
///========= @LJ ============= @LJ ============ @LJ ============= @LJ ============== @LJ ============= @LJ =============== @LJ ==========

};

} // abinitio
} // protocols

#endif
