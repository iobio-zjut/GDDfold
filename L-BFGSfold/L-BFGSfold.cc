// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ClassicAbinitio.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka

// Unit Headers
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/simple_moves/SymmetricFragmentMover.hh>

// Package Headers
#include <protocols/simple_moves/GunnCost.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/WhileMover.hh>
#include <protocols/abinitio/AllResiduesChanged.hh>


//for minmover
#include <protocols/jd2/JobDistributor.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/Jump.hh>
#include <utility/excn/Exceptions.hh>



//for cenrot
#include <protocols/moves/CompositionMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
//#include <protocols/simple_moves/BackboneMover.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/ozstream.hh>
#include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#ifdef WIN32
#include <ctime>
#endif

//debug

#include <protocols/moves/MonteCarlo.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

///====================================================headers===================================================
//#include <protocols/abinitio/LJImplementation.fwd.hh>
//#include <protocols/abinitio/LJImplementation.hh>

#include <protocols/abinitio/LJAngleRotation.fwd.hh>
#include <protocols/abinitio/LJAngleRotation.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Ramachandran.hh>

#include <protocols/simple_moves/FragmentMover.fwd.hh>

///=========================== @CGLFold ===============================
///@note 用户添加程序
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<algorithm>
//#include<unordered_map>
#include<set>
#include<math.h>
#include<cstdlib>
///@brief make_pair头文件 pair
#include<utility>

#include <numeric/random/random.hh>
///@brief 产生随机数
#include<stdlib.h>
#include<time.h>
///@brief 计算rmsd的函数所在的头文件
#include <core/scoring/rms_util.hh>
///@brief 控制输出格式
#include<iomanip>
///@brief 读取天然态蛋白质
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
///@brief 输出pdb文件
#include <basic/options/keys/out.OptionKeys.gen.hh>
///@brief 计算二级结构元素的距离
//#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/Atom.hh>
#include <numeric/xyzVector.hh>

#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>

#include <core/kinematics/MoveMap.hh>
#include <utility/pointer/owning_ptr.hh>
#include <memory>

#include<stdlib.h>
//#include "spline.h"
#include <numeric>

///ctreate a new direct
#include<sys/stat.h>
#include<sys/types.h>
#include"Fragment_assembly.hh"
#include<fstream>
#include<sstream>
using namespace std;
//using namespace SplineSpace;

core::pose::Pose initPose;
core::pose::Pose nativePose;
ofstream Local_acc;
///=========================== @CGLFold ===============================

static basic::Tracer tr( "protocols.abinitio" );

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/*!
@detail call this:
ClassicAbinitio::register_options() before devel::init().
Derived classes that overload this function should also call Parent::register_options()
*/

// This method of adding options with macros is a pain in the ass for people
// trying to nest ClassicAbinitio as part of other protocols. If you don't call
// ClassicAbinitio::register_options() in your main function, you get a really
// unintuitive segfault as the options system doesn't know about the options
// listed below. The solution is to call register_options() in your main method
// before devel::init(), which is really ugly as the main method shouldn't need
// to know what protocols are called, and it's prone to error because it's an
// easy thing to forget.
// This should get some more thought before it becomes the standard way to add options.

void protocols::abinitio::ClassicAbinitio::register_options() {
	Parent::register_options();
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( OptionKeys::abinitio::increase_cycles );
	option.add_relevant( OptionKeys::abinitio::smooth_cycles_only );
	option.add_relevant( OptionKeys::abinitio::debug );
	option.add_relevant( OptionKeys::abinitio::skip_convergence_check );
	option.add_relevant( OptionKeys::abinitio::log_frags );
	option.add_relevant( OptionKeys::abinitio::only_stage1 );
	option.add_relevant( OptionKeys::abinitio::end_bias );
	option.add_relevant( OptionKeys::abinitio::symmetry_residue );
	option.add_relevant( OptionKeys::abinitio::vdw_weight_stage1 );
	option.add_relevant( OptionKeys::abinitio::override_vdw_all_stages );
	option.add_relevant( OptionKeys::abinitio::recover_low_in_stages );
	option.add_relevant( OptionKeys::abinitio::close_chbrk );
}


namespace protocols {
namespace abinitio {
  
    Real ClassicAbinitio::TrialEnergy = 0.0;
  
///============================== @CGLFold ===============================================
  typedef id::TorsionType TorsionType;
  typedef std::pair< Size, TorsionType > MoveMapTorsionID;
///============================== @CGLFold ===============================================
  
//little helper function
bool contains_stageid( utility::vector1< ClassicAbinitio::StageID > vec, ClassicAbinitio::StageID query ) {
	return find( vec.begin(), vec.end(), query) != vec.end();
}




/// @detail  large (stage1/stage2)
/// small(stage2/stage3/stage4)
/// smooth_small ( stage3/stage4)
ClassicAbinitio::ClassicAbinitio(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small,
	int  /*dummy otherwise the two constructors are ambiguous */
) :
	brute_move_small_(std::move( brute_move_small )),
	brute_move_large_( brute_move_large ),
	smooth_move_small_(std::move( smooth_move_small ))
{
	BaseClass::type( "ClassicAbinitio" );
	get_checkpoints().set_type("ClassicAbinitio");
	// std::cerr << "ClassicAbinitio::constructor has stubbed out...(fatal) see code file";
	// runtime_assert( 0 ); //---> needs more implementation to use this constructor: e.g. read out movemap from FragmentMover...
	movemap_ = brute_move_large->movemap();
	//  set_defaults( pose ); in constructor virtual functions are not called
	bSkipStage1_ = false;
	bSkipStage2_ = false;

	close_chbrk_ = false;

	stage4_cycles_pack_rate_ = 0.25;
}

ClassicAbinitio::ClassicAbinitio(
	core::fragment::FragSetCOP fragset_small,
	core::fragment::FragSetCOP fragset_large,
	core::kinematics::MoveMapCOP movemap
)  :
	movemap_( movemap )
{
	BaseClass::type( "ClassicAbinitio" );
	get_checkpoints().set_type("ClassicAbinitio");
	using namespace basic::options;
	//simple_moves::ClassicFragmentMoverOP bms, bml, sms;
	using simple_moves::FragmentCostOP;
	using simple_moves::ClassicFragmentMover;
	using simple_moves::SymmetricFragmentMover;
	using simple_moves::SmoothFragmentMover;
	using simple_moves::SmoothSymmetricFragmentMover;
	using simple_moves::GunnCost;
	if ( option[ OptionKeys::abinitio::log_frags ].user() ) {
		if ( !option[ OptionKeys::abinitio::debug ] ) utility_exit_with_message( "apply option abinitio::log_frags always together with abinitio::debug!!!");
		bms = simple_moves::ClassicFragmentMoverOP( new simple_moves::LoggedFragmentMover( fragset_small, movemap ) );
		bml = simple_moves::ClassicFragmentMoverOP( new simple_moves::LoggedFragmentMover( fragset_large, movemap ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothFragmentMover( fragset_small, movemap, FragmentCostOP( new GunnCost ) ) );
	} else if ( option[ OptionKeys::abinitio::symmetry_residue ].user() ) {
		Size const sr (  option[ OptionKeys::abinitio::symmetry_residue ] );
		bms = simple_moves::ClassicFragmentMoverOP( new SymmetricFragmentMover( fragset_small, movemap, sr ) );
		bml = simple_moves::ClassicFragmentMoverOP( new SymmetricFragmentMover( fragset_large, movemap, sr ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothSymmetricFragmentMover( fragset_small, movemap, FragmentCostOP( new GunnCost ), sr ) );
	} else {
		bms = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_small, movemap ) );
		bml = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_large, movemap ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothFragmentMover ( fragset_small, movemap, FragmentCostOP( new GunnCost ) ) );
	}

	bms->set_end_bias( option[ OptionKeys::abinitio::end_bias ] ); //default is 30.0
	bml->set_end_bias( option[ OptionKeys::abinitio::end_bias ] );
	sms->set_end_bias( option[ OptionKeys::abinitio::end_bias ] );

	brute_move_small_ = bms;
	brute_move_large_ = bml;
	smooth_move_small_ = sms;

	using namespace core::pack::task;
	//init the packer
	pack_rotamers_ = minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover() );
	TaskFactoryOP main_task_factory( new TaskFactory );
	main_task_factory->push_back( operation::TaskOperationCOP( new operation::RestrictToRepacking ) );
	//main_task_factory->push_back( new operation::PreserveCBeta );
	pack_rotamers_->task_factory(main_task_factory);

	bSkipStage1_ = false;
	bSkipStage2_ = false;

	close_chbrk_ = false;

	stage4_cycles_pack_rate_ = 0.25;
}

/// @details Call parent's copy constructor and perform a shallow
/// copy of all the data.  NOTE: Shallow copy is only to preserve
/// behavior pre 9/7/2009 when the compiler-provided copy constructor
/// was being invoked.
ClassicAbinitio::ClassicAbinitio( ClassicAbinitio const & src ) :
	//utility::pointer::ReferenceCount(),
	Parent( src )
{
	stage1_cycles_ = src.stage1_cycles_;
	stage2_cycles_ = src.stage2_cycles_;
	stage3_cycles_ = src.stage3_cycles_;
	stage4_cycles_ = src.stage4_cycles_;
	stage5_cycles_ = src.stage5_cycles_;
	score_stage1_ = src.score_stage1_;
	score_stage2_ = src.score_stage2_;
	score_stage3a_ = src.score_stage3a_;
	score_stage3b_ = src.score_stage3b_;
	score_stage4_ = src.score_stage4_;
	score_stage4rot_ = src.score_stage4rot_;
	score_stage5_ = src.score_stage5_;
	apply_large_frags_ = src.apply_large_frags_;
	short_insert_region_ = src.short_insert_region_;
	just_smooth_cycles_ = src.just_smooth_cycles_;
	bQuickTest_ = src.bQuickTest_;
	close_chbrk_ = src.close_chbrk_;
	temperature_ = src.temperature_;
	movemap_ = src.movemap_;
	mc_ = src.mc_;
	brute_move_small_ = src.brute_move_small_;
	brute_move_large_ = src.brute_move_large_;
	smooth_move_small_ = src.smooth_move_small_;
	trial_large_ = src.trial_large_;
	trial_small_ = src.trial_small_;
	smooth_trial_small_ = src.smooth_trial_small_;
	total_trials_ = src.total_trials_;
	bSkipStage1_ = src.bSkipStage1_;
	bSkipStage2_ = src.bSkipStage2_;
	bSkipStage3_ = src.bSkipStage3_;
	bSkipStage4_ = src.bSkipStage4_;
	bSkipStage5_ = src.bSkipStage5_;
	recover_low_stages_ = src.recover_low_stages_;
}

/// @brief Explicit destructor is needed to destroy all the OPs
/// The compiler does all the work, but it requires that we place
/// the destructor in the .cc file.
ClassicAbinitio::~ClassicAbinitio() = default;

/// @brief setup moves, mc-object, scores
/// @details can't call this from constructor; virtual functions don't operate until construction has completed.

void
ClassicAbinitio::init( core::pose::Pose const& pose ) {
	// Parent::init( pose );
	set_defaults( pose );
	// bInitialized_ = true;
}

/// @brief ClassicAbinitio has virtual functions... use this to obtain a new instance
moves::MoverOP
ClassicAbinitio::clone() const
{
	return moves::MoverOP( new ClassicAbinitio( *this ) );
}

///======================================================================================================
///=========== @CGLFold ========================== @CGLFold ======================== @CGLFold ===========
void 
ClassicAbinitio::read_parameters(){
	ifstream informParam("./parameter_list");
	string line;
	while (getline(informParam,line)){
		if ( line[0] != '#' && line[0] != ' ' ){
			istringstream line_data( line );
			string parameterName;
			Real parameterValue;
			line_data >> parameterName >> parameterValue;
			parametersMap.insert(make_pair(parameterName, parameterValue));
		}
	}
	informParam.close();
}

void ClassicAbinitio::get_parameters(map<string, Real>& parametersMap){
	if (parametersMap.find("NP") != parametersMap.end()){
	    NP_ = parametersMap["NP"];
	    NP_reciprocal = (Real)1 / NP_;
		std::cout << "**********************NP=" << NP_ << std::endl;
	}
	else{
	    cout << "====================================parameter 'NP' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("G_LP") != parametersMap.end()){
		G_LP_ = parametersMap["G_LP"];
		std::cout << "**********************G_LP=" << G_LP_ << std::endl;
	}
	else{
		cout << "====================================parameter 'G_LP' is not found in parameters!!!" << endl;
		exit(0);
	}
	
	if (parametersMap.find("G_global") != parametersMap.end()){
		G_global_ = parametersMap["G_global"];
		std::cout << "**********************G_global=" << G_global_ << std::endl;
	}
	else{
		cout << "====================================parameter 'G_global' is not found in parameters!!!" << endl;
		exit(0);
	}
	
	if (parametersMap.find("G_local") != parametersMap.end()){
		G_local_ = parametersMap["G_local"];
		std::cout << "**********************G_local=" << G_local_ << std::endl;
	}
	else{
		cout << "====================================parameter 'G_local' is not found in parameters!!!" << endl;
		exit(0);
	}
	
	if (parametersMap.find("K") != parametersMap.end()){
		Knumber = parametersMap["K"];
		std::cout << "**********************K=" << Knumber << std::endl;
	}
	else{
		cout << "====================================parameter 'G_LP' is not found in parameters!!!" << endl;
		exit(0);
	}
}

void 
ClassicAbinitio::read_contact(){
	ifstream read_contact("./Native.disprofile");
	string contact_line;
	while (getline(read_contact, contact_line))
	{
	    Contact contact;
	    istringstream Data(contact_line);
	    
	    Data >> contact.residue1;
	    Data >> contact.residue2;
	    Data >> contact.confidence;
	    Contacts.push_back(contact);
	    
	    Cscore_full += -pow(8.0, contact.confidence);
	}
	if ( Contacts.size() == 0 ){
		cout << "ERROR!!!\n" << "ERROR : Read Contact failed!!!\n" << "ERROR!!!" << endl;
		exit(0);
	}
}


void 
ClassicAbinitio::read_distance(){
	vector<double> X(36);
	
	X[0]=2.25; X[1]=2.75; X[2]=3.25; X[3]=3.75; X[4]=4.25;
	X[5]=4.75; X[6]=5.25; X[7]=5.75; X[8]=6.25; X[9]=6.75;
	X[10]=7.25; X[11]=7.75; X[12]=8.25; X[13]=8.75; X[14]=9.25;
	X[15]=9.75; X[16]=10.25; X[17]=10.75; X[18]=11.25; X[19]=11.75;
	X[20]=12.25; X[21]=12.75; X[22]=13.25; X[23]=13.75; X[24]=14.25;
	X[25]=14.75; X[26]=15.25; X[27]=15.75; X[28]=16.25; X[29]=16.75;
	X[30]=17.25; X[31]=17.75; X[32]=18.25; X[33]=18.75; X[34]=19.25;
	X[35]=19.75;
	
	ifstream LT("distance_profile.txt");  //读取文件值
    string line;                          //设置行
	Size line_num( 0 );                   //设置行值
	Size res1( 0 );                       //存残基位置
	Size res2( 0 );
	multimap<Real, pair<pair<Size, Size>, vector<Real> > > splines_map;  //一行之中概率最大值，残基对位置，细分后的区间值
	
	while ( getline(LT, line) ){
		vector<double> Y;
		vector<string> words;
		
		if ( line[0] != ' ' ){            //判断是否是空行
			++line_num;
			if (line_num % proteinLength == 0){  //计算残基位置
				res1 = line_num / proteinLength;
				res2 = proteinLength;
			}
			else{
				res1 = line_num / proteinLength;
				++res1;
				res2 = line_num % proteinLength;
			}
			words.push_back(line);         //将这一行的值存入
		}
	
		Real p_sum( 0 );                   //一行所有值的和
		Real p_max( 0 );                   //一行所有值的最大值
		vector<string>::const_iterator it = words.begin();  //迭代
		//cout << "res1 " << res1 << " " << "res2 "<< res2 << " ";
		while (it != words.end()){
			
			istringstream line_str(*it);
			string word;
			
			while (line_str >> word){
				double d;
				stringstream s;
				s << word;
				s >> d;                    //类型转换
				Y.push_back(d);
				//cout << d << " ";
				p_sum = p_sum + d;
				if (d > p_max)
					p_max = d;
			}
			++it;
		}
		//cout << endl;
		
		if( p_sum < 0.5){                 //如果在2-20中概率和小于0.5则不考虑这个残基对的参考价值
			continue;
		}
		
		
		//tk::spline s;                     //调函数进行三次样条拟合
		//s.set_points(X,Y);
		
		//区间化
		//cout << "r1 " << res1 << " " << "r2 "<< res2 << " ";
		vector< Real > dis_map;           //存放细分后的区间值
		for(int i=0; i<40; i++){
			if( i < 4){
				dis_map.push_back( 0 );   //0-19共20个值，都设置为0.000001
			}else{
				double x=0.5*i;
				pair< Size, Size > dis_pair( x, x+0.5 );
				//cout << Y[i-4] << " ";
				dis_map.push_back( Y[i-4] );  
			}
			
		}
		//cout << endl;
		dis_map.push_back( 0 );
		pair<Size, Size> res_pair(res1, res2);
		pair<pair<Size, Size>, vector<Real> > splines( make_pair(res_pair, dis_map) );
		
		splines_map.insert( make_pair(p_max, splines) );	

	}		
	LT.close();
	
	for ( multimap<Real, pair<pair<Size, Size>, vector<Real> > >::iterator iter = splines_map.end();;){
		--iter;
		
		Size res1( iter->second.first.first );
		Size res2( iter->second.first.second );
		
		bool save = true;
		if ( fabs( res1 - res2 ) < 6 )           //去除残基间距小于6的残基对
			save = false;
		else{
			for ( Size j = 0; j < Distance_splines.size(); ++j ){
				Size e_res1( Distance_splines[j].first.first );
				Size e_res2( Distance_splines[j].first.second );
				if ( sqrt( pow((res1 - e_res1), 2) + pow((res2 - e_res2), 2) ) <= 2 ||
					 sqrt( pow((res1 - e_res2), 2) + pow((res2 - e_res1), 2) ) <= 2	||  res1 < res2){
					save = false;
					break;
				}
			}
		}
	    
	    if ( save )
		    Distance_splines.push_back( iter->second );
		
		if ( iter == splines_map.begin() )
			break;
	}
	
	cout << "Inital distance pair: " << splines_map.size() << endl;
	cout << "filtered distance pair: " << Distance_splines.size() << endl;
	
/*	
	ifstream read_dist("./filtered_dmap.txt");
	string dist_line;
	while (getline(read_dist, dist_line))
	{
	    Distance distance;
	    istringstream Data(dist_line);
	    
	    Data >> distance.res1 >> distance.res2 >> distance.dist >> distance.pro >>distance.mean >> distance.var;
	    
	    Distance_list.push_back(distance);
	}
	cout << "# Number of distance pair: " << Distance_list.size() << endl;
	if ( Distance_list.size() == 0 ){
		cout << "ERROR!!!\n" << "ERROR : Read Distance failed!!!\n" << "ERROR!!!" << endl;
		exit(0);
	}
	*/
}


void
ClassicAbinitio::read_SS(){
	ifstream read_SS("./ssmap");
	string ss_line;
	bool flag( 0 );
	while (getline(read_SS, ss_line))
	{
	    SS_unit ss;
	    istringstream Data(ss_line);
	    Size index;
	    char acid_type, ss_type;
	    Data >> index >> acid_type >> ss_type >> ss.L >> ss.H >> ss.E;
	    
	    if (index == 1 )
		    flag = true;
	    if ( flag ){
		    if ( ss_type == 'C' )
			    SS.push_back( 'L' );
		    else
			    SS.push_back( ss_type );
		    SS_map.push_back( ss );
		    
		    SSscore_full += max( max(ss.L, ss.H), ss.E);
	    }
	}
	cout << SS << "   " << SSscore_full << endl;
	if ( SS_map.size() != nativePose.total_residue() ){
		cout << "ERROR!!!\n" << "ERROR : Read Secondary Structure failed, the number of SS not equal to residue!!!\n" << "ERROR!!!" << endl;
		exit(0);
	}
}

Real 
ClassicAbinitio::SS_score(core::pose::Pose& pose){
	Real SSscore(0);
	string SS_pose = pose.secstruct();
	for (Size i = 0; i < SS_map.size(); ++i ){
		if (SS_pose[i] == 'L')
			SSscore += SS_map[i].L;
		else if (SS_pose[i] == 'H')
			SSscore += SS_map[i].H;
		else if (SS_pose[i] == 'E')
			SSscore += SS_map[i].E;
	}
	return SSscore;
}



Real 
ClassicAbinitio::Contact_score(core::pose::Pose& pose){
	Real Cscore( 0 );
	for (Size i = 0; i < Contacts.size(); ++i){
		Size residue1( Contacts[i].residue1 );
		Size residue2( Contacts[i].residue2 );
		Real confidence ( Contacts[i].confidence );

		//6距离以内不要
		if((residue2-residue1)<=6)continue;

		Real distance( 0 );
		if ( pose.residue(residue1).name3() == "GLY" || pose.residue(residue2).name3() == "GLY" ){
			numeric::xyzVector<Real> CA1 = pose.residue(residue1).xyz("CA");
			numeric::xyzVector<Real> CA2 = pose.residue(residue2).xyz("CA");
			distance = CA1.distance(CA2);
		}else{
			numeric::xyzVector<Real> CB1 = pose.residue(residue1).xyz("CB");
			numeric::xyzVector<Real> CB2 = pose.residue(residue2).xyz("CB");
			distance = CB1.distance(CB2);
		}
		//Cscore += log(abs(distance-confidence)+0.1);
		Cscore += pow(distance-confidence,2);

	}
	return Cscore;
}

Real 
ClassicAbinitio::Distance_score(core::pose::Pose& pose){
	Real Dscore( 0 );
	for (Size i = 0; i < Distance_splines.size(); ++i){
		Size residue1( Distance_splines[i].first.first );
		Size residue2( Distance_splines[i].first.second );
		vector< Real > distance_spline( Distance_splines[i].second );
		
		Real distance( 0 );
		if ( pose.residue(residue1).name3() == "GLY" || pose.residue(residue2).name3() == "GLY" ){
			numeric::xyzVector<Real> CA1 = pose.residue(residue1).xyz("CA");
			numeric::xyzVector<Real> CA2 = pose.residue(residue2).xyz("CA");
			distance = CA1.distance(CA2);
		}else{
			numeric::xyzVector<Real> CB1 = pose.residue(residue1).xyz("CB");
			numeric::xyzVector<Real> CB2 = pose.residue(residue2).xyz("CB");
			distance = CB1.distance(CB2);
		}
		
		Real distance_num( 0 );
		Size distance_index( 0 );
		
		//if( fmod( distance, 0.5 )==0 ){
		distance_index = distance/0.5;

		distance_num = distance_spline[distance_index];
		
		if( distance_index > 40 ){
			distance_num = 0;
		}
		
		//cout << "res1:" << " " << residue1 << "res2:" << " " << residue2 << "distance:" << " " << distance << "p:" << " " << distance_num << endl;
		
		Dscore += -log( max(distance_num, 0.000001) );		
	}
	
	return Dscore;
	
/*	
	///@brief Calculation of rotated point coordinates	
	Real Dscore( 0 );
	for (Size i = 0; i < Distance_list.size(); ++i){
		Size residue1( Distance_list[i].res1 );
		Size residue2( Distance_list[i].res2 );
		Real dis ( Distance_list[i].dist );
		Real residue_pro ( Distance_list[i].pro );
		Real residue_var ( Distance_list[i].var );
		
		
		Real distance( 0 );
		if ( pose.residue(residue1).name3() == "GLY" || pose.residue(residue2).name3() == "GLY" ){
			numeric::xyzVector<Real> CA1 = pose.residue(residue1).xyz("CA");
			numeric::xyzVector<Real> CA2 = pose.residue(residue2).xyz("CA");
			distance = CA1.distance(CA2);
		}else{
			numeric::xyzVector<Real> CB1 = pose.residue(residue1).xyz("CB");
			numeric::xyzVector<Real> CB2 = pose.residue(residue2).xyz("CB");
			distance = CB1.distance(CB2);
		}
		
		
		//Dscore += log( fabs( distance - dis)  + 1 ) * residue_pro;
		Dscore += log( fabs( distance - dis)  + 1 )/ residue_var;
	}
	return Dscore;
*/
}

Real 
ClassicAbinitio::Local_distance_score(core::pose::Pose& pose,vector<pair<pair<Size, Size>, vector<Real> > >& restraint_dist_map){	

	Real total_score( 0 );
	for (Size i = 0; i < restraint_dist_map.size(); ++i){
		
		Size r1( restraint_dist_map[i].first.first );
		Size r2( restraint_dist_map[i].first.second );
		
		Real distance( 0 );
		if ( pose.residue(r1).name3() == "GLY" || pose.residue(r2).name3() == "GLY" ){
			numeric::xyzVector<Real> CA1 = pose.residue(r1).xyz("CA");
			numeric::xyzVector<Real> CA2 = pose.residue(r2).xyz("CA");
			
			distance = CA1.distance(CA2);
		}else{
			numeric::xyzVector<Real> CB1 = pose.residue(r1).xyz("CB");
			numeric::xyzVector<Real> CB2 = pose.residue(r2).xyz("CB");
			
			distance = CB1.distance(CB2);
		}
		
		//total_score +=  restraint_dist_map[i].second[2] * exp( fabs(restraint_dist_map[i].second[0] - distance)) ;
		total_score += ( restraint_dist_map[i].second[2] * log( fabs(restraint_dist_map[i].second[0] - distance) +1 )/ restraint_dist_map[i].second[1] );
	}
	return total_score;
}

Real 
ClassicAbinitio::Local_distance_score1(core::pose::Pose& pose,vector<pair<pair<Size, Size>, vector<Real> > >& restraint_dist_map){	

	Real total_score( 0 );
	for (Size i = 0; i < restraint_dist_map.size(); ++i){
		
		Size r1( restraint_dist_map[i].first.first );
		Size r2( restraint_dist_map[i].first.second );
		
		Real distance( 0 );
		if ( pose.residue(r1).name3() == "GLY" || pose.residue(r2).name3() == "GLY" ){
			numeric::xyzVector<Real> CA1 = pose.residue(r1).xyz("CA");
			numeric::xyzVector<Real> CA2 = pose.residue(r2).xyz("CA");
			cout << "Gr1" <<  "("  << CA1.x() << "," << CA1.y() << ","<< CA1.z()<<  ")"
			<< "Gr2" <<  "("  << CA2.x() << "," << CA2.y() << ","<< CA2.z()<<  ")"	;
			distance = CA1.distance(CA2);
		}else{
			numeric::xyzVector<Real> CB1 = pose.residue(r1).xyz("CB");
			numeric::xyzVector<Real> CB2 = pose.residue(r2).xyz("CB");
			cout << "Gr1" <<  "("  << CB1.x() << "," << CB1.y() << ","<< CB1.z()<<  ")"
			<< "Gr2" <<  "("  << CB2.x() << "," << CB2.y() << ","<< CB2.z()<<  ")"	;
			distance = CB1.distance(CB2);
		}
		
		//total_score +=  restraint_dist_map[i].second[2] * exp( fabs(restraint_dist_map[i].second[0] - distance)) ;
		total_score += ( restraint_dist_map[i].second[2] * log( fabs(restraint_dist_map[i].second[0] - distance) +1 )/ restraint_dist_map[i].second[1] );
	}
	cout << endl;
	return total_score;
}


Real 
ClassicAbinitio::Average_Cscore_of_Populatin(vector< Real >& population_Cscore){
	Real total_Cscore( 0 );
	for (Size i = 0; i < NP_; ++i)
		total_Cscore += population_Cscore[i];
	Real average_Cscore( total_Cscore * NP_reciprocal );
	
	return average_Cscore;
}


Real 
ClassicAbinitio::Average_Cscore_of_Populatin(Real& old_Cscore, Real& new_Cscore, Real& average_Cscore){
	return average_Cscore + (new_Cscore - old_Cscore) * NP_reciprocal;
}

//LT
Real 
ClassicAbinitio::Average_Dscore_of_Populatin(vector< Real >& population_Dscore){
	Real total_Dscore( 0 );
	for (Size i = 0; i < NP_; ++i)
		total_Dscore += population_Dscore[i];
	Real average_Dscore( total_Dscore * NP_reciprocal );
	
	return average_Dscore;
}


Real 
ClassicAbinitio::Average_Dscore_of_Populatin(Real& old_Dscore, Real& new_Dscore, Real& average_Dscore){
	return average_Dscore + (new_Dscore - old_Dscore) * NP_reciprocal;
}

bool 
ClassicAbinitio::boltzmann_accept(Real const& targetEnergy, Real const& trialEnergy){
	if ( trialEnergy <= targetEnergy )
		return true;
	else{
		if ( exp( -(trialEnergy - targetEnergy) * 0.5 ) < numeric::random::rg().uniform() )
			return false;
		else
			return true;
	}
}

bool 
ClassicAbinitio::boltzmann_accept(const Real& targetEnergy, const Real& trialEnergy, const Real& recipocal_KT){
	/// recipocal_KT = 1 / KT
	if ( trialEnergy <= targetEnergy )
		return true;
	else{
		if ( exp( -log(abs(trialEnergy - targetEnergy)+0.1) * recipocal_KT ) < numeric::random::rg().uniform() )
			return false;
		else
			return true;
	}
}

void 
ClassicAbinitio::Output_SPICKER_All_Data(core::pose::Pose& pose, const Real& energy){
	++Spicker_All_num;
	if ( Spicker_All_num > 20000 ){
		Spicker_All_num = 1;
		SPICKER_All_Data.close();
		++Spicker_All_file_num;
		string file_name;
		stringstream ss;
		ss << "./output_files/SPICKER_data/spicker.data" << Spicker_All_file_num;
		ss >> file_name;
		SPICKER_All_Data.open( file_name.c_str() );
	}
	Size Num_resdue( pose.total_residue() );
	SPICKER_All_Data << right << setw(8) << Num_resdue
		    << right << setw(10) << setprecision(3) << fixed << energy
		    << right << setw(8) << Spicker_All_num
		    << right << setw(8) << Spicker_All_num
		    << endl;
	for ( Size r = 1; r <= Num_resdue; ++r ){
		numeric::xyzVector<Real> CA_ = pose.residue(r).xyz("CA");
		
		SPICKER_All_Data << right << setw(10) << setprecision(3) << fixed << CA_.x()
			    << right << setw(10) << setprecision(3) << fixed << CA_.y()
			    << right << setw(10) << setprecision(3) << fixed << CA_.z()
			    << endl;
	}
}

void 
ClassicAbinitio::SPICKER_Demand_All(core::pose::Pose& nativePose){
	SPICKER_All_Data.close();
  
	ofstream tra("./output_files/SPICKER_data/tra.in");
	tra << " " << Spicker_All_file_num << " 1 1" << endl;
	for (Size n = 1; n <= Spicker_All_file_num; ++n)
		tra << "spicker.data" << n << endl;
	tra.close();
	
	ofstream rmsinp("./output_files/SPICKER_data/rmsinp");
	rmsinp << 1 << "  " << proteinLength << "    ! these two numbers indicate region for RMSD calculations [1," << proteinLength << "]"
		<< "\n" << proteinLength << "       ! protein length";
	rmsinp.close();
	
	ofstream seq("./output_files/SPICKER_data/seq.dat");
	for (Size r = 1; r <= proteinLength; ++r)
		seq << right << setw(5) << r 
		    << right << setw(6) << nativePose.residue(r).name3() 
		    << endl;
	seq.close();
	
	ofstream CA("./output_files/SPICKER_data/CA");
	for (Size r = 1; r <= proteinLength; ++r){
		numeric::xyzVector<Real> CA_ = nativePose.residue(r).xyz("CA");
	  
		CA << "ATOM"
		  << right << setw(7) << r
		  << right << setw(4) << "CA"
		  << right << setw(5) << nativePose.residue(r).name3()
		  << right << setw(6) << r
		  << right << setw(12) << CA_.x()
		  << right << setw(8) << CA_.y()
		  << right << setw(8) << CA_.z()
		  << endl;
	}
	CA.close();
}


vector< core::pose::Pose > 
ClassicAbinitio::generate_population_random(const Size& population_size, const Size& seq_length){
	
	vector< core::pose::Pose > population;
	core::pose::Pose init_pose(initPose);
	Real new_pose_phi(0);
	Real new_pose_psi(0);
	
	for (Size i = 0; i < population_size; ++i)
	{
		core::pose::Pose new_pose( init_pose );
		for (Size j = 1; j <= seq_length; ++j) 
		{
			new_pose_phi = numeric::random::rg().random_range(-180, 180);
			new_pose_psi = numeric::random::rg().random_range(-180, 180);
			
			norm( new_pose_phi );
			norm( new_pose_psi );
			
			new_pose.set_phi( j, new_pose_phi );
			new_pose.set_psi( j, new_pose_psi );
		}
		population.push_back(new_pose);
	}
	
	return population;
}

void
ClassicAbinitio::learn_period(vector<core::pose::Pose> &population, Size const &population_size, const Size& seq_length, const Size& G_length, scoring::ScoreFunction const & scorefxn){
	
	Size g = 0;
	Real score = 0;
	PopulationRmsd = calculate_population_Rmsd( population );
	PopulationEnergy = calculate_population_energy( population, scorefxn );
	bool crowd_flag = false;
	int count = 0;
	int accept = 0;
	ofstream Ave_acc;
	Ave_acc.open("./Ave_acc_LP.csv");
	ofstream Every;
	Every.open("./every_energy.csv");
	ofstream Ave_energy;
	Ave_energy.open("./Ave_energy.csv");
	
	while( g < G_length )
	{
		std::cout << "=======generation " << g << "=======" << std::endl;
		for(Size i = 0; i < population_size; i++)
		{
			core::pose::Pose trial( population[i] );
			mutation_rand_1( population, trial, population_size, seq_length, i );
			score = scorefxn(trial);
			
			Size Index_parent;
			vector<Real> dist_parent;
			
			for(Size s = 0; s < population_size; s++)
			{
				dist_parent.push_back( core::scoring::CA_rmsd( trial, population[s] ) );
			}
			
			vector<Real>::iterator min_parent = min_element(dist_parent.begin(),dist_parent.end()); 
			Index_parent=distance(dist_parent.begin(), min_parent);
				
			if( score <= PopulationEnergy[Index_parent] )
			{
				population[Index_parent] = trial;
				PopulationEnergy[Index_parent] = score;
				PopulationRmsd[Index_parent] = core::scoring::CA_rmsd( trial, nativePose );
				accept++;
				
				Every << "ESELP," << g << "," << Index_parent << "," << PopulationRmsd[Index_parent] << "," << PopulationEnergy[Index_parent] << std::endl;
			}
				
			//std::cout << "trial " << i << " try to parents ";
			//std::cout << std::endl;
			//std::cout << "Finally trial " << i << " replace the parent " << Index_parent << std::endl;
		}
		
		Ave_energy << "ESELP," << g << "," << Average_Dscore_of_Populatin( PopulationEnergy ) << "," << Average_Dscore_of_Populatin( PopulationRmsd ) << endl;
		Ave_acc <<  Average_Dscore_of_Populatin( PopulationEnergy ) << ","
				<<  Average_Dscore_of_Populatin( PopulationRmsd ) << ","<<accept;
		Ave_acc << endl;
		std::cout << "In generation " << g << " have " << accept << " pose update" << std::endl;
		accept = 0;
		g++;
	}
	Output_last_generation( population, PopulationEnergy, 'P');
	Ave_acc.close();
	Every.close();
}

vector<core::pose::Pose>
ClassicAbinitio::k_mediods( const Size& k_number, vector<core::pose::Pose> &population, Size const &population_size, vector<Real> &state_tag){

	state_tag.resize(population_size,0);
	vector<Real> medoids_tag;
	vector<core::pose::Pose> medoids;
	vector< vector<Real> > distmatrix;
	vector<Real> dist(population_size,0);
	for(Size i=0; i<population_size; i++)
		distmatrix.push_back(dist);
	
	std::cout << "Index of medoids in initial: ";
	Size i=0;
	while( i<k_number )
	{
		medoids_tag.push_back(numeric::random::rg().random_range(0, population_size-1));
		std::cout << medoids_tag[i] << "\t";
		i++;
	}
	std::cout << std::endl;
	distmatrix = compute_distmatrix( population, population, population_size, population_size );
	
	bool clusterChanged = true;
	Size num=0;//update number
	while(clusterChanged && num<100)//clusterChanged && num<1000
	{
		clusterChanged = false;

		//step one : find the nearest mediods of each point
		//cout<<"find the nearest mediods of each point : "<<endl;
		for(Size i=0; i<population_size; i++)//for each individual
		{
			Size minIndex = 0;
			Size minDist = 100000;
			Size j=0;
			while(j<k_number)//for each medoids
			{
				if( distmatrix[i][medoids_tag[j]] < minDist )
				{
					minDist = distmatrix[i][medoids_tag[j]]; 
					minIndex = j;//record the nearest cluster of current individual
				}
				j++;
			}
			//update assignment of each individual
			//std::cout << i << " belong to medoids " << minIndex << std::endl;
			state_tag[i] = minIndex;
			std::cout << state_tag[i] << ", ";
		}
		std::cout << std::endl;
		//step two : update the mediods
		std::cout << "update the mediods:" << std::endl;
		vector<Real> medoids_tag2;
		vector<Real> tag_cluster;//record index in each cluster
		for(Size med=0; med<k_number; med++)//for each cluster
		{
			Size minerror;
			for(Size i=0; i<population_size; i++)//for each individual
			{
				if( state_tag[i] == med )//if same cluster
				{
					//std::cout << i << " belong to mediods " << med << std::endl;
					tag_cluster.push_back(i);//record index in each cluster
				}
			}
			std::cout << med << " cluster have the number is " << tag_cluster.size() << std::endl;
			std::cout << med << " cluster have data are :" << std::endl;
			for(Size i=0; i<tag_cluster.size(); i++)
				std::cout << tag_cluster[i] << "\t";
			std::cout << std::endl;
			
			Size index=0;
			for(Size i=0; i<tag_cluster.size(); i++)
			{
				Size error=0;//save absolutely error, assign zero to ensure that the next data is not affected

				for(Size j=0; j<tag_cluster.size(); j++)
				{
					error += distmatrix[tag_cluster[i]][tag_cluster[j]];
				}
				if(i==0)
				{
					minerror = error;
					index = tag_cluster[i];
				}
				if(error<minerror)
				{
					minerror = error;
					index = tag_cluster[i];
				}
			}
			//cout<<endl;
			medoids_tag2.push_back(index);
			tag_cluster.clear();//ensure that the next cluster is not affected
			
		}//for(each cluster)
		if(medoids_tag != medoids_tag2)
		{
			clusterChanged = true;
			//cout<<num<<" time is true"<<endl;
		}
		else
		{
			//cout<<num<<" time is false"<<endl;
		}
		medoids_tag.assign(medoids_tag2.begin(),medoids_tag2.end());//2 to 1
		
		std::cout << "medoids: " << std::endl;
		for(Size i=0; i<k_number; i++)
			std::cout << medoids_tag[i] << "\t";
		std::cout << std::endl;
		
		medoids_tag2.clear();
		num++;
		//cout<<"the number of cluster: "<<mediods_tag.size()<<endl;
	}//while
	
	std::cout << "update times: " << num << std::endl;
	std::cout << "Finally medoids: " << std::endl;
	for(Size i=0; i<k_number; i++)
	{
		std::cout << medoids_tag[i] << "\t";
		medoids.push_back( population[medoids_tag[i]] );
	}
	std::cout << std::endl;
	
	return medoids;
}

vector< vector<Real> >
ClassicAbinitio::compute_distmatrix(vector<core::pose::Pose> &population1, vector<core::pose::Pose> &population2, Size const &population_size1, Size const &population_size2){
	
	vector<Real> dist;
	vector< vector<Real> > distmatrix;
	for(Size i=0; i<population_size1; i++)
	{
		core::pose::Pose first( population1[i] );
		for(Size j=0; j<population_size2; j++)
		{
			core::pose::Pose second( population2[j] );
			dist.push_back( core::scoring::CA_rmsd( first, second ) );
		}
		distmatrix.push_back(dist);
		dist.clear();
	}
	return distmatrix;
}

void
ClassicAbinitio::entropy_switch_strategy(vector<core::pose::Pose> &population, Size const &population_size, const Size& seq_length, const Size& G_global, const Size& G_local, scoring::ScoreFunction const & scorefxn, vector<Real> &state_tag, const Size& k_number, vector<core::pose::Pose> &mediods_set){
	
	Size g = 0;
	Real score = 0;
	Real rand_num = 0;
	Real odds_max = 0;
	PopulationRmsd = calculate_population_Rmsd( population );
	PopulationEnergy = calculate_population_energy( population, scorefxn );
	vector<double> entropy;
	entropy.push_back(0);
	vector<Real> tag_parent(population_size,0);
	tag_parent = state_tag;
	vector<Real> tag_child(population_size,0);
	
	vector<Real> trans_temp(k_number,0);
	vector< vector<Real> > trans;
	vector<Real> frequency_temp(k_number,0);
	vector< vector<Real> > frequency;
	vector<Real> record_temp(k_number,0);
	vector< vector<Real> > record;
	for(Size i=0; i<k_number; i++)
	{
		trans.push_back(trans_temp);
		frequency.push_back(frequency_temp);
		record.push_back(record_temp);
	}
	
	ofstream Ave_entropy;
	Ave_entropy.open("./Ave_entropy.csv");
	ofstream Ave_acc;
	Ave_acc.open("./Ave_acc_main.csv");
	ofstream cluster_process;
	cluster_process.open("./cluster_process.txt");
	ofstream Every("./every_energy.csv", ios::app);
	ofstream Ave_energy("./Ave_energy.csv", ios::app);
	ofstream global_FE;
	global_FE.open("./global_FE.txt");
	
	cluster_process << "Inital cluster: " << endl;
	for(Size i=0; i<population_size; i++)
	{
		cluster_process << state_tag[i] << ", ";
	}
	cluster_process << endl;
	
	bool isFinish = false;
	Real accept_num = 0;
	Real accept_num_probability = 0;
	
	while( g < G_global )
	{
		
		for(Size i = 0; i < population_size; i++)
		{
			core::pose::Pose trial( population[i] );
			Real a(G_global);
			Real b(g);
			
			//main algorithm
			if( g==0 || g==1 ) //g==0 || g==1
			{
				std::cout << "generation " << g << " Entropy = " << entropy[g] << ", pose " << i << ": mutation_rand_1" << std::endl;
				mutation_rand_1( population, trial, population_size, seq_length, i );
			}
			else
			{
				if( entropy[g] == 0 )
				{
					Output_last_generation( population, PopulationEnergy, 'G');
					//std::cout << "generation " << g << " Entropy = " << entropy[g] << ", pose " << i << ": minmover" << std::endl;
					//cluster_process << "generation " << g << " Entropy=" << entropy[g] << ", pose " << i << ": minmover" << std::endl;
					minmover( population, population_size, seq_length, G_local, scorefxn );
					isFinish = true;
					break; //out of i<NP
				}
				else
				{
					rand_num = numeric::random::rg().uniform();
					if( rand_num < entropy[g] )
					{
						//std::cout << "generation " << g << " Entropy = " << entropy[g] << ", pose " << i << ": mutation_rand_1" << std::endl;
						cluster_process << "generation " << g << " Entropy=" << entropy[g] << ", random=" << rand_num << ", pose " << i << ": mutation_rand_1" << std::endl;
						mutation_rand_1( population, trial, population_size, seq_length, i );
						//convex_mutation2(population, PopulationEnergy, trial, population_size, proteinLength, i, b/a);
					}
					else
					{
						//std::cout << "generation " << g << " Entropy = " << entropy[g] << ", pose " << i << ": mutation_2" << std::endl;
						cluster_process << "generation " << g << " Entropy=" << entropy[g] << ", random=" << rand_num << ", pose " << i << ": mutation_2" << std::endl;
						//mutation_best_1( population, trial, population_size, seq_length, i );
						//convex_mutation2(population, PopulationEnergy, trial, population_size, proteinLength, i, b/a);
						mutation_current2rand_1( population, trial, population_size, seq_length, i );
					}
				}
			}
			score = scorefxn(trial);
			if( score <= PopulationEnergy[i] )
			{
				population[i] = trial;
				PopulationEnergy[i] = score;
				PopulationRmsd[i] = core::scoring::CA_rmsd( trial, nativePose );
				//std::cout << "accept" << std::endl;
				cluster_process << "accept" << std::endl;
				accept_num ++;
				
				Every << "global," << g << "," << i << "," << PopulationRmsd[i] << "," << PopulationEnergy[i] << std::endl;
			}
			else
			{
				rand_num = numeric::random::rg().uniform();
				//std::cout << "rand=" << rand_num << ", expdE=" << exp( (PopulationEnergy[i] - score) * entropy[g]) << ", dE=" << PopulationEnergy[i] - score << ", entropy=" << entropy[g] << std::endl;
				cluster_process << "rand=" << rand_num << ", expdE=" << exp( (PopulationEnergy[i] - score) * entropy[g]) << ", dE=" << PopulationEnergy[i] - score << ", entropy=" << entropy[g] << std::endl;
				if( rand_num < exp( (PopulationEnergy[i] - score) * entropy[g]) )
				{
					population[i] = trial;
					PopulationEnergy[i] = score;
					PopulationRmsd[i] = core::scoring::CA_rmsd( trial, nativePose );
					//std::cout << "probability accept" << std::endl;
					cluster_process << "probability accept" << std::endl;
					accept_num ++;
					accept_num_probability ++;
				}
				else
				{
					//std::cout << "reject" << std::endl;
					cluster_process << "reject" << std::endl;
				}
			}
		}
		
		if(isFinish == true)
			break;
		
		std::cout << "generation " << g << " accept " << accept_num << ", probability accept " << accept_num_probability << std::endl;
		cluster_process << "generation " << g << " accept " << accept_num << ", probability accept " << accept_num_probability << std::endl;
		accept_num = 0;
		accept_num_probability = 0;
		
		vector<Real> state_num(k_number,0);
		for(Size i=0; i<k_number; i++)
		{
			for(Size j=0; j<population_size; j++)
			{
				if(tag_parent[j]==i)
					state_num[i]++;
			}
			std::cout << i << " cluster have " << state_num[i] << " pose" << std::endl;
			cluster_process << i << " cluster have " << state_num[i] << " pose" << endl;
		}
			
		cluster(population, population_size, state_tag, k_number, mediods_set);
		tag_child = state_tag;
		std::cout << "cluster state_tag: ";
		cluster_process << "generation " << g << " cluster state_tag: " << endl;
		for(Size i=0; i<population_size; i++)
		{
			std::cout << state_tag[i] << ", ";
			cluster_process << state_tag[i] << ", ";
		}
		std::cout << std::endl;
		cluster_process << endl;
		
		Real sum = 0;
		Real sum_trans = 0;
		Real sum_notrans = 0;
		for(Size s=0; s<population_size; s++)
		{
			trans[ tag_parent[s] ][ tag_child[s] ] ++;
			frequency[ tag_parent[s] ][ tag_child[s] ] = ( trans[ tag_parent[s] ][ tag_child[s] ] )/(population_size);
			record[ tag_parent[s] ][ tag_child[s] ] = frequency[ tag_parent[s] ][ tag_child[s] ] * log( frequency[ tag_parent[s] ][ tag_child[s] ] );
		}
		cluster_process << "trans: " << endl;
		for(Size i=0; i<k_number; i++)
		{
			for(Size j=0; j<k_number; j++)
			{
				cluster_process << trans[i][j] << "\t";
			}
			cluster_process << endl;
		}
		cluster_process << "frequency: " << endl;
		for(Size i=0; i<k_number; i++)
		{
			for(Size j=0; j<k_number; j++)
			{
				cluster_process << frequency[i][j] << "\t";
			}
			cluster_process << endl;
		}
		cluster_process << "record: " << endl;
		for(Size i=0; i<k_number; i++)
		{
			for(Size j=0; j<k_number; j++)
			{
				cluster_process << record[i][j] << "\t";
			}
			cluster_process << endl;
		}
		
		for(Size i=0; i<k_number; i++)
		{
			for(Size j=0; j<k_number; j++)
			{
				if(frequency[i][j]!=0)
					sum += frequency[i][j] * log(frequency[i][j]);

				if(i!=j)
				{
					if(frequency[i][j]!=0)
						sum_trans += frequency[i][j] * log(frequency[i][j]);
					else
						sum_trans += 0;
				}
				else
				{
					if(frequency[i][j]!=0)
						sum_notrans += frequency[i][j] * log(frequency[i][j]);
					else
						sum_notrans += 0;
				}
			}
		}
		Real Entropy = - sum;
		Real Entropy_trans = - sum_trans;
		Real Entropy_notrans = - sum_notrans;
		Real odds = Entropy_trans / Entropy_notrans;
		
		if(g==0)
			odds_max = odds;
		Real norm_odds = odds/odds_max;
		
		Real E_Max = 0, p_Max = 0;
		p_Max = 1/(k_number * k_number * 1.0);
		for(Size i=0; i<(k_number*k_number); i++)
			E_Max += ( - p_Max * log(p_Max) );
		Real norm_Entropy = Entropy/E_Max;
		
		entropy.push_back(norm_odds); //norm_odds,or norm_Entropy
		/*
		if(g==0)
		{
			Output_last_generation( population, PopulationEnergy, '0');
		}*/
		
		std::cout << "Entropy_" << g << " = " << norm_Entropy << ",  " << Entropy << ", " << E_Max << std::endl;
		std::cout << "norm_odds=" << norm_odds << ", odds_max=" << odds_max << ", odds=" << odds << ", Entropy_trans=" << Entropy_trans << ", Entropy_notrans=" << Entropy_notrans << std::endl;
		cluster_process << "Entropy_" << g << " = " << norm_Entropy << ",  " << Entropy << ", " << E_Max << std::endl;
		cluster_process << "norm_odds=" << norm_odds << ", odds_max=" << odds_max << ", odds=" << odds << ", Entropy_trans=" << Entropy_trans << ", Entropy_notrans=" << Entropy_notrans << std::endl;
		//Ave_entropy << norm_Entropy << "," << odds << "," << norm_odds << endl;
		Ave_entropy << "global," << g << "," << norm_odds << "," << odds << "," << norm_Entropy << endl;
		Ave_energy  << "global," << g << "," << Average_Dscore_of_Populatin( PopulationEnergy ) << "," << Average_Dscore_of_Populatin( PopulationRmsd ) << endl;
		Ave_acc <<  Average_Dscore_of_Populatin( PopulationEnergy ) << ","
				<<  Average_Dscore_of_Populatin( PopulationRmsd ) ;
		Ave_acc << endl;
		
		g++;
		
		tag_parent.clear();
		tag_parent = tag_child;
		tag_child.clear();
		for(Size i=0; i<k_number; i++)
		{
			for(Size j=0; j<k_number; j++)
			{
				trans[i][j] = 0;
				frequency[i][j] = 0;
				record[i][j] = 0;
			}
		}
		
	}
	
	if( g == G_global)
		minmover( population, population_size, seq_length, G_local, scorefxn );
	
	global_FE << g <<endl;
	
	//Output_last_generation( population, PopulationEnergy, 'G');
	Ave_entropy.close();
	Ave_acc.close();
	cluster_process.close();
	Every.close();
	global_FE.close();
	
}

void
ClassicAbinitio::cluster(vector<core::pose::Pose> &population, Size const &population_size, vector<Real> &state_tag, const Size& k_number, vector<core::pose::Pose> &mediods_set){
	
	//std::cout << "cluster state_tag: ";
	for(Size i=0; i<population_size; i++)
	{
		vector<Real> dist_parent;
		for(Size j=0; j<k_number; j++)
		{
			dist_parent.push_back( core::scoring::CA_rmsd(population[i],mediods_set[j]) );
		}
		vector<Real>::iterator min_parent = min_element(dist_parent.begin(),dist_parent.end()); 
		state_tag[i]=distance(dist_parent.begin(), min_parent);
		//std::cout << state_tag[i] << ", ";
	}
	//std::cout << std::endl;
}

void
ClassicAbinitio::minmover(vector<core::pose::Pose> &population, Size const &population_size, const Size& seq_length, const Size& G_local, scoring::ScoreFunction const & scorefxn){
	
	using namespace core;
	using namespace protocols;
	using namespace protocols::moves;
	using namespace core::scoring;

	using namespace pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	ofstream Every("./every_energy.csv", ios::app);
	ofstream Ave_energy("./Ave_energy.csv", ios::app);
	clock_t Initime = clock();
	
	std::string min_type =  option[ OptionKeys::run::min_type ]();
	core::Real min_tol =  option[ OptionKeys::run::min_tolerance ]();
	core::optimization::MinimizerOptions options( min_type, min_tol, true, false );
	core::kinematics::MoveMap final_mm;
	final_mm.set_chi(  false  );
	final_mm.set_bb(  true  );
	final_mm.set_jump(  false );

	for (Size m=0; m<population_size; m++)
	{
		cout << "# add constraints " << m << endl;
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( population[m] );
	}

	PopulationRmsd = calculate_population_Rmsd( population );
	PopulationEnergy = calculate_population_energy( population, scorefxn );

	Size accept_count(0);

	cout << "====================================== Start algorithm ======================================" << endl;
	for(Size g = 1; g <= G_local; ++g)
	{
		cout << "#  Generation: " << g << endl;
		accept_count=0;

		for(Size m = 0; m < population_size; m++)
		{
			Real a(G_local);
			Real b(g);

			core::pose::Pose Ctrial( population[m] );

			Real pre_enegy(scorefxn(Ctrial));

//------------------------------------mutation-------------------------------------
			//convex_mutation1(population, Ctrial,NP_, proteinLength, m);
			//mutation_rand_1(population, Ctrial,NP_, proteinLength, m);
			convex_mutation2(population, PopulationEnergy, Ctrial, population_size, seq_length, m,b/a);
			Real muta_enegy(scorefxn(Ctrial));

//------------------------------------minmover-------------------------------------

			core::optimization::AtomTreeMinimizer().run( Ctrial, final_mm, scorefxn, options );
			Real min_enegy(scorefxn(Ctrial));

//------------------------------------selection-------------------------------------

			if( min_enegy < pre_enegy )
			{
				population[m] = Ctrial;
				PopulationRmsd[m]=core::scoring::CA_rmsd( Ctrial, nativePose );
				PopulationEnergy[m]=min_enegy;
				accept_count++;
				
				Every << "local," << g << "," << m << "," << PopulationRmsd[m] << "," << PopulationEnergy[m] << endl;

			}

		}

		Ave_energy << "local," << g << "," << Average_Dscore_of_Populatin( PopulationEnergy ) << "," << Average_Dscore_of_Populatin( PopulationRmsd ) << endl;
		Ave_acc <<  Average_Dscore_of_Populatin( PopulationEnergy ) << ","
				<<  Average_Dscore_of_Populatin( PopulationRmsd ) << ","<<accept_count;
		Ave_acc << endl;
		std::cout << "G:" << g << "---accept_count---" << accept_count << std::endl;

	}
	Output_last_generation( population, PopulationEnergy, 'L');
	
	Ave_energy.close();
	Every.close();

	clock_t endtime = clock();
	tr << "Time about local stage: " << (double(endtime) - Initime )/( CLOCKS_PER_SEC ) << " seconds." << std::endl;
	
	
}

vector< core::pose::Pose > 
ClassicAbinitio::generate_population(const Size& population_size){
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;
	
	vector< core::pose::Pose > population;
	core::pose::Pose init_pose(initPose);
	//core::scoring::constraints::add_constraints_from_cmdline_to_pose( init_pose );
	bool success(true);
	
	cout << "#  Generating the initial population." << endl;
	
	if ( !bSkipStage1_ ) {
	//	stage3_to_stage4 = false;
		PROF_START( basic::STAGE1 );
		clock_t starttime = clock();
		if ( !prepare_stage1( init_pose ) ) {
			cout << "ERROR!!!\n" << "ERROR : prepare_stage1 failed!!!\n" << "ERROR!!!" << endl;
			exit(0);
		}
	//	tr.Info <<  "\n===================================================================\n";
	//	tr.Info <<  "   Stage 1                                                         \n";
	//	tr.Info <<  "   Folding with score0 for max of " << stage1_cycles() << std::endl;

		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( init_pose, "stage0" );
		}
		if ( !get_checkpoints().recover_checkpoint( init_pose, get_current_tag(), "stage_1", false /* fullatom*/, true /*fold tree */ ) ) {
			ConstraintSetOP orig_constraints(NULL);
			orig_constraints = init_pose.constraint_set()->clone();
			
///=======================================用户添加程序=====================================================================
			for (Size i = 0; i < population_size; ++i){
			    core::pose::Pose new_pose( init_pose );
			    mc_->reset( new_pose );
			    success = do_stage1_cycles( new_pose );
			    recover_low( new_pose, STAGE_1 );
			    population.push_back(new_pose);
			}
///======================================================================================================================
			if ( tr.Info.visible() ) current_scorefxn().show( tr, init_pose );
			mc().show_counters();
			total_trials_+=mc().total_trials();
			mc().reset_counters();

			init_pose.constraint_set( orig_constraints ); // restore constraints - this is critical for checkpointing to work
			get_checkpoints().checkpoint( init_pose, get_current_tag(), "stage_1", true /*fold tree */ );
		} //recover checkpoint
		get_checkpoints().debug( get_current_tag(), "stage_1", current_scorefxn()( init_pose ) );

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE1 );
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
	//		tr.Info << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
			output_debug_structure( init_pose, "stage1" );
		}
	//	tr.Info <<  "\n===================================================================\n";
	//	tr.Info <<  "  Finished Abinitio Starge1!!!" << endl;
	//	tr.Info <<  endl;
	} //skipStage1
	if ( !success ) {
		cout << "ERROR!!!\n" << "ERROR : prepare_stage1 failed!!!\n" << "ERROR!!!" << endl;
		exit(0);
	}
	if ( !bSkipStage2_ ) {
	//	stage3_to_stage4 = true;
	//	
	//	tr.Info <<  "\n===================================================================\n";
	//	tr.Info <<  "   Stage 2                                                         \n";
	//	tr.Info <<  "   Folding with score1 for " << stage2_cycles() << std::endl;

		PROF_START( basic::STAGE2 );
		clock_t starttime = clock();
		if ( close_chbrk_ ) {
			Real const setting( 0.25 );
			set_score_weight( scoring::linear_chainbreak, setting, STAGE_2 );
	//		tr.Info <<  " Chain_break score assigned " << std::endl;
		}
		if ( !prepare_stage2( init_pose ) )  {
			cout << "ERROR!!!\n" << "ERROR : prepare_stage1 failed!!!\n" << "ERROR!!!" << endl;
			exit(0);
		}
		if ( !get_checkpoints().recover_checkpoint( init_pose, get_current_tag(), "stage_2", false /* fullatom */, true /*fold tree */ ) ) {
			ConstraintSetOP orig_constraints(NULL);
			orig_constraints = init_pose.constraint_set()->clone();	
			
///=======================================用户添加程序===============================stage2================================
			for (Size i = 0; i < NP_; ++i){
				mc_->reset( population[i] );
				success = do_stage2_cycles( population[i] );
				recover_low( population[i], STAGE_2 );
			}
///======================================================================================================================                 
			if  ( tr.visible() ) current_scorefxn().show( tr, init_pose );
			mc().show_counters();
			total_trials_+=mc().total_trials();
			mc().reset_counters();
			init_pose.constraint_set( orig_constraints ); // restore constraints - this is critical for checkpointing to work
			get_checkpoints().checkpoint( init_pose, get_current_tag(), "stage_2", true /*fold tree */ );
		}
		get_checkpoints().debug( get_current_tag(), "stage_2", current_scorefxn()( init_pose ) );

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE2 );
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( init_pose, "stage2" );
	//		tr << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
		}
	//	tr.Info <<  "\n===================================================================\n";
	//	tr.Info <<  "  Finished Abinitio Starge2!!!" << endl;
	//	tr.Info <<  endl;
	} //bSkipStage2
	if ( !success ) {
		cout << "ERROR!!!\n" << "ERROR : prepare_stage1 failed!!!\n" << "ERROR!!!" << endl;
		exit(0);
	}



	return population;
}

vector< Real > 
ClassicAbinitio::calculate_population_energy_2(vector< core::pose::Pose >& population){
	vector<Real> population_energy;
	
	for (Size i = 0; i < NP_; ++i){
		TrialEnergy = (*score_stage2_)(population[i]);
		if ( find_if( population_energy.begin(), population_energy.end(), isSimilar ) == population_energy.end() ){
			population_energy.push_back(TrialEnergy);
		}
		else{
		//	cout << "The TrialPose is Similar to The Pose in Population!!!" << endl;
 		//	cout << "Restarting a New Trajectory :" << endl;
			bool flag = false;
			while ( !flag ){
				core::pose::Pose new_pose( generate_a_new_pose() );
				
				TrialEnergy = (*score_stage2_)(new_pose);
				
				if ( find_if( population_energy.begin(), population_energy.end(), isSimilar ) == population_energy.end() ){
					population[i] = new_pose;
					population_energy.push_back(TrialEnergy);
					flag = true;
				}
			}
		}
	}
	return population_energy;
}

vector< Real > 
ClassicAbinitio::calculate_population_energy(vector< core::pose::Pose >& population,scoring::ScoreFunction const & scorefxn){
	vector<Real> population_energy;
	
	for (Size i = 0; i < NP_; ++i){
		TrialEnergy = scorefxn(population[i]);
		population_energy.push_back(TrialEnergy);
	}
	return population_energy;
}

vector< Real > 
ClassicAbinitio::calculate_population_Cscore(vector< core::pose::Pose >& population){
	vector<Real> population_Cscore;
	
	for (Size i = 0; i < NP_; ++i)
		population_Cscore.push_back( Contact_score( population[i] ) );
	
	return population_Cscore;
}

vector< Real > 
ClassicAbinitio::calculate_population_Dscore(vector< core::pose::Pose >& population){
	vector<Real> population_Dscore;
	
	for (Size i = 0; i < NP_; ++i)
		population_Dscore.push_back( Distance_score( population[i] ) );
	
	return population_Dscore;
}


vector< Real > 
ClassicAbinitio::calculate_population_Rmsd(vector< core::pose::Pose >& population){
	vector<Real> population_Rmsd;
	
	for (Size i = 0; i < NP_; ++i)
		population_Rmsd.push_back( core::scoring::CA_rmsd( population[i], nativePose ) );
	
	return population_Rmsd;
}

vector< Real > 
ClassicAbinitio::calculate_population_SSscore(vector< core::pose::Pose >& population){
	vector<Real> population_SSscore;
	
	for (Size i = 0; i < NP_; ++i)
		population_SSscore.push_back( SS_score( population[i] ) );
	
	return population_SSscore;
}


core::pose::Pose 
ClassicAbinitio::generate_a_new_pose(){
	core::pose::Pose new_pose( initPose );
	
	bool success(true);
	if ( !bSkipStage1_ ){
		if ( !prepare_stage1( new_pose ) ) {
			cout << "ERROR!!!\n" << "ERROR : prepare_stage1 failed!!!\n" << "ERROR!!!" << endl;
			exit(0);
		}
	//	stage3_to_stage4 = false;
		success = do_stage1_cycles( new_pose );
		recover_low( new_pose, STAGE_1 );
		
		if ( !success ) {
			cout << "ERROR!!!\n" << "ERROR : do_stage1_cycles failed!!!\n" << "ERROR!!!" << endl;
			exit(0);
		}
	}
	
	if ( !bSkipStage2_ ){
		if ( !prepare_stage2( new_pose ) ) {
			cout << "ERROR!!!\n" << "ERROR : prepare_stage2 failed!!!\n" << "ERROR!!!" << endl;
			exit(0);
		}
	//	stage3_to_stage4 = true;
		success = do_stage2_cycles( new_pose );
		recover_low( new_pose, STAGE_2 );
		
		if ( !success ) {
			cout << "ERROR!!!\n" << "ERROR : do_stage2_cycles failed!!!\n" << "ERROR!!!" << endl;
			exit(0);
		}
	}
	
	return new_pose;
}

void ClassicAbinitio::SS_search(vector< core::pose::Pose >& population, vector<Real>& population_energy, vector<Real>& population_SSscore){
	cout << "#=================== SS_search =====================" << endl;
	cout << SS << "   " << SSscore_full << endl;
	for ( Size i = 0; i < NP_; ++i){
		for ( Size k = 0; k < 4000; ++k){
			core::pose::Pose trial_pose( population[i] );
			FragAssem_->apply( trial_pose );
			Real trialEnergy = (*score_stage2_)(trial_pose);
			if ( boltzmann_accept(population_energy[i], trialEnergy) ){
				Real trial_SSsccore = SS_score( trial_pose );
				if ( boltzmann_accept(population_SSscore[i], trial_SSsccore, 1) ){
					population[i] = trial_pose;
					population_energy[i] = trialEnergy;
					population_SSscore[i] = trial_SSsccore;
				}
			}
		}
		cout << population[i].secstruct() << "   " << population_SSscore[i] << "   " << population_SSscore[i]/SSscore_full << endl;
	}
}

/*
void
ClassicAbinitio::calculate_distance(core::pose::Pose &pose ){

	vector<vector<Real> > Distance_fragment;
	vector<Real> w_Distance_fragment;
	Real distance;
	core::pose::Pose trial_pose = pose;
	Size c = 0;
	
	//读取片段库的二面角并转化为distance输出

	fragmet_assembly_mover large_frag_assembly;
	large_frag_assembly.read_fragment_library("aat000_09_05.200_v1_3"); 
	Size protein_lenth = pose.total_residue();
	
	for (Size i = 1; i <= pose.total_residue() - 8; ++i){  	//遍历序列位置
		
		core::pose::Pose trial_pose = pose;
		
		vector<Real> w_Distance_fragment;
		for (Size j = 0; j < 200; ++j){	
			
			large_frag_assembly.specific2_position_fragement_insertion(trial_pose, i, j);
			
			if ( trial_pose.residue(i).name3() == "GLY" || trial_pose.residue(i+8).name3() == "GLY" ){
				numeric::xyzVector<Real> CA1 = trial_pose.residue(i).xyz("CA");
				numeric::xyzVector<Real> CA2 = trial_pose.residue(i+8).xyz("CA");
				distance = sqrt((CA2.x() - CA1.x())*(CA2.x() - CA1.x()) + (CA2.y() - CA1.y())*(CA2.y() - CA1.y()) + (CA2.z() - CA1.z())*(CA2.z() - CA1.z()) );
			}else{
				numeric::xyzVector<Real> CB1 = trial_pose.residue(i).xyz("CB");
				numeric::xyzVector<Real> CB2 = trial_pose.residue(i+8).xyz("CB");
				distance = sqrt((CB2.x() - CB1.x())*(CB2.x() - CB1.x()) + (CB2.y() - CB1.y())*(CB2.y() - CB1.y()) + (CB2.z() - CB1.z())*(CB2.z() - CB1.z()) );
			}
			
			w_Distance_fragment.push_back(distance); 
			
		}		
		Distance_fragment.push_back( w_Distance_fragment );	
	}
	
	
	vector<Real>pd;
	string line1,line2;
	Size res1,res2;
	Real lines1,lines2,lines3,dis,pr,mean,var;
	Real two_pi = sqrt(2 * 3.1415925);
	Real sigma = 1.5;
	stringstream s;	
	
	ofstream window;   //输出片段概率
	window.open("./output_files/window.txt");
	
	//读取预测的距离分布
	vector<Real>dist;
	vector<vector<Real>>si;
	ifstream sigma_p ("./alldist.txt");  // 读取trrosetta预测的distance值
	while (getline(sigma_p, line1)){
		istringstream lineData( line1 );	
		lineData >> res1 >> res2 >> dis >> pr >> mean >> var;
		if (res2 - res1 == 8){
			vector<Real>num;
			num.push_back(res1);
			num.push_back(var);	
			si.push_back(num);
			dist.push_back(dis);
			resnum.push_back(res1);
		}
		if (res1 - res2 == 8){
			vector<Real>num;
			num.push_back(res2);
			num.push_back(var);	
			si.push_back(num);
			dist.push_back(dis);
			resnum.push_back(res2);
		}
		
	}
	

	
	for (Size i = 0; i < dist.size(); ++i){
		Real sum = 0;
		vector<Real>d1;
		for (Size j = 0; j < Distance_fragment[i].size(); ++j){             /////////////////////////////判断是否是有预测的distance
			for(Size k = 0; k < si.size(); ++k){
				if (si[k][0] == j){
					sigma = si[k][1];
				}
			}
			
			Real mu = dist[i];
			Size sum1,Gaussian2;
			Real Gaussian1,Gaussian;
			
			Gaussian = (1/(two_pi * sigma)) *  (exp(-pow(Distance_fragment[i][j] - mu,2)/(2 * pow(sigma,2)))) ;
			
			sum += Gaussian;
			
			d1.push_back(sum);
			
		}
		
		vector<Real>d11;
		for (Size k = 0; k < d1.size(); ++k){
			d11.push_back( d1[k] / d1[199] );
		}
		
		
		for (Size m = 0; m < d11.size(); ++m){
			if (m == 0){
				window << " Window   " << i << "  prob   " << d11[0] ;
			}
			else{
				window << i << " : " << d11[m] - d11[m - 1] << endl;
			}
		}
			
		d2.push_back(d11);
		
	}
	
	// ofstream output;           //输出预测的distance
	// output.open("./output_files/distance.csv");
	// for (Size i = 0; i < Distance_fragment.size(); i++){
		// for (Size j = 0; j < Distance_fragment[i].size(); j++){
				// output << Distance_fragment[i][j] << endl;	
		// }		
	// }

}
*/

int ClassicAbinitio::findElement(vector<Real> v, Real key){
	int len = v.size();

	for(int i=0; i<len; i++){

		if(v[i] == key){
			return i;
		}
	}
	return -1;
}
	/*                    =============        根据预测的distance按高斯分布概率选择片段              ==============           */
void ClassicAbinitio::choose_frag(core::pose::Pose &pose){
	
	Size protein_lenth = pose.total_residue();
	Real index(-1);
	int cout11=0;
	while(index==-1){
			insert_po = numeric::random::rg().random_range(1, protein_lenth - 8);
			index=findElement(resnum,insert_po);
			cout11++;
			if(cout11>100)break;
	}
	

//	Size insert_p = (int)(d2[insert_po-1][199]);
	Real insert_fragment( numeric::random::rg().uniform());
	
	          //根据偏差选取片段

	for(Size j = 0; j < d2[index].size(); ++j){
		if ( j == 199 || ( d2[index][j] <= insert_fragment && insert_fragment < d2[index][j+1] ) ){
	
			ind = j;
			
			break;
		}	
		
	}
	// cout << "(" << insert_po << "," << ind << ")   ";
	
	
	///  提取预测的distance中序列间隔等于8的残基部分
	// read_distance();
	
	// ofstream pd;
	// pd.open("./output_files/predict_distance.txt");
	// for (Size k = 0; k < Distance_splines.size(); ++k){
		// // cout << Distance_splines[k] << endl;
		// pd << Distance_splines[k] << endl;
	// }
	
// }
}


void 
ClassicAbinitio::Output_last_generation(vector< core::pose::Pose >& population, vector< Real >& population_Dscore, char c){
    

	Size min_D_index( 0 );
	Size min_RMSD_index( 0 );
	Real min_D( 100000);
	Real min_RMSD( 100000);
	string name= "last_generation";
	string name2= "model_with_min_Energy";
	string name3= "model_with_min_RMSD";
	string csv=".csv";  
	string pdb=".pdb";
	ofstream output(name+c+csv);
	//ofstream out("./output_files/last.csv");
	vector<Real> PopulationRmsd( calculate_population_Rmsd(population));
	
	for ( Size p = 0; p < NP_; ++p ){
		if ( population_Dscore[p] < min_D ){
			min_D = population_Dscore[p];
			min_D_index = p;
		}
		 output << p + 1 << "," << setprecision(2) << fixed << population_Dscore[p] << "," 
		 << setprecision(2) << fixed << core::scoring::CA_rmsd( population[p], nativePose ) <<  endl; 
		
		if ( PopulationRmsd[p] < min_RMSD ){
			min_RMSD = PopulationRmsd[p];
			min_RMSD_index = p;
		}
    }
	cout << "min_RMSD = " << core::scoring::CA_rmsd( population[min_D_index], nativePose ) << endl;
		
	population[min_D_index].dump_pdb(name2+c+pdb);
	population[min_RMSD_index].dump_pdb(name3+c+pdb);
}

void
ClassicAbinitio::norm(Real &angle){
	if(angle>=180){
		angle=angle-360;
	}
	if (angle<=-180)
	{
		angle=360+angle;/* code */
	}
}


void 
ClassicAbinitio::mutation_rand_1(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
	int rand1, rand2, rand3;
	core::pose::Pose pose1, pose2, pose3;
	double scale_F = 0.1;
	Real posemutation_phi(0);
	Real posemutation_psi(0);

	do rand1 = numeric::random::rg().random_range(0, pop_size - 1); while(rand1 == target_id);
	do rand2 = numeric::random::rg().random_range(0, pop_size - 1); while(rand2 == target_id || rand2 == rand1);
	do rand3 = numeric::random::rg().random_range(0, pop_size - 1); while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);

	pose1 = mypose[rand1];
	pose2 = mypose[rand2];
	pose3 = mypose[rand3];
	posemutation = mypose[target_id];
	
	//std::cout << setw(4) << "random select pose: target=" << target_id << ", rand1=" << rand1 << ", rand2=" << rand2 << ", rand3=" << rand3 << std::endl;
	
	for(int k = 1; k <= seq_length; k++)
	{
		if( numeric::random::rg().uniform() > 0.5 ) 
			continue;

		posemutation_phi = pose1.phi(k) + scale_F * ( pose2.phi(k) - pose3.phi(k) );
		posemutation_psi = pose1.psi(k) + scale_F * ( pose2.psi(k) - pose3.psi(k) );

		norm( posemutation_phi );
		norm( posemutation_psi );

		posemutation.set_phi( k, posemutation_phi );
		posemutation.set_psi( k, posemutation_psi );
	}
}

void 
ClassicAbinitio::mutation_current2rand_1(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
	int rand1, rand2, rand3;
	core::pose::Pose pose1, pose2, pose3, pose_current;
	double scale_F = 0.1;
	Real posemutation_phi(0);
	Real posemutation_psi(0);

	do rand1 = numeric::random::rg().random_range(0, pop_size - 1); while(rand1 == target_id);
	do rand2 = numeric::random::rg().random_range(0, pop_size - 1); while(rand2 == target_id || rand2 == rand1);
	do rand3 = numeric::random::rg().random_range(0, pop_size - 1); while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);

	pose1 = mypose[rand1];
	pose2 = mypose[rand2];
	pose3 = mypose[rand3];
	pose_current = mypose[target_id];
	
	//std::cout << setw(4) << "random select pose: target=" << target_id << ", rand1=" << rand1 << ", rand2=" << rand2 << ", rand3=" << rand3 << std::endl;
	
	for(int k = 1; k <= seq_length; k++)
	{
		if( numeric::random::rg().uniform() > 0.5 ) 
			continue;

		posemutation_phi = pose_current.phi(k) + scale_F * ( pose1.phi(k) - pose_current.phi(k) ) + scale_F * ( pose2.phi(k) - pose3.phi(k) );
		posemutation_psi = pose_current.psi(k) + scale_F * ( pose1.psi(k) - pose_current.psi(k) ) + scale_F * ( pose2.psi(k) - pose3.psi(k) );

		norm( posemutation_phi );
		norm( posemutation_psi );

		posemutation.set_phi( k, posemutation_phi );
		posemutation.set_psi( k, posemutation_psi );
	}
}


void 
ClassicAbinitio::mutation_best_1(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
	int rand2=0, rand3=0, best=0;
	core::pose::Pose pose2, pose3, pose_best;
	double scale_F = 0.5;
	Real posemutation_phi(0);
	Real posemutation_psi(0);

	do rand2 = numeric::random::rg().random_range(0, pop_size - 1); while(rand2 == target_id);
	do rand3 = numeric::random::rg().random_range(0, pop_size - 1); while(rand3 == target_id || rand3 == rand2);
	
	double min = PopulationEnergy[0];
	for (int i = 1; i < pop_size; i++)
	{
		if(min > PopulationEnergy[i])
		{
			min = PopulationEnergy[i];
			best = i;
		}
	}

	pose2 = mypose[rand2];
	pose3 = mypose[rand3];
	pose_best = mypose[best];
	posemutation = mypose[target_id];
	
	//std::cout << setw(4) << "random select pose: target=" << target_id << ", rand1=" << rand1 << ", rand2=" << rand2 << ", rand3=" << rand3 << std::endl;
	
	for(int k = 1; k <= seq_length; k++)
	{
		if( numeric::random::rg().uniform() > 0.5 ) 
			continue;

		posemutation_phi = pose_best.phi(k) + scale_F * ( pose2.phi(k) - pose3.phi(k) );
		posemutation_psi = pose_best.psi(k) + scale_F * ( pose2.psi(k) - pose3.psi(k) );

		norm( posemutation_phi );
		norm( posemutation_psi );

		posemutation.set_phi( k, posemutation_phi );
		posemutation.set_psi( k, posemutation_psi );
	}
}

void 
ClassicAbinitio::convex_mutation1(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
	int rand1, rand2, rand3, rand_fra1, rand_fra2;
	pose::Pose pose1, pose2, pose3;
  
	do rand1 = rand()%pop_size; while(rand1 == target_id);
	do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
	do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  
	rand_fra1 = rand()%(seq_length-3)+1;
	do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  
	pose1 = mypose[rand1];
	pose2 = mypose[rand2];
	pose3 = mypose[rand3];
	for(int k = 0 ; k<3 ; k++)
	{
		pose3.set_phi(rand_fra1 + k,pose1.phi(rand_fra1 + k));
		pose3.set_psi(rand_fra1 + k,pose1.psi(rand_fra1 + k));
	}

	for(int k = 0 ; k<3 ; k++)
	{
		pose3.set_phi(rand_fra2 + k,pose2.phi(rand_fra2 + k));
		pose3.set_psi(rand_fra2 + k,pose2.psi(rand_fra2 + k));  
	}

	posemutation = pose3;
}


void 
ClassicAbinitio::convex_mutation2(vector< core::pose::Pose >& mypose,vector<Real>& population_Cscore,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id,Real rate)
{
  int rand1, rand2, rand3, rand_fra1, rand_fra2,rand_start,rand_end,tmp;
  pose::Pose pose1, pose2, pose3;
  // do rand_start=rand()%seq_length;while(rand_start==0);
  // do rand_end = rand()%seq_length;while(abs(rand_start-rand_end)>=20||rand_end==0);
  // if(rand_start>rand_end){
  // 	tmp=rand_start;
  // 	rand_start=rand_end;
  // 	rand_end=tmp;
  // }
  int i=0;
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1||((population_Cscore[rand1]>population_Cscore[rand2])&&(++i<100)));
 // do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  Real aCR(rand()%101/(double)100);
  pose1 = mypose[rand1];
  pose2 = mypose[rand2];
  pose3 = mypose[target_id];
  Real pose3_phi(0);
  Real pose3_psi(0);
  for(int k = 1 ; k<=seq_length ; k++){
  	aCR=rand()%101/(double)100;
  	if(aCR>0.5)continue;
  	pose3_phi=pose3.phi(k)+0.1*(pose1.phi(k)-pose2.phi(k));
  	pose3_psi=pose3.psi(k)+0.1*(pose1.psi(k)-pose2.psi(k));
  	norm(pose3_phi);
  	norm(pose3_psi);

    pose3.set_phi(k ,pose3_phi);
    pose3.set_psi(k ,pose3_psi);
  }
 

  posemutation = pose3;
}


bool ClassicAbinitio::Distance_accept(Real target_Dscore, Real trial_Dscore){
  	if( trial_Dscore < target_Dscore ){
    	//std::cout<<"          accept          "<<std::endl;
    	return true;
  	}
  	else
    	return boltzmann_accept(target_Dscore, trial_Dscore, KTc_);
}
void 
ClassicAbinitio::convex_crossover(core::pose::Pose& target, core::pose::Pose& posemutation, int seq_length){
  int cr_pos;
  cr_pos = rand()%(seq_length-3)+1;
    
  for(int k = 0 ; k<3 ; k++){
    posemutation.set_phi(cr_pos + k, target.phi(cr_pos + k));
    posemutation.set_psi(cr_pos + k, target.psi(cr_pos + k));
  }
}
bool 
ClassicAbinitio::evaluation(core::pose::Pose& pre_pose,core::pose::Pose& pro_pose,Real rate){

	Real aCR(0);
	Real pre_score( 0 );
	Real pro_score( 0 );
	Real pre_distance( 0 );
	Real pro_distance( 0 );
	for (Size i = 0; i < Contacts.size(); ++i){
		// aCR= rand()%101/(double)100;
		// if(aCR>rate)continue;


		Size residue1( Contacts[i].residue1 );
		Size residue2( Contacts[i].residue2 );
		Real confidence ( Contacts[i].confidence );

		//6距离以内不要
		if((residue2-residue1)<6)continue;

		
		if ( pre_pose.residue(residue1).name3() == "GLY" || pre_pose.residue(residue2).name3() == "GLY" ){
			numeric::xyzVector<Real> pre_CA1 = pre_pose.residue(residue1).xyz("CA");
			numeric::xyzVector<Real> pre_CA2 = pre_pose.residue(residue2).xyz("CA");
			pre_distance = pre_CA1.distance(pre_CA2);

			numeric::xyzVector<Real> pro_CA1 = pro_pose.residue(residue1).xyz("CA");
			numeric::xyzVector<Real> pro_CA2 = pro_pose.residue(residue2).xyz("CA");
			pro_distance = pro_CA1.distance(pro_CA2);
		}else{
			numeric::xyzVector<Real> pro_CB1 = pro_pose.residue(residue1).xyz("CB");
			numeric::xyzVector<Real> pro_CB2 = pro_pose.residue(residue2).xyz("CB");
			pro_distance = pro_CB1.distance(pro_CB2);

			numeric::xyzVector<Real> pre_CA1 = pre_pose.residue(residue1).xyz("CB");
			numeric::xyzVector<Real> pre_CA2 = pre_pose.residue(residue2).xyz("CB");
			pre_distance = pre_CA1.distance(pre_CA2);
		}
		//Cscore += log(abs(distance-confidence)+0.1);
		pre_score += pow(pre_distance-confidence,2);
		pro_score += pow(pro_distance-confidence,2);

	}
	return pre_score>pro_score?true:false;
}
void 
ClassicAbinitio::gassu_mutation(vector< core::pose::Pose >& mypose,core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
	/*
  	int rand1, rand2, rand3, rand_fra1, rand_fra2,rand_start,rand_end,tmp;
  	pose::Pose pose1, pose2, pose3;

  	std::default_random_engine gen(rand());
  	std::normal_distribution<double> dis(0,1);

  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  Real aCR(rand()%101/(double)100);
  pose1 = mypose[rand1];
  pose2 = mypose[rand2];
  pose3 = mypose[target_id];
  Real pose3_phi(0);
  Real pose3_psi(0);
  for(int k = 1 ; k<=seq_length ; k++){
  	aCR=rand()%101/(double)100;
  	if(aCR>0.1)continue;
  	pose3_phi=pose3.phi(k)+dis(gen);
  	pose3_psi=pose3.psi(k)+dis(gen);
  	norm(pose3_phi);
  	norm(pose3_psi);

    pose3.set_phi(k ,pose3_phi);
    pose3.set_psi(k ,pose3_psi);
  }
 

  posemutation = pose3;
  */
}


void ClassicAbinitio::apply( pose::Pose & pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;

	Parent::apply( pose );
	if ( option[ OptionKeys::run::dry_run ]() ) return;

///======= @Distance-Guide-Fold ============ @Distance-Guide-Fold =========== @Distance-Guide-Fold ========
///===== @BEGIN =========== @BEGIN ============ @BEGIN ============ @BEGIN ============ @BEGIN ============
	
	cout << "===========================================================================================" << endl;
	cout << "============================== Distance Guide Fold ========================================" << endl;
	cout << "===========================================================================================" << endl;
	///using Rosetta stage1 for populatin initiolization
	bSkipStage1_ = false;
	bSkipStage2_ = true;
	bSkipStage3_ = true;
	bSkipStage4_ = true;

	// pose.dump_pdb("init.pdb");
	// return;
	core::import_pose::pose_from_file( nativePose, option[ in::file::native ]() , core::import_pose::PDB_file);

	//Instantiate a angle search object
	//LJAngleRotationOP LJAR( new LJAngleRotation );
	
	//read and extract parameters
	read_parameters();
	get_parameters(parametersMap);

	proteinLength = pose.total_residue();
	std::cout << "**********************proteinLength=" << proteinLength << std::endl;
	
	//read_contact();
	initPose = pose;
	srand(20); 

	 //initPose.dump_pdb("init.pdb");
	
	system("mkdir -p output_files/SPICKER_data");
	
	///@note Spicker Parameters
	Spicker_All_num = 0;
	Spicker_All_file_num = 1;
	SPICKER_All_Data.open("./output_files/SPICKER_data/spicker.data1");
	
	Ave_acc.open("./Ave_acc.csv");
	
	using namespace core;
	using namespace protocols;
	using namespace protocols::moves;
	using namespace core::scoring;

	using namespace pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::scoring::ScoreFunctionOP score_function_;
	score_function_ = core::scoring::get_score_function(false);
	
///=========================================================================================================	
	///@brief create a initial population and calculate its energy
	cout << "================================ start Initialization =========================================" << endl;
	clock_t starttime = clock();
	//vector<core::pose::Pose> population( generate_population_random(NP_, proteinLength) );
	vector<core::pose::Pose> population( generate_population(NP_) );  //Rosetta stage1
	clock_t endInitime = clock();
	tr << "Time about Initialization: " << (double(endInitime) - starttime )/( CLOCKS_PER_SEC ) << " seconds." << std::endl;
	cout << "============================== Finished Initialization ========================================" << endl;
	
	cout << "==================================== Start add constraints ===================================" << endl;
	for (Size m=0; m<NP_; m++)
	{
		cout << "# add constraints " << m << endl;
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( population[m] );
	}
	vector<Real> PopulationRmsd( calculate_population_Rmsd(population));
	vector<Real> PopulationEnergy( calculate_population_energy(population,*score_function_) );
	Output_last_generation( population, PopulationEnergy, '0');
	
	std::string min_type =  option[ OptionKeys::run::min_type ]();
	core::Real min_tol =  option[ OptionKeys::run::min_tolerance ]();
	core::optimization::MinimizerOptions options( min_type, min_tol, true, false );
	core::kinematics::MoveMap final_mm;
	final_mm.set_chi(  false  );
	final_mm.set_bb(  true  );
	final_mm.set_jump(  false );
	
	for(Size m = 0; m < NP_; m++)
	{
		core::pose::Pose Ctrial( population[m] );
		Real pre_enegy((*score_function_)(Ctrial));
		
		core::optimization::AtomTreeMinimizer().run( Ctrial, final_mm, *score_function_, options );
		Real min_enegy((*score_function_)(Ctrial));
		
		if( min_enegy < pre_enegy )
		{
			population[m] = Ctrial;
			PopulationRmsd[m]=core::scoring::CA_rmsd( Ctrial, nativePose );
			PopulationEnergy[m]=min_enegy;
		}
	}
	Output_last_generation( population, PopulationEnergy, 'O');
	clock_t endtime_minmover = clock();
	tr << "Time about algorithm: " << (double(endtime_minmover) - endInitime )/( CLOCKS_PER_SEC ) << " seconds." << std::endl;
	
	cout << "==================================== start Spicker ============================================" << endl;
	string path = "./output_files/";
	string pdb = ".pdb"; 
	for(Size i=0; i<NP_; i++)
	{
		Output_SPICKER_All_Data(population[i], PopulationEnergy[i]);
		population[i].dump_pdb( path + std::to_string(i) + pdb );
	}
	SPICKER_Demand_All(nativePose);
	
	clock_t endtime = clock();
	tr << "Time about all: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC ) << " seconds." << std::endl;
	

	return;
}// ClassicAbinitio::apply( pose::Pose & pose )


std::string
ClassicAbinitio::get_name() const {
	return "ClassicAbinitio";
}

//@brief return FramgentMover for smooth_small fragment insertions (i.e., stage4 moves)
simple_moves::FragmentMoverOP
ClassicAbinitio::smooth_move_small() {
	return smooth_move_small_;
}

//@brief return FragmentMover for small fragment insertions ( i.e., stage3/4 moves )
simple_moves::FragmentMoverOP
ClassicAbinitio::brute_move_small() {
	return brute_move_small_;
}

//@brief return FragmentMover for large fragment insertions (i.e., stage1/2 moves )
simple_moves::FragmentMoverOP
ClassicAbinitio::brute_move_large() {
	return brute_move_large_;
}

//@brief change the movemap ( is propagated to mover-objects )
//@detail overload if your extension stores additional moves as member variables
void
ClassicAbinitio::set_movemap( core::kinematics::MoveMapCOP mm )
{
	movemap_ = mm;
	if ( smooth_move_small_ ) smooth_move_small_->set_movemap( mm );
	if ( brute_move_small_  ) brute_move_small_ ->set_movemap( mm );
	if ( brute_move_large_  ) brute_move_large_ ->set_movemap( mm );
}

//@brief set new instances of FragmentMovers
void
ClassicAbinitio::set_moves(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small
)
{
	smooth_move_small_ = smooth_move_small;
	brute_move_small_  = brute_move_small;
	brute_move_large_  = brute_move_large;
	update_moves();
}

//@brief returns current movemap
core::kinematics::MoveMapCOP
ClassicAbinitio::movemap() {
	return movemap_;
}

//@detail read cmd_line options and set default versions for many protocol members: trials/moves, score-functions, Monte-Carlo
void ClassicAbinitio::set_defaults( pose::Pose const& pose ) {
	temperature_ = 2.0;
	bSkipStage1_ = false;
	bSkipStage2_ = false;
	bSkipStage3_ = false;
	bSkipStage4_ = false;
	bSkipStage5_ = true; //vats is turned off by default
	set_default_scores();
	set_default_options();
	set_default_mc( pose, *score_stage1_ );
	update_moves();
}

//@detail called to notify about changes in Movers: new movemap or Moverclass
void ClassicAbinitio::update_moves() {
	/* set apply_large_frags_ and
	short_insert_region_
	*/
	/* what about move-map ? It can be set manually for all Fragment_Moves .. */
	// set_move_map();
	set_trials();
}

//@detail create instances of TrialMover for our FragmentMover objects
void ClassicAbinitio::set_trials() {
	// setup loop1
	runtime_assert( brute_move_large_ != nullptr );
	trial_large_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_large_, mc_ ) );
	//trial_large_->set_keep_stats( true );
	trial_large_->keep_stats_type( moves::accept_reject );

	runtime_assert( brute_move_small_ != nullptr );
	trial_small_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_small_, mc_ ) );
	//trial_small_->set_keep_stats( true );
	trial_small_->keep_stats_type( moves::accept_reject );

	runtime_assert( smooth_move_small_ != nullptr );
	smooth_trial_small_ = moves::TrialMoverOP( new moves::TrialMover( smooth_move_small_, mc_ ) );
	//smooth_trial_small_->set_keep_stats( true );
	smooth_trial_small_->keep_stats_type( moves::accept_reject );

	//build trial_pack mover
	moves::SequenceMoverOP combo_small( new moves::SequenceMover() );
	combo_small->add_mover(brute_move_small_);
	combo_small->add_mover(pack_rotamers_);
	trial_small_pack_ = moves::TrialMoverOP( new moves::TrialMover(combo_small, mc_) );
	moves::SequenceMoverOP combo_smooth( new moves::SequenceMover() );
	combo_smooth->add_mover(smooth_move_small_);
	combo_smooth->add_mover(pack_rotamers_);
	smooth_trial_small_pack_ = moves::TrialMoverOP( new moves::TrialMover(combo_smooth, mc_) );
}

//@detail sets Monto-Carlo object to default
void ClassicAbinitio::set_default_mc(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn
) {
	set_mc( moves::MonteCarloOP( new moves::MonteCarlo( pose, scorefxn, temperature_ ) ) );
}

//@detail sets Monto-Carlo object
void ClassicAbinitio::set_mc( moves::MonteCarloOP mc_in ) {
	mc_ = mc_in;
	if ( trial_large_ ) trial_large_->set_mc( mc_ );
	if ( trial_small_ ) trial_small_->set_mc( mc_ );
	if ( smooth_trial_small_ ) smooth_trial_small_->set_mc( mc_ );
}

//@detail override cmd-line setting for "increase_cycling"
void ClassicAbinitio::set_cycles( Real increase_cycles ) {
	stage1_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage2_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage3_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage4_cycles_ = static_cast< int > (4000 * increase_cycles);
	stage5_cycles_ = static_cast< int > (50000* increase_cycles);//vats

	using namespace basic::options;
	if ( option[ OptionKeys::abinitio::only_stage1 ]() ) {
		stage2_cycles_ = 0;
		stage3_cycles_ = 0;
		stage4_cycles_ = 0;
		bSkipStage2_ = bSkipStage3_ = /*bSkipStage3_ =*/ true;  // Was bSkipStage4_ meant? ~Labonte
	}
}

void ClassicAbinitio::set_default_scores() {
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	tr.Debug << "creating standard scoring functions" << std::endl;

	if ( option[ OptionKeys::abinitio::stage1_patch ].user() ) {
		score_stage1_  = ScoreFunctionFactory::create_score_function( "score0", option[ OptionKeys::abinitio::stage1_patch ]() );
	} else {
		score_stage1_  = ScoreFunctionFactory::create_score_function( "score0" );
	}

	if ( option[ OptionKeys::abinitio::stage2_patch ].user() ) {
		score_stage2_  = ScoreFunctionFactory::create_score_function( "score1", option[ OptionKeys::abinitio::stage2_patch ]() );
	} else {
		score_stage2_  = ScoreFunctionFactory::create_score_function( "score1" );
	}

	if ( option[ OptionKeys::abinitio::stage3a_patch ].user() ) {
		score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2", option[ OptionKeys::abinitio::stage3a_patch ]() );
	} else {
		score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2" );
	}

	if ( option[ OptionKeys::abinitio::stage3b_patch ].user() ) {
		score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5", option[ OptionKeys::abinitio::stage3b_patch ]() );
	} else {
		score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5" );
	}

	if ( option[ OptionKeys::abinitio::stage4_patch ].user() ) {
		score_stage4_  = ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage4_patch ]() );
	} else {
		score_stage4_  = ScoreFunctionFactory::create_score_function( "score3" );
	}

	//loading the cenrot score
	score_stage4rot_ = ScoreFunctionFactory::create_score_function( "score4_cenrot_relax" );
	//score_stage4rot_->set_weight(core::scoring::cen_rot_dun, 0.0);
	score_stage4rot_sc_ = ScoreFunctionFactory::create_score_function( "score4_cenrot_repack" );
	//score_stage4rot_sc_->set_weight(core::scoring::cen_rot_dun, 1.0);

	if ( option[ OptionKeys::abinitio::stage5_patch ].user() ) { //vats
		score_stage5_  = ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage5_patch ]() );
	} else {
		score_stage5_  = ScoreFunctionFactory::create_score_function( "score3" );
	}


	if ( option[ OptionKeys::abinitio::override_vdw_all_stages ] ) {
		set_score_weight( scoring::vdw, option[ OptionKeys::abinitio::vdw_weight_stage1 ], ALL_STAGES );
	}
}


/// @brief sets a score weight for all stages of abinitio
void ClassicAbinitio::set_score_weight( scoring::ScoreType type, Real setting, StageID stage ) {
	tr.Debug << "set score weights for ";
	if ( stage == ALL_STAGES ) tr.Debug << "all stages ";
	else tr.Debug << "stage " << (stage <= STAGE_3a ? stage : ( stage-1 ) ) << ( stage == STAGE_3b ? "b " : " " );
	tr.Debug << scoring::name_from_score_type(type) << " " << setting << std::endl;
	if ( score_stage1_  && ( stage == STAGE_1  || stage == ALL_STAGES ) ) score_stage1_ ->set_weight(type, setting);
	if ( score_stage2_  && ( stage == STAGE_2  || stage == ALL_STAGES ) ) score_stage2_ ->set_weight(type, setting);
	if ( score_stage3a_ && ( stage == STAGE_3a || stage == ALL_STAGES ) ) score_stage3a_->set_weight(type, setting);
	if ( score_stage3b_ && ( stage == STAGE_3b || stage == ALL_STAGES ) ) score_stage3b_->set_weight(type, setting);
	if ( score_stage4_  && ( stage == STAGE_4  || stage == ALL_STAGES ) ) score_stage4_ ->set_weight(type, setting);
	if ( score_stage4rot_  && ( stage == STAGE_4  || stage == ALL_STAGES ) ) score_stage4rot_ ->set_weight(type, setting);
	if ( score_stage5_  && ( stage == STAGE_5  || stage == ALL_STAGES ) ) score_stage5_ ->set_weight(type, setting);//vats
}

//@brief currently used score function ( depends on stage )
scoring::ScoreFunction const& ClassicAbinitio::current_scorefxn() const {
	return mc().score_function();
}

//@brief set current scorefunction
void ClassicAbinitio::current_scorefxn( scoring::ScoreFunction const& scorefxn ) {
	mc().score_function( scorefxn );
}

//@brief set individual weight of current scorefunction --- does not change the predefined scores: score_stageX_
void ClassicAbinitio::set_current_weight( core::scoring::ScoreType type, core::Real setting ) {
	scoring::ScoreFunctionOP scorefxn ( mc().score_function().clone() );
	scorefxn->set_weight( type, setting );
	mc().score_function( *scorefxn ); //trigger rescore
}

void ClassicAbinitio::set_default_options() {
	bSkipStage1_ = bSkipStage2_ = bSkipStage3_ = bSkipStage4_ = false;
	bSkipStage5_ = true; //vats turned off by default
	using namespace basic::options;
	just_smooth_cycles_ = option[ OptionKeys::abinitio::smooth_cycles_only ]; // defaults to false
	bQuickTest_ = basic::options::option[ basic::options::OptionKeys::run::test_cycles ]();

	if ( bQuickTest() ) {
		set_cycles( 0.001 );
	} else {
		set_cycles( option[ OptionKeys::abinitio::increase_cycles ] ); // defaults to factor of 1.0
	}

	if ( just_smooth_cycles_ ) {
		bSkipStage1_ = bSkipStage2_ = bSkipStage3_ = bSkipStage5_ = true;
	}
	if ( option[ OptionKeys::abinitio::only_stage1 ] ) {
		bSkipStage2_ = bSkipStage3_ = bSkipStage4_ = bSkipStage5_= true;
	}

	if ( option[ OptionKeys::abinitio::include_stage5 ] ) {
		bSkipStage5_ = false;
	}

	apply_large_frags_   = true;  // apply large frags in phase 2!

	// in rosetta++ switched on in fold_abinitio if contig_size < 30 in pose_abinitio never
	short_insert_region_ = false;  // apply small fragments in phase 2!

	if ( option[ OptionKeys::abinitio::recover_low_in_stages ].user() ) {
		for ( int it : option[ OptionKeys::abinitio::recover_low_in_stages ]() ) {
			if ( it == 1 ) recover_low_stages_.push_back( STAGE_1 );
			else if ( it == 2 ) recover_low_stages_.push_back( STAGE_2 );
			else if ( it == 3 ) {
				recover_low_stages_.push_back( STAGE_3a );
				recover_low_stages_.push_back( STAGE_3b );
			} else if ( it == 4 ) recover_low_stages_.push_back( STAGE_4 );
		}
	} else {
		recover_low_stages_.clear();
		recover_low_stages_.push_back( STAGE_1 );
		recover_low_stages_.push_back( STAGE_2 );
		recover_low_stages_.push_back( STAGE_3a );
		recover_low_stages_.push_back( STAGE_3b );
		recover_low_stages_.push_back( STAGE_4 );
		recover_low_stages_.push_back( STAGE_5 );
	}

	close_chbrk_ = option[ OptionKeys::abinitio::close_chbrk ];

}


/// @brief (helper) functor class which keeps track of old pose for the
/// convergence check in stage3 cycles
/// @detail
/// calls of operator ( pose ) compare the
class hConvergenceCheck;
using hConvergenceCheckOP = utility::pointer::shared_ptr<hConvergenceCheck>;

class hConvergenceCheck : public moves::PoseCondition {
public:
	hConvergenceCheck()= default;
	void reset() { ct_ = 0; bInit_ = false; }
	void set_trials( moves::TrialMoverOP trin ) {
		trials_ = trin;
		runtime_assert( trials_->keep_stats_type() < moves::no_stats );
		last_move_ = 0;
	}
	bool operator() ( const core::pose::Pose & pose ) override;
private:
	pose::Pose very_old_pose_;
	bool bInit_{ false };
	Size ct_ = 0;
	moves::TrialMoverOP trials_;
	Size last_move_;
};

// keep going --> return true
bool hConvergenceCheck::operator() ( const core::pose::Pose & pose ) {
	if ( !bInit_ ) {
		bInit_ = true;
		very_old_pose_ = pose;
		return true;
	}
	runtime_assert( trials_ != nullptr );
	tr.Trace << "TrialCounter in hConvergenceCheck: " << trials_->num_accepts() << std::endl;
	if ( numeric::mod(trials_->num_accepts(),100) != 0 ) return true;
	if ( (Size) trials_->num_accepts() <= last_move_ ) return true;
	last_move_ = trials_->num_accepts();
	// change this later to this: (after we compared with rosetta++ and are happy)
	// if ( numeric::mod(++ct_, 1000) != 0 ) return false; //assumes an approx acceptance rate of 0.1

	// still here? do the check:

	core::Real converge_rms = core::scoring::CA_rmsd( very_old_pose_, pose );
	very_old_pose_ = pose;
	if ( converge_rms >= 3.0 ) {
		return true;
	}
	// if we get here thing is converged stop the While-Loop
	tr.Info << " stop cycles in stage3 due to convergence " << std::endl;
	return false;
}


bool ClassicAbinitio::do_stage1_cycles( pose::Pose &pose ) {
	AllResiduesChanged done( pose, brute_move_large()->insert_map(), *movemap() );
	moves::MoverOP trial( stage1_mover( pose, trial_large() ) );

	// FragmentMoverOP frag_mover = brute_move_large_;
	// fragment::FragmentIO().write("stage1_frags_classic.dat",*frag_mover->fragments());

	Size j;
	for ( j = 1; j <= stage1_cycles(); ++j ) {
		trial->apply( pose ); // apply a large fragment insertion, accept with MC boltzmann probability
		if ( done(pose) ) {
			tr.Info << "Replaced extended chain after " << j << " cycles." << std::endl;
			mc().reset( pose ); // make sure that we keep the final structure
			return true;
		}
	}
	tr.Warning << "extended chain may still remain after " << stage1_cycles() << " cycles!" << std::endl;
	done.show_unmoved( pose, tr.Warning );
	mc().reset( pose ); // make sure that we keep the final structure
	return true;
}

bool ClassicAbinitio::do_stage2_cycles( pose::Pose &pose ) {

	//setup cycle
	moves::SequenceMoverOP cycle( new moves::SequenceMover() );
	if ( apply_large_frags_   ) cycle->add_mover( trial_large_->mover() );
	if ( short_insert_region_ ) cycle->add_mover( trial_small_->mover() );

	Size nr_cycles = stage2_cycles() / ( short_insert_region_ ? 2 : 1 );
	moves::TrialMoverOP trials( new moves::TrialMover( cycle, mc_ptr() ) );
	moves::RepeatMover( stage2_mover( pose, trials ), nr_cycles ).apply(pose);

	//is there a better way to find out how many steps ? for instance how many calls to scoring?
	return true; // as best guess
}

/*! @detail stage3 cycles:
nloop1 : outer iterations
nloop2 : inner iterations
stage3_cycle : trials per inner iteration
every inner iteration we switch between score_stage3a ( default: score2 ) and score_stage3b ( default: score 5 )

prepare_loop_in_stage3() is called before the stage3_cycles() of trials are started.

first outer loop-iteration is done with TrialMover trial_large()
all following iterations with trial_small()

start each iteration with the lowest_score_pose. ( mc->recover_low() -- called in prepare_loop_in_stage3() )

*/
bool ClassicAbinitio::do_stage3_cycles( pose::Pose &pose ) {
	using namespace ObjexxFCL;

	// interlaced score2 / score 5 loops
	// nloops1 and nloops2 could become member-variables and thus changeable from the outside
	int nloop1 = 1;
	int nloop2 = 10; //careful: if you change these the number of structures in the structure store changes.. problem with checkpointing
	// individual checkpoints for each stage3 iteration would be a remedy. ...

	if ( short_insert_region_ ) {
		nloop1 = 2;
		nloop2 = 5;
	}

	hConvergenceCheckOP convergence_checker ( nullptr );
	if ( !option[ basic::options::OptionKeys::abinitio::skip_convergence_check ] ) {
		convergence_checker = hConvergenceCheckOP( new hConvergenceCheck );
	}

	moves::TrialMoverOP trials = trial_large();
	int iteration = 1;
	for ( int lct1 = 1; lct1 <= nloop1; lct1++ ) {
		if ( lct1 > 1 ) trials = trial_small(); //only with short_insert_region!
		for ( int lct2 = 1; lct2 <= nloop2; lct2++, iteration++  ) {
			tr.Debug << "Loop: " << lct1 << "   " << lct2 << std::endl;

			if ( !prepare_loop_in_stage3( pose, iteration, nloop1*nloop2 ) ) return false;

			if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2),
					false /*fullatom */, true /*fold tree */ ) ) {


				tr.Debug << "  Score stage3 loop iteration " << lct1 << " " << lct2 << std::endl;
				if ( convergence_checker ) {
					moves::TrialMoverOP stage3_trials = stage3_mover( pose, lct1, lct2, trials );
					convergence_checker->set_trials( stage3_trials ); //can be removed late
					moves::WhileMover( stage3_trials, stage3_cycles(), convergence_checker ).apply( pose );
				} else {    //no convergence check -> no WhileMover
					moves::RepeatMover( stage3_mover( pose, lct1, lct2, trials ), stage3_cycles() ).apply( pose );
				}

				if ( numeric::mod( (int)iteration, 2 ) == 0 || iteration > 7 ) recover_low( pose, STAGE_3a );
				recover_low( pose, STAGE_3b );

				get_checkpoints().checkpoint( pose, get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2), true /*fold tree */ );
			}//recover_checkpoint
			get_checkpoints().debug( get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2), current_scorefxn()( pose ) );

			//   structure_store().push_back( mc_->lowest_score_pose() );
		} // loop 2
	} // loop 1
	return true;
}


// interlaced score2 / score 5 loops
/*! @detail stage4 cycles:
nloop_stage4: iterations
stage4_cycle : trials per  iteration

first iteration: use trial_small()
following iterations: use trial_smooth()
only trial_smooth() if just_smooth_cycles==true

prepare_loop_in_stage4() is called each time before the stage4_cycles_ of trials are started.

start each iteration with the lowest_score_pose. ( mc->recover_low()  in prepare_loop_in_stage4()  )

*/
bool ClassicAbinitio::do_stage4_cycles( pose::Pose &pose ) {
	Size nloop_stage4 = 3;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[corrections::score::cenrot]() ) nloop_stage4=2;

	for ( Size kk = 1; kk <= nloop_stage4; ++kk ) {
		tr.Debug << "prepare ..." << std::endl ;
		if ( !prepare_loop_in_stage4( pose, kk, nloop_stage4 ) ) return false;

		if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk), false /* fullatom */, true /* fold_tree */ ) ) {
			moves::TrialMoverOP trials;
			if ( kk == 1 && !just_smooth_cycles_ ) {
				trials = trial_small();
			} else {
				tr.Debug << "switch to smooth moves" << std::endl;
				trials = trial_smooth();
			}

			tr.Debug << "start " << stage4_cycles() << " cycles" << std::endl;
			moves::RepeatMover( stage4_mover( pose, kk, trials ), stage4_cycles() ).apply(pose);
			tr.Debug << "finished" << std::endl;
			recover_low( pose, STAGE_4 );

			get_checkpoints().checkpoint( pose, get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk), true /*fold tree */ );
		}
		get_checkpoints().debug( get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk),  current_scorefxn()( pose ) );

		//don't store last structure since it will be exactly the same as the final structure delivered back via apply
		//  if( kk < nloop_stage4 ) // <-- this line was missing although the comment above was existant.
		//   structure_store().push_back( mc_->lowest_score_pose() );
	}  // loop kk

	if ( option[corrections::score::cenrot] ) {
		//switch to cenrot model
		tr.Debug << "switching to cenrot model ..." << std::endl;
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot(chemical::CENTROID_ROT);
		to_cenrot.apply(pose);

		//init pose
		(*score_stage4rot_)( pose );
		pack_rotamers_->score_function(score_stage4rot_sc_);
		pack_rotamers_->apply(pose);

		mc_->reset(pose);
		replace_scorefxn( pose, STAGE_4rot, 0 );
		//mc_->set_temperature(1.0);
		//mc_->set_autotemp(true, 1.0);

		//debug
		//tr.Debug << "starting_energy: " << (*score_stage4rot_)( pose ) << std::endl;
		//tr.Debug << "starting_temperature: " << mc_->temperature() << std::endl;

		for ( Size rloop=1; rloop<=3; rloop++ ) {
			//change vdw weight
			switch (rloop) {
			case 1 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)/9.0);
				break;
			case 2 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)*3.0);
				break;
			case 3 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)*3.0);
				break;
			}

			//stage4rot
			//for (Size iii=1; iii<=100; iii++){
			//pose::Pose startP = pose;
			//tr << "temperature: " << mc_->temperature() << std::endl;
			moves::RepeatMover( stage4rot_mover( pose, rloop, trial_smooth() ), stage4_cycles()/100 ).apply(pose);
			//tr << "delta_rms: " << core::scoring::CA_rmsd( startP, pose ) << std::endl;
			//}
		}
	}

	return true;
}

bool ClassicAbinitio::do_stage5_cycles( pose::Pose &pose ) {//vats

	Size nmoves = 1;
	core::kinematics::MoveMapOP mm_temp( new core::kinematics::MoveMap( *movemap() ) );
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( mm_temp, temperature_, nmoves) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 5.0 );

	moves::TrialMoverOP trials( new moves::TrialMover( small_mover, mc_ptr() ) );
	moves::RepeatMover( stage5_mover( pose, trials ), stage5_cycles() ).apply( pose );

	// moves::MoverOP trial( stage5_mover( pose, small_mover ) );
	// Size j;
	// for( j = 1; j <= stage5_cycles(); ++j ) {
	//  trial->apply( pose );
	// }
	mc().reset( pose );
	return true;

}


moves::TrialMoverOP
ClassicAbinitio::stage1_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}


moves::TrialMoverOP
ClassicAbinitio::stage2_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}


moves::TrialMoverOP
ClassicAbinitio::stage3_mover( pose::Pose &, int, int, moves::TrialMoverOP trials ) {
	return trials;
}

moves::TrialMoverOP
ClassicAbinitio::stage4_mover( pose::Pose &, int, moves::TrialMoverOP trials ) {
	return trials;
}

moves::TrialMoverOP
ClassicAbinitio::stage4rot_mover( pose::Pose &, int, moves::TrialMoverOP trials ) {
	if ( trials == trial_small_ ) {
		return trial_small_pack_;
	} else {
		return smooth_trial_small_pack_;
	}
}

moves::TrialMoverOP //vats
ClassicAbinitio::stage5_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}

void ClassicAbinitio::recover_low( core::pose::Pose& pose, StageID stage ){
	if ( contains_stageid( recover_low_stages_, stage ) ) {
		mc_->recover_low( pose );
	}
}

// anything you want to have done before the stages ?
void ClassicAbinitio::replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real /*intra_stage_progress */ ) {
	// must assume that the current pose is the one to be accepted into the next stage! (this change was necessary for
	// checkpointing to work correctly.

	//intra_stage_progress = intra_stage_progress;
	if ( score_stage1_  && ( stage == STAGE_1 ) ) current_scorefxn( *score_stage1_ );
	if ( score_stage2_  && ( stage == STAGE_2 ) ) current_scorefxn( *score_stage2_ );
	if ( score_stage3a_ && ( stage == STAGE_3a) ) current_scorefxn( *score_stage3a_ );
	if ( score_stage3b_ && ( stage == STAGE_3b) ) current_scorefxn( *score_stage3b_ );
	if ( score_stage4_  && ( stage == STAGE_4 ) ) current_scorefxn( *score_stage4_ );
	if ( score_stage4rot_  && ( stage == STAGE_4rot ) ) current_scorefxn( *score_stage4rot_ );
	if ( score_stage5_  && ( stage == STAGE_5 ) ) current_scorefxn( *score_stage5_ );//vats
	Real temperature( temperature_ );
	if ( stage == STAGE_5 ) temperature = 0.5;
	mc_->set_autotemp( true, temperature );
	mc_->set_temperature( temperature ); // temperature might have changed due to autotemp..
	mc_->reset( pose );
}


moves::TrialMoverOP ClassicAbinitio::trial_large() {
	return ( apply_large_frags_ ? trial_large_ : trial_small_ );
}

moves::TrialMoverOP ClassicAbinitio::trial_small() {
	return trial_small_;
}

moves::TrialMoverOP ClassicAbinitio::trial_smooth() {
	return smooth_trial_small_;
}

// prepare stage1 sampling
bool ClassicAbinitio::prepare_stage1( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_1, 0.5 );
	mc_->set_autotemp( false, temperature_ );
	// mc_->set_temperature( temperature_ ); already done in replace_scorefxn
	// mc_->reset( pose );
	(*score_stage1_)( pose );
	/// Now handled automatically.  score_stage1_->accumulate_residue_total_energies( pose ); // fix this
	return true;
}

bool ClassicAbinitio::prepare_stage2( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_2, 0.5 );

	(*score_stage2_)(pose);
	/// Now handled automatically.  score_stage2_->accumulate_residue_total_energies( pose );
	return true;
}


bool ClassicAbinitio::prepare_stage3( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_3a, 0 );
	//score for this stage is changed in the do_stage3_cycles explicitly
	if ( option[ templates::change_movemap ].user() && option[ templates::change_movemap ] == 3 ) {
		kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
		new_mm->set_bb( true );
		set_movemap( new_mm ); // --> store it in movemap_ --> original will be reinstated at end of apply()
	}
	return true;
}


bool ClassicAbinitio::prepare_stage4( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_4, 0 );
	(*score_stage4_)( pose );
	/// Now handled automatically.  score_stage4_->accumulate_residue_total_energies( pose ); // fix this

	if ( option[ templates::change_movemap ].user() && option[ templates::change_movemap ] == 4 ) {
		kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
		new_mm->set_bb( true );
		tr.Debug << "option: templates::change_movemap ACTIVE: set_movemap" << std::endl;
		set_movemap( new_mm ); // --> store it in movemap_ --> original will be reinstated at end of apply()
	}
	return true;
}

bool ClassicAbinitio::prepare_stage5( core::pose::Pose &pose ) {//vats
	// temperature_ = 0.5; //this has to be reset to original temperature!!!
	// no special if-statement in replace_scorefxn...OL
	replace_scorefxn( pose, STAGE_5, 0 );
	(*score_stage5_)( pose );
	return true;
}


bool ClassicAbinitio::prepare_loop_in_stage3( core::pose::Pose &pose/*pose*/, Size iteration, Size total ){
	// interlace score2/score5

	Real chbrk_weight_stage_3a = 0;
	Real chbrk_weight_stage_3b = 0;

	if ( numeric::mod( (int)iteration, 2 ) == 0 || iteration > 7 ) {
		Real progress( iteration );
		chbrk_weight_stage_3a = 0.25 * progress;
		tr.Debug << "select score_stage3a..." << std::endl;
		recover_low( pose, STAGE_3a );
		replace_scorefxn( pose, STAGE_3a, 1.0* iteration/total );
	} else {
		Real progress( iteration );
		chbrk_weight_stage_3b = 0.05 * progress;
		tr.Debug << "select score_stage3b..." << std::endl;
		recover_low( pose, STAGE_3b );
		replace_scorefxn( pose, STAGE_3b, 1.0* iteration/total );
	}

	if ( close_chbrk_ ) {

		set_score_weight( scoring::linear_chainbreak, chbrk_weight_stage_3a , STAGE_3a );
		set_score_weight( scoring::linear_chainbreak, chbrk_weight_stage_3b , STAGE_3b );

	}


	return true;
}

bool ClassicAbinitio::prepare_loop_in_stage4( core::pose::Pose &pose, Size iteration, Size total ){
	replace_scorefxn( pose, STAGE_4, 1.0* iteration/total );

	Real chbrk_weight_stage_4 (iteration*0.5+2.5);

	if ( close_chbrk_ ) {
		set_current_weight( scoring::linear_chainbreak, chbrk_weight_stage_4 );
	}

	return true;
}

//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
moves::MonteCarloOP
ClassicAbinitio::mc_ptr() {
	return mc_;
}


void ClassicAbinitio::output_debug_structure( core::pose::Pose & pose, std::string prefix ) {
	using namespace core::io::silent;

	mc().score_function()( pose );
	Parent::output_debug_structure( pose, prefix );

	if ( option[ basic::options::OptionKeys::abinitio::explicit_pdb_debug ]() ) {
		pose.dump_pdb( prefix + get_current_tag() + ".pdb" );
	}

	if ( option[ basic::options::OptionKeys::abinitio::log_frags ].user() ) {
		std::string filename = prefix + "_" + get_current_tag() + "_" + std::string( option[ basic::options::OptionKeys::abinitio::log_frags ]() );
		utility::io::ozstream output( filename );
		auto& log_frag = dynamic_cast< simple_moves::LoggedFragmentMover& > (*brute_move_large_);
		log_frag.show( output );
		log_frag.clear();
	}

} // ClassicAbinitio::output_debug_structure( core::pose::Pose & pose, std::string prefix )

} //abinitio
} //protocols


