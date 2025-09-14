#pragma once

#include "SP_QSP_shared/ABM_Base/Coord3D.h"
#include "SP_QSP_shared/ABM_Base/BaseAgent.h"
#include "SP_QSP_shared/ABM_Base/RNG.h"
#include "../../core/AgentEnum.h"
#include <boost/serialization/nvp.hpp>
#include <cmath>

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

/*! This class contains functions that
	populate voxels randomly.
	Can be used during grid initialization or when new 
	voxels enters grid.
*/
class VoxelContentGen
{
private:
typedef BaseAgent::AgentType AgentType;
typedef BaseAgent::AgentState AgentState;

public:
	VoxelContentGen();
	~VoxelContentGen();
	//! populate one voxel with cell agents
	bool get_type_state(const Coord3D& c, RNG& rng, 
		AgentType& type, AgentState& state, int& div)const;
	//! setup variables, equilibrium IC
	void setup(bool stationary, double cancer_prob, 
		int xlim, int ylim, int zlim); 
	//! setup variables, mandate IC
	void setup(bool stationary, int xlim, int ylim, int zlim, int x0, int y0, int z0);
	//circular setup of tumor
	void setup(bool stationary, double center_x, double center_y, double center_z, double radius);
	// setup cell cdf
	void set_cell_cdf(bool stationary);
private:
	friend class boost::serialization::access;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);

	//! same density everywhere
	bool _stationary;
	//! min/max xyz for initial population
	int _x_min;
	int _y_min;
	int _z_min;
	int _x_max;
	int _y_max;
	int _z_max;
	//circular tumor initialization
	bool _circular;
	double _radius;
	int _center_x;
	int _center_y;
	int _center_z;
	//! Position of fibroblast
	bool _include_fib;
	int _x_fib;
	int _y_fib;
	int _z_fib;
	//! probability of subtypes and division numbers
	std::vector<double> _celltype_cdf;
};

template<class Archive>
inline void VoxelContentGen::serialize(Archive & ar, const unsigned int /* version */){
	ar & BOOST_SERIALIZATION_NVP(_stationary);
	ar & BOOST_SERIALIZATION_NVP(_x_min);
	ar & BOOST_SERIALIZATION_NVP(_y_min);
	ar & BOOST_SERIALIZATION_NVP(_z_min);
	ar & BOOST_SERIALIZATION_NVP(_x_max);
	ar & BOOST_SERIALIZATION_NVP(_y_max);
	ar & BOOST_SERIALIZATION_NVP(_z_max);
	ar & BOOST_SERIALIZATION_NVP(_celltype_cdf);
}

};// end of namespace
};