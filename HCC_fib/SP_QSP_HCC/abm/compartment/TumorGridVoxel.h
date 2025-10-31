#pragma once

#include "SP_QSP_shared/ABM_Base/AgentGridVoxel.h"
#include "../../core/AgentEnum.h"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
//#include <boost/serialization/export.hpp>

//class CancerCell;

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

class TumorGridVoxel : public AgentGridVoxel
{
public:
	TumorGridVoxel(){};
	TumorGridVoxel(int x, int y, int z);
	virtual ~TumorGridVoxel();
	//! if this voxel is open to Agent of given type
	virtual bool isOpenToType(AgentType t)const;
	//! get distance to (0,0,0)
	double getDistToOrigin(void)const;
	//std::string getLoc(void) const { return _str_loc; };
	//! Get vasculature density
	double getVasDensity() const { return _vas_density; };
	void setVasDensity(double new_density);

	//! Get ECM density
	double getECMDensity() const { return _ECM_density; };
	void setECMDensity(double new_density);

	//! Get ECM density
	double getFibField() const { return _Fib_field; };
	void setFibField(double new_field);

protected:
	double _distToOrigin;
	//std::string _str_loc;
	virtual int getVoxelContent()const;
private:
	friend class boost::serialization::access;
	double _vas_density;
	double _ECM_density;
	double _Fib_field;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);
};

// derived class objects serialized through pointers need to be exported
//BOOST_CLASS_EXPORT_KEY(TumorGridVoxel);

template<class Archive>
inline void TumorGridVoxel::serialize(Archive & ar, const unsigned int version){
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AgentGridVoxel);
	ar & BOOST_SERIALIZATION_NVP(_distToOrigin);
	ar & BOOST_SERIALIZATION_NVP(_vas_density);
	ar & BOOST_SERIALIZATION_NVP(_ECM_density);
}

};
};
