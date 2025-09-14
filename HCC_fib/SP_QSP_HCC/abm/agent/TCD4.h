#pragma once

#include "Cell_Tumor.h"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include "../../pde/DiffuseGrid.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

class TCD4:
	public Cell_Tumor
{
public:
	TCD4(){};
	TCD4(SpatialCompartment* c );
	TCD4(const TCD4& c);
	virtual ~TCD4();

	virtual CellAgent* createCellCopy() const { return new TCD4(*this); };

	//! print cancer cell information
	virtual std::string toString() const;

	//! PD1_PDL1 bond in synapse
	static double get_CTLA4_Ipi(double Ipi);

	static int getTCD4Life();
	//void setDead(void);

	//! step function for TCD4 cell
	bool agent_movement_step(double t, double dt, Coord& c);
	bool agent_state_step(double t, double dt, Coord& c);
	//! set CD4 cell to Th (T helper cell)
	void setTh();
	//! set CD4 cell to Treg (T regulatory cell)
	void setTreg();
	//! move sources (IL10)
	void move_all_source_sink(void)const;
	//! remove all sources
	void remove_all_source_sink(void);
	//! TGFB release time remaining, sec
	double _TGFB_release_remain;

	int get_divide_cd_TCD4_exp() const { return _divide_cd_TCD4_exp; };
	//! get cell agent type, in this case cancer cell.
	virtual AgentType getType() const { return AgentTypeEnum::CELL_TYPE_TCD4; };
	//! return source pointer
	// BioFVMSinkSource*& get_source_IFNg(void) { return _source_IFNg; };
	BioFVMSinkSource*& get_source_IL_2(void) { return _source_IL_2; };
	BioFVMSinkSource*& get_source_TGFB(void) { return _source_TGFB; };
	BioFVMSinkSource*& get_source_IL_10(void) { return _source_IL_10; };

private:
	friend class boost::serialization::access;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);
	//! flag that divison is to occur
	bool _divide_flag;
	//! accumulative IL2 exposure, sec*ng/mL 
	double _IL2_exposure;
	//! cool down timer for division
	int _divide_cd_TCD4_exp;
	//! number of remaining division
	int _divide_limit_TCD4_exp;
	// total number of CTLA4 molecules on TCD4 cells
	double _CTLA4;
	//BioFVMSinkSource* _source_IL_10;
	//BioFVMSinkSource* _source_IFNg;
	BioFVMSinkSource* _source_IL_2;
	BioFVMSinkSource* _source_IL_10;
	BioFVMSinkSource* _source_TGFB;
};

//BOOST_CLASS_EXPORT_KEY(TCD4)

template<class Archive>
inline void TCD4::serialize(Archive & ar, const unsigned int /* version */){
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Cell_Tumor);
	ar & BOOST_SERIALIZATION_NVP(_divide_cd_TCD4_exp);
	ar & BOOST_SERIALIZATION_NVP(_divide_limit_TCD4_exp);	
	// ar & BOOST_SERIALIZATION_NVP(_source_IFNg);
	ar & BOOST_SERIALIZATION_NVP(_source_IL_2);
	ar& BOOST_SERIALIZATION_NVP(_source_IL_10);
	ar& BOOST_SERIALIZATION_NVP(_source_TGFB);
}

};
};
