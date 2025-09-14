#pragma once

#include "Cell_Tumor.h"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

class Mac:
	public Cell_Tumor
{
public:
	Mac(){};
	Mac(SpatialCompartment* c );
	Mac(const Mac& c);
	virtual ~Mac();

	virtual CellAgent* createCellCopy() const { return new Mac(*this); };

	//! print cancer cell information
	virtual std::string toString() const;

	//! step function for cancer cell
	// virtual void agentStep(double t, double dt, AgentStep & as);
	//! get cell agent type, in this case cancer cell.
	virtual AgentType getType() const { return AgentTypeEnum::CELL_TYPE_MAC; };
	//virtual void agentStep
	virtual bool agent_movement_step(double t, double dt, Coord& c);
	virtual bool agent_state_step(double t, double dt, Coord& c);
	virtual void agent_state_scan(void);

	//! PD1_PDL1 bond in synapse
	static double get_PD1_PDL1(double PDL1, double Nivo);
	//! get suppression from PD1_PDL1 bond
	static double get_H(double bond, double n, double k50);
	//! calculate probability of Cancer cell killed by Macrophage
	static double get_kill_prob(double supp, double q);

	//! move all cytokine source
	void move_all_source_sink(void) const;
	//! remove all source sink
	void remove_all_source_sink(void);
	//! set macrophage to M1
	void setM1();
	//! set macrophage to M2 
	void setM2();
	//return source pointer
	BioFVMSinkSource*& get_source_IFNg(void) { return _source_IFNg; };
	BioFVMSinkSource*& get_source_IL_12(void) { return _source_IL_12; };
	BioFVMSinkSource*& get_source_TGFB(void) { return _source_TGFB; };
	BioFVMSinkSource*& get_source_IL_10(void) { return _source_IL_10; };
	BioFVMSinkSource*& get_source_VEGFA(void) { return _source_VEGFA; };
	BioFVMSinkSource*& get_sink_CCL2(void) { return _sink_CCL2; };

	int getMacLife();
private:
	//! cytokine source
	BioFVMSinkSource* _source_IFNg;
	BioFVMSinkSource* _source_IL_12;
	BioFVMSinkSource* _source_TGFB;
	BioFVMSinkSource* _source_IL_10;
	BioFVMSinkSource* _source_VEGFA;
	BioFVMSinkSource* _sink_CCL2;

	//! neighbor CD
	double _max_neighbor_PDL1;

	friend class boost::serialization::access;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);
};

//BOOST_CLASS_EXPORT_KEY(Mac)

template<class Archive>
inline void Mac::serialize(Archive & ar, const unsigned int /* version */){
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Cell_Tumor);
	ar& BOOST_SERIALIZATION_NVP(_source_IFNg);
	ar& BOOST_SERIALIZATION_NVP(_source_IL_12);
	ar& BOOST_SERIALIZATION_NVP(_source_TGFB);
	ar& BOOST_SERIALIZATION_NVP(_source_IL_10);
	ar& BOOST_SERIALIZATION_NVP(_source_VEGFA);
	ar& BOOST_SERIALIZATION_NVP(_sink_CCL2);
}
};
};
