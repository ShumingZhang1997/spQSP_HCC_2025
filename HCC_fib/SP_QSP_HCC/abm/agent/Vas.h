#pragma once

#include "Cell_Tumor.h"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>


namespace SP_QSP_IO {
namespace SP_QSP_HCC {

typedef std::set<Coord3D> CoordSet;
typedef std::vector<Coord3D> CoordVec;
typedef std::vector<Coord3D> ShapeCoords;
typedef std::vector<ShapeCoords> ShapeVec;

class Vas :
	public Cell_Tumor
{
public:
	Vas() {};
	Vas(SpatialCompartment* c);
	Vas(const Vas& c);
	virtual ~Vas();

	virtual CellAgent* createCellCopy() const { return new Vas(*this); };

	//virtual void agentStep
	virtual bool agent_movement_step(double t, double dt, Coord& c);
	virtual bool agent_state_step(double t, double dt, Coord& c);
	virtual void agent_state_scan(void) {};

	//! print cancer cell information
	virtual std::string toString() const;

	//! remove all sources
	void remove_all_source_sink(void);
	//! step function for cancer cell
	//virtual void agentStep(double t, double dt, AgentStep & as);
	//! get cell agent type, in this case cancer cell.
	virtual AgentType getType() const { return AgentTypeEnum::CELL_TYPE_VAS; };

	//! set tip cells
	void set_tip();
	void set_stalk();
	void set_phalanx();
	void make_mature();
	Vas* get_upstream_neighbor() { return _upstream_neighbor; };
	void set_upstream_neighbor(Vas* upstream_neighbor);

	unsigned int get_tip_id() { return _tip_id; };
	void set_tip_id(unsigned int tip_id);
	void set_branch();
	void set_moveDirection();
	// update vasculature density
	BioFVMSinkSource*& get_source_O2(void)  { return _source_O2; };
	BioFVMSinkSource*& get_sink_VEGFA(void) { return _sink_VEGFA; };
private:
	friend class boost::serialization::access;
	//! boost serialization
	template<class Archive>
	void serialize(Archive& ar, const unsigned int /*version*/);
	// limits of rounds of cell division
	int _divide_limit;
	// cool down of cell division
	int _divide_cd;

	Coord3D _moveDirection;
	bool _tumble;
	//whether a cell become a branching cell
	bool _branch;
	Vas* _upstream_neighbor;

	unsigned int _tip_id;


	BioFVMSinkSource* _source_O2;
	BioFVMSinkSource* _sink_VEGFA;
};

//BOOST_CLASS_EXPORT_KEY(Vas)

template<class Archive>
inline void Vas::serialize(Archive& ar, const unsigned int /* version */) {
	ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(Cell_Tumor);
	ar& BOOST_SERIALIZATION_NVP(_divide_limit);
	ar& BOOST_SERIALIZATION_NVP(_divide_cd);
	ar& BOOST_SERIALIZATION_NVP(_source_O2);
	ar& BOOST_SERIALIZATION_NVP(_sink_VEGFA);
}

};
};
