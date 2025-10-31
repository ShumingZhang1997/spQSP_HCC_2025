#pragma once

#include "Cell_Tumor.h"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>


namespace SP_QSP_IO{
namespace SP_QSP_HCC{

class Fib:
	public Cell_Tumor
{
public:
	Fib(){};
	Fib(SpatialCompartment* c );
	Fib(const Fib& c);
	virtual ~Fib();

	virtual CellAgent* createCellCopy() const { return new Fib(*this); };

	//! move sources (ECM)
	void move_all_source_sink(void)const;

	//! remove all sources
	void remove_all_source_sink(void);

	//! print cancer cell information
	virtual std::string toString() const;

	virtual bool agent_movement_step(double t, double dt, Coord& c);
	virtual bool agent_state_step(double t, double dt, Coord& c);
	//virtual bool agent_state_step(double t, double dt, Coord& c);
	//! step function for cancer cell
	//virtual void agentStep(double t, double dt, AgentStep & as);
	//! get cell agent type, in this case cancer cell.
	virtual AgentType getType() const { return AgentTypeEnum::CELL_TYPE_FIB; };

	//get and set the pointer of the linked fibroblast
	Fib* getNext() { return _next; };
	Fib* getPrevious() { return _previous; };
	void setNext(Fib* next);
	void setPrevious(Fib* previous);
	//set fibroblast state to CAFs
	void setCAF();


private:
	friend class boost::serialization::access;


	Fib* _previous;
	Fib* _next;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);
};

//BOOST_CLASS_EXPORT_KEY(Fib)

template<class Archive>
inline void Fib::serialize(Archive & ar, const unsigned int /* version */){
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Cell_Tumor);
	ar& BOOST_SERIALIZATION_NVP(_previous);
	ar& BOOST_SERIALIZATION_NVP(_next);

}

};
};
