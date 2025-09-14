#pragma once

#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"
#include "TumorGridVoxel.h"
#include "VoxelContentGen.h"
#include "../agent/TCell.h"
#include "../agent/CancerCell.h"
#include "../agent/Mac.h"
#include "../agent/Fib.h"
#include "../agent/TCD4.h"
#include "../agent/MDSC.h"
#include "../agent/Vas.h"
#include "../../pde/DiffuseGrid.h"
#include "../../core/Stats.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
//#include <boost/serialization/export.hpp>

#include <algorithm> 
#include <set>

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

typedef std::vector<Coord3D> CellSource;
typedef boost::property_tree::ptree icProperty;

//! Generic tumor compartment
class Tumor : public SpatialCompartment
{
public:
enum TumExVar{
	TUMEX_CC,
	TUMEX_CC_DEATH,
	TUMEX_CC_T_KILL,
	TUMEX_CC_MAC_KILL,
	TUMEX_TEFF_REC,
	TUMEX_TH_REC,
	TUMEX_TREG_REC,
	TUMEX_mAPC_TRANSMIG,
	TUMEX_VAR_NUM
};

public:
	//! Need default constructor for serialization to work
	Tumor() {};
	Tumor(int x, int y, int z);
	virtual ~Tumor();

	// add one source of entry
	void add_immune_source();
	// reset all entry point
	void reset_immune_source();
	// add vasculature of entry
	void add_init_vas(const Coord3D& c);

	//! simulate for one slice
	virtual void timeSlice(unsigned long slice);

	//! recruit T Cells from invasive front 
	//void recruitTCellsFront(double min, double max);

	//! get stats
	const Stats& get_stats(void) const{ return _stats; };

	//! default header for extra remark column when writing cell grid to file
	virtual std::string getExtraRemarkHeader() const;
	//! print vaculature density grid to file
	std::string printStromaToFile();
	//! print grid snapshot to file
	std::string printGridToFile() const;
	//! print grid snapshot to file
	std::string printGradientToFile(int axis);
	//! print cell ODE stats to file
	void printCellOdeToFile(unsigned long slice) const;
	//! print grid snapshot to screen
	void printGridToScreen(unsigned long slice) const;
	//! get chemical grid
	DiffuseGrid & get_chem_grid(void) { return _chem; };
	//! get concentration of chemokine i.
	double get_chem(const Coord3D&c, chem_ID i)const; 
	// return ECM density given tumor coordinate
	double get_coord_ECM(const Coord3D& c) const;
	//! get 3D gradient vectors at voxel c.
	std::vector< std::vector<double> > get_gradient_voxel(const Coord3D& c);

	//! add to abm var exchange counter
	void inc_abm_var_exchange(TumExVar v){ _var_abm_to_qsp[v] += 1.0;};
	//! return variables needed for QSP module 
	const std::vector<double>& get_var_exchange(void);
	//! update ABM module with variables from QSP 
	void update_abm_with_qsp(const std::vector<double>&);

	


	// GETTER FUNCTIONS:
	// 
	//! get QSP cancer cell counts
	double get_QSP_cancer_counts(void)const { return _concentration_cc; };
	//! tumor volume
	double get_Tum_Vol(void)const{ return _tumor_volume; };
	//! Teff concentration in the tumor
	double get_Teff(void)const{ return _concentration_t_eff_tum; };
	//! TCD4 concentration in the tumor
	double get_Treg(void)const{ return _concentration_t_reg_tum; };
	//! TH concentration in the tumor
	double get_Th(void)const{ return _concentration_t_h_tum; };
	//! nivo concentration
	double get_Nivo(void)const{ return _concentration_nivo; };
	//! ipi concentration
	double get_Ipi(void)const { return _concentration_ipi; };
	//! ent concentration
	double get_Ent(void)const{ return _concentration_ent; };
	//! cabo concentration
	double get_Cabo(void)const { return _concentration_cabo; };
	// R_cabo
	double get_R_Cabo(void)const { return _R_cabo; };
	// cancer counts
	double get_cancer_counts(void)const { return _cancer_counts; };
	// cancer associated fibroblast counts
	double get_caf(void)const { return _concentration_caf; };
	// SETTER Function

	//! create one initial cell and configure to explicity initial condicitons
	void TEST_AddOneCell(AgentType type, AgentState state, const Coord3D& c) {
		createOneInitCell(type, state, c);
	};

	//! set allow shift
	void set_allow_shift(bool allowed){ _allow_shift_grid = allowed; };
	

	VoxelContentGen _voxel_ic;

protected:
	//! initialize agent grid 
	virtual bool initAgentGrid();
private:

	friend class boost::serialization::access;
	template<class Archive>
	//! boost serialization
	void serialize(Archive & ar, const unsigned int /*version*/);

	void time_slice_recruitment(void);
	void time_slice_movement(double t, double dt);
	void time_slice_state_change(double t, double dt);
	//! scan all agents for the last round in a slice
	void time_slice_final_scan(void);
	void time_slice_molecular(double t);

	//! initialize compartment: setup environmen��t
	void initEnvironment();
	//! initialize compartment: setup initial cells
	void initCell(std::string filename);
	//! initial cell: center cancer cell
	void init_cell_single_center(void);
	//! initial cell: fill grid 
	void init_cell_fill_grid(void);
	//! create one random cell
	Cell_Tumor* populate_voxel_random(const Coord3D&);
	//! check open voxel specifically for fibroblast
	bool getOpenNeighborForFib(Coord3D start, Coord3D& result);
	//! create a cluster of initial cells centered around provided coordinate
	void createClusterInitCell(icProperty &ic);
	//! create one initial cell and configure to explicity initial condicitons
	CellAgent* createOneInitCell(AgentType type, AgentState state, const Coord3D& crd);
	//! Teff/TCD4 recruitment probability at each entry point
	double get_T_recruitment_prob(double c, double base) const;
	// calculate recruitment probability for both macrophage and MDSCs (Myeloid cells) for each entry point
	double get_Myeloid_recruitment_prob(double c, double base) const;
	// update vasculature 
	void update_vas();
	// update ECM in the tumor
	void update_ECM();
	// update cabozantinib resistance
	void update_R_cabo();
	//! adjust camera center by shifting grid
	void shift_adjust_center(void);
	//! get the vector of camera shift
	Coord3D get_cam_shift(void);
	//! shift contents of the grid by crd
	void shift_grid(Coord3D& crd);
	// calculate distance of closest target (cells, agents, etc.)
	Grid3D<double> distances_to_nearest(std::function<bool(Coord3D const&)> callback);
	// make a new unique id
	unsigned int make_unique_id();
	// stats
	Stats _stats;
	// Diffusion grid
	DiffuseGrid _chem;
	//! list of init vas
	CellSource _init_vas;
	//! list of T cell sources
	CellSource _t_source;
	//! list of MDSC sources
	CellSource _mdsc_source;
	//! list of macrophage sources
	CellSource _mac_source;
	//! list of APC sources
	CellSource _APC_source;
	//! Dummy T cell for initialization
	TCell * _tInitDummy;
	//! Dummy MDSC for initialization
	MDSC * _MDSCInitDummy;	
	//! Dummy Cancer cell for initialization
	CancerCell * _cInitDummy;
	Mac * _macInitDummy;
	Fib * _fibInitDummy;
	TCD4 * _tcd4InitDummy;
	Vas* _VasInitDummy;

	//! voxel dimension, no need to serialize
	//double _voxel_size;

	/* variables used in grid shifting */

	//! allow shifting
	bool _allow_shift_grid;
	//! mass center anchor
	Coord3D _center_target;
	//! temporary AgentGrid
	AgentGrid _agGrid_temp;

	/* following are variables communicated between tumor and blood */

	//! variable values passed from QSP to ABM (no need to serialize)
	std::vector<double> _var_abm_to_qsp;
	//! concentration of cancer cells(no need to serialize).
	double _concentration_cc;
	//! concentration of cytotoxic t cells in blood(Unit: SI; no need to serialize).
	double _concentration_t_cyt;
	//! concentration of CD4 treg cells in blood(Unit: SI; no need to serialize).
	double _concentration_t_reg;
	//! concentration of CD4 th cells in blood(Unit: SI; no need to serialize).
	double _concentration_t_h;
	//! concentration of cytotoxic t cells in the tumor(Unit: SI; no need to serialize).
	double _concentration_t_eff_tum;	
	//! concentration of regulatory t cells in the tumor(Unit: SI; no need to serialize).
	double _concentration_t_reg_tum;
	//! concentration of th cells in the tumor(Unit: SI; no need to serialize).
	double _concentration_t_h_tum;
	//! concentration of macrophage m1 in the tumor(Unit: SI; no need to serialize).
	double _concentration_m1;
	//! concentration of macrophage m2 in the tumor(Unit: SI; no need to serialize).
	double _concentration_m2;
	//! concentration of mdsc in the tumor(Unit: SI; no need to serialize).
	double _concentration_mdsc;
	//! concentration of dead cancer cells(Unit: SI; no need to serialize).
	double _concentration_cx;
	//! concentration of exhausted t cells(Unit: SI; no need to serialize).
	double _concentration_t_exh;			
	//! concentration of nivo in Tumor(Unit: SI; no need to serialize).
	double _concentration_nivo;
	//! concentration of ipi in Tumor(Unit: SI; no need to serialize).
	double _concentration_ipi;
	//! concentration of ent in Tumor(Unit: SI; no need to serialize).
	double _concentration_ent;
	//! concentration of cabo in Tumor(Unit: SI; no need to serialize).
	double _concentration_cabo;
	//! cancer associated fibroblast counts
	double _concentration_caf;
	//! tumor volume(Unit: SI; no need to serialize).
	double _tumor_volume;	
	//! fraction of cancer cell to carrying capacity
	double _f_tum_cap;
	//! maximum number of cancer cell in the QSP model
	double _cmax;
	//! Resistance of cabozantinib
	double _R_cabo;
	//! cancer cell count
	double _cancer_counts;
	
	unsigned int _next_unique_id;

};

//BOOST_CLASS_EXPORT_KEY(Tumor);

template<class Archive>
inline void Tumor::serialize(Archive & ar, const unsigned int /* version */) {
	ar.template register_type<CancerCell>();
	ar.template register_type<TCell>();
	ar.template register_type<MDSC>();
	ar.template register_type<Mac>();
	ar.template register_type<Fib>();
	ar.template register_type<TCD4>();
	ar.template register_type<Vas>();
	ar.template register_type<TumorGridVoxel>();
	//ar.template register_type<BioFVMSinkSource>();
	// this is needed because cell._compartment is of type SpatialCompartment*
	ar.template register_type<Tumor>();
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpatialCompartment);
	ar & BOOST_SERIALIZATION_NVP(_voxel_ic);
	ar & BOOST_SERIALIZATION_NVP(_stats);
	ar & BOOST_SERIALIZATION_NVP(_t_source);
	ar & BOOST_SERIALIZATION_NVP(_init_vas);
	ar & BOOST_SERIALIZATION_NVP(_mdsc_source);
	ar & BOOST_SERIALIZATION_NVP(_mac_source);
	ar & BOOST_SERIALIZATION_NVP(_chem);
	ar & BOOST_SERIALIZATION_NVP(_allow_shift_grid);
	ar & BOOST_SERIALIZATION_NVP(_center_target);
	/*
	ar & BOOST_SERIALIZATION_NVP(_agGrid_temp);
	ar & BOOST_SERIALIZATION_NVP(_tInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_cInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_macInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_fibInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_tcd4InitDummy);
	*/
}

};
};

