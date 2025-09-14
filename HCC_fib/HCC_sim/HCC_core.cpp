#include "HCC_core.h"

//#include "HCC/SP_QSP_HCC/core/Param.h"
#include "HCC/SP_QSP_HCC/core/GlobalUtilities.h"
#include "InitialCondition.h"

#include <algorithm>    // std::max

extern FileOutputHub output_hub;

extern RNG rng;

//extern SP_QSP_IO::Param params;

namespace SP_QSP_IO {
	namespace SP_QSP_HCC {
		extern Param params;
	}
};
static auto& params = SP_QSP_IO::SP_QSP_HCC::params;

extern InitialCondition ic;
extern std::string initialCellFileName_core;
extern std::string initialCellFileName_margin;

typedef SP_QSP_IO::Coord Coord;

HCC_Core::HCC_Core()
: _tumor_core(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X),
	params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y),
	params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z))
, _tumor_margin(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X),
	params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y),
	params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z))
, _lymph()
{


}

HCC_Core::HCC_Core(std::string core_antigen_file, std::string margin_antigen_file)
	: _tumor_core(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X),
		params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y),
		params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z))
	, _tumor_margin(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X),
		params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y),
		params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z))
	, _lymph()
{
	_tumor_core.loadAntigenData(core_antigen_file);
	_tumor_margin.loadAntigenData(margin_antigen_file);

}

HCC_Core::~HCC_Core()
{
}

/*! Setup QSP module.
*/
void HCC_Core::setup_qsp(CancerVCT::Param& p){
	_lymph.setup_param(p);
	params.update_from_qsp();
}
/*! initialize compartments: randomly populate voxels
	This function is called only when creating the model for the first time,
	and not reading from saved state. Objects created in this function should
	be already in the serialization and can be loaded directly.
*/
void HCC_Core::initializeSimulation(void){

	// rule to randomly populate voxel during initlaization or grid shifting
	_tumor_margin.set_allow_shift(ic.getVal(IC_MARGIN_GRID_SHIFT));
	_tumor_margin._voxel_ic.setup(ic.getVal(IC_MARGIN_STATIONARY),
		ic.getVal(IC_DENSITY_MARGIN_CANCER),
		params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X),
		params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y),
		ic.getVal(IC_MARGIN_CANCER_BOUNDARY));

	_tumor_core.set_allow_shift(ic.getVal(IC_CORE_GRID_SHIFT));
	_tumor_core._voxel_ic.setup(ic.getVal(IC_CORE_STATIONARY),
			ic.getVal(IC_DENSITY_CORE_CANCER),
			0, 0, 0);


	//t cell sources
	// margin:
	{
		std::vector<Coord> c_normal, c_tumor;
		unsigned int nr_source_normal, nr_source_tumor;
		_tumor_margin.for_each_grid_coord(true, true, true, [&](Coord&c){
			if (c.z > ic.getVal(IC_MARGIN_CANCER_BOUNDARY)){
				c_normal.push_back(c);
			}
			else{
				c_tumor.push_back(c);
			}
		});
		nr_source_normal = int(c_normal.size()*ic.getVal(IC_MARGIN_NORMAL_VAS_FOLD)
			* params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_REC_PORT_PROB));
		rng.shuffle_first_k(c_normal, nr_source_normal);
		nr_source_tumor = int(c_tumor.size()*ic.getVal(IC_MARGIN_TUMOR_VAS_FOLD)
			* params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_REC_PORT_PROB));
		rng.shuffle_first_k(c_tumor, nr_source_tumor);

		for (size_t i = 0; i < nr_source_normal; i++)
		{
			//_tumor_margin.add_lymphocyte_source(c_normal[i]);
			//_tumor_margin.add_mdsc_source(c_normal[i]);
		}
		for (size_t i = 0; i < nr_source_tumor; i++)
		{
			_tumor_margin.add_lymphocyte_source(c_tumor[i]);
			_tumor_margin.add_mdsc_source(c_tumor[i]);
		}
		std::cout << "margin nr sources: tumor: " << nr_source_tumor 
			<< ", normal: " << nr_source_normal << std::endl;
	}
	// core:
	{
		std::vector<Coord> c_tumor;
		unsigned int nr_source_tumor;
		_tumor_margin.for_each_grid_coord(true, true, true, [&](Coord&c){
			c_tumor.push_back(c);
		});
		nr_source_tumor = int(c_tumor.size()*ic.getVal(IC_CORE_TUMOR_VAS_FOLD)
			* params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_REC_PORT_PROB));
		rng.shuffle_first_k(c_tumor, nr_source_tumor);
		for (size_t i = 0; i < nr_source_tumor; i++)
		{
			_tumor_core.add_lymphocyte_source(c_tumor[i]);
			_tumor_core.add_mdsc_source(c_tumor[i]);
		}
		std::cout << "core nr sources: tumor: " << nr_source_tumor << std::endl;
	}

	std::string s;
	_tumor_core.initCompartment(s);
	_tumor_margin.initCompartment(s);
}

/*! initialize compartments: create initial cells from file input specifications
	This function is called only when creating the model for the first time,
	and not reading from saved state. Objects created in this function should
	be already in the serialization and can be loaded directly.
*/
void HCC_Core::initializeSimulation(std::string core, std::string margin){

	_tumor_core.set_allow_shift(false);
	_tumor_margin.set_allow_shift(false);

	_tumor_core.initCompartment(core);
	_tumor_margin.initCompartment(margin);
}

void HCC_Core::timeSlice(const long slice){
	
	const double dt = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_SEC_PER_TIME_SLICE);
	const double t0 = slice * dt;
	std::cout << slice << std::endl;
	// std::cout << "RNG check (" << slice << ") START : " << rng.get_unif_01() << std::endl;

	/* update cancer number and blood concentration */
	auto& qsp_var = _lymph.get_var_exchange();
	double lymphCC = qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_C];

	/*Update the ABM for All teff*/
	auto& qsp_var_teff = _lymph.get_var_exchange_teff();

	/* if QSP halted, skip*/
	std::cout << "lymph CC: " << lymphCC << std::endl;
	double abm_min_cc = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_C1_MIN);

	if (lymphCC > abm_min_cc)
	{
		_tumor_core.update_abm_with_qsp(qsp_var, qsp_var_teff);
		_tumor_margin.update_abm_with_qsp(qsp_var, qsp_var_teff);

		std::cout << "nivo: " << qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_NIVO] << std::endl;

		/*
		for (auto& v : qsp_var)
		{
		std::cout << v << ", ";
		}
		std::cout << std::endl;
		*/

		/* ABM time step */
		_tumor_core.timeSlice(slice);
		//std::cout << "RNG check (" << slice << ") CORE: " << rng.get_unif_01() << std::endl;
		_tumor_margin.timeSlice(slice);
		//std::cout << "RNG check (" << slice << ") MARGI : " << rng.get_unif_01() << std::endl;

		/* update QSP variables */
		auto& abm_var_0 = _tumor_core.get_var_exchange();
		auto& abm_var_1 = _tumor_margin.get_var_exchange();

		size_t abm_var_len = abm_var_0.size();
		auto abm_var = std::vector<double>(abm_var_len, 0);

		double w = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_WEIGHT_QSP);

		double fraction_margin = slice > params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_RESECT_TIME_STEP)
			? ic.getVal(IC_FRACTION_MARGIN)
			: ic.getVal(IC_FRACTION_MARGIN_RES);

		double fraction_core = 1 - fraction_margin;

		double tumCC_core = abm_var_0[SP_QSP_IO::SP_QSP_HCC::Tumor::TUMEX_CC];
		double tumCC_margin = abm_var_1[SP_QSP_IO::SP_QSP_HCC::Tumor::TUMEX_CC];

		double abm_scaler_core = (1 - w) / w * lymphCC / (tumCC_core + abm_min_cc )* fraction_core;
		double abm_scaler_margin = (1 - w) / w * lymphCC / (tumCC_margin  + abm_min_cc) * fraction_margin;

		std::cout << "scalor:\n" << "core:" << abm_scaler_core 
			<< "\nmargin: " << abm_scaler_margin << std::endl;

		for (size_t i = 0; i < abm_var_len; i++)
		{
			abm_var[i] = abm_var_0[i] * abm_scaler_core
				+ abm_var_1[i] * abm_scaler_margin;
		}

		_lymph.update_qsp_var(abm_var);

	}

	/* QSP time step */
	_lymph.time_step(t0, dt);
	//std::cout << "RNG check (" << slice << ") QSP: " << rng.get_unif_01() << std::endl;
	return;

}

void HCC_Core::write_stats_header(void) const {

	auto& statsStream_core = output_hub.getStatsFstream(0);
	statsStream_core  << _tumor_core.get_stats().writeHeader();
	auto& statsStream_margin = output_hub.getStatsFstream(1);
	statsStream_margin  << _tumor_margin.get_stats().writeHeader();
	return;
}

void HCC_Core::write_stats_slice(unsigned long slice)const{
	{
		auto& statsStream = output_hub.getStatsFstream(0);
		statsStream << _tumor_core.get_stats().writeSlice(slice);
		statsStream.flush();
	}
	{
		auto& statsStream = output_hub.getStatsFstream(1);
		statsStream << _tumor_margin.get_stats().writeSlice(slice);
		statsStream.flush();
	}
	return;
}


void HCC_Core::write_QSP(unsigned long slice, bool header)const{
	auto& stream = output_hub.get_lymph_blood_QSP_stream();
	if (header){
		stream << "time" << _lymph.getQSPHeaders() << std::endl;
	}
	else{
		stream << slice << _lymph << std::endl;
	}
	return;
}

void HCC_Core::writeOde(unsigned long slice){
	_tumor_core.printCellOdeToFile(slice);
}


void HCC_Core::briefStats(unsigned long slice){
	std::cout << "Time: " << slice << std::endl;
	{
		const auto& stats = _tumor_core.get_stats();
		std::cout << "Core: " << "nrCell: " << _tumor_core.getNrCell()
			<< ", CD8: " << stats.getTCell()
			<< ", Treg: " << stats.getTreg()
			<< ", MDSC: " << stats.getMDSC()
			<< ", Cancer cell:" << stats.getCancerCell() << std::endl;
	}
	{
		const auto& stats = _tumor_margin.get_stats();
		std::cout << "Margin: " << "nrCell: " << _tumor_margin.getNrCell()
			<< ", CD8: " << stats.getTCell()
			<< ", Treg: " << stats.getTreg()
			<< ", MDSC: " << stats.getMDSC()
			<< ", Cancer cell:" << stats.getCancerCell() << std::endl;
	}
}

/*! Print grid info to file.
    \param [in] slice
	\param [in] option: 1. only cellular scale; 2. only molecular scale; 3. both scales
*/
void HCC_Core::writeGrids(unsigned long slice, unsigned int option){
	if (option == 1 || option == 3)
	{
		{
			std::ofstream& snap = output_hub.getNewGridToSnapshotStream(slice, "cell_core_");
			snap << _tumor_core.compartment_cells_to_string();
			snap.close();
		}
		{
			std::ofstream& snap = output_hub.getNewGridToSnapshotStream(slice, "cell_margin_");
			snap << _tumor_margin.compartment_cells_to_string();
			snap.close();
		}
	}
	if (option == 2 || option == 3)
	{
		{
			std::ofstream&  snap = output_hub.getNewGridToSnapshotStream(slice, "grid_core_");
			snap << _tumor_core.printGridToFile();
			snap.close();

		}
		{
			std::ofstream&  snap = output_hub.getNewGridToSnapshotStream(slice, "grid_margin_");
			snap << _tumor_margin.printGridToFile();
			snap.close();
		}
	}
}

