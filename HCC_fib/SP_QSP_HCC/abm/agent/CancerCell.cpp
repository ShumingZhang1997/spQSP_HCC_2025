//#include <boost/serialization/export.hpp>
#include "CancerCell.h"

//BOOST_CLASS_EXPORT_IMPLEMENT(CancerCell)
#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"
#include "../compartment/Tumor.h"
#include "../compartment/LymphCentral.h"
#include "../../core/GlobalUtilities.h"


namespace SP_QSP_IO{
namespace SP_QSP_HCC{

using std::string;
using std::stringstream;

static int CancerCellSize = 1;

CancerCell::CancerCell()
	:_count_neighbor_Teff(0)
{
}

CancerCell::CancerCell(SpatialCompartment* c)
	:Cell_Tumor(c)
	, _stemID(getID())
	, _divideCD(int(params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE)) + .5)
	, _divideFlag(true)
	, _divideCountRemaining(0)
	, _count_neighbor_Teff(0)
	, _count_neighbor_MacM1(0)
	, _source_CCL2(NULL)
	, _source_TGFB(NULL)
	, _source_VEGFA(NULL)
	, _sink_O2(NULL)
	, _sink_IFNg(NULL)
{
	_state = AgentStateEnum::CANCER_STEM;
}


CancerCell::CancerCell(const CancerCell& c)
	:Cell_Tumor(c)
	, _stemID(getID())
	, _divideCD(c._divideCD)
	, _divideFlag(c._divideFlag)
	, _divideCountRemaining(c._divideCountRemaining)
	, _count_neighbor_Teff(0)
	, _count_neighbor_MacM1(0)
	, _source_CCL2(NULL)
	, _source_TGFB(NULL)
	, _source_VEGFA(NULL)
	, _sink_O2(NULL)
	, _sink_IFNg(NULL)
{
	if (_state == AgentStateEnum::CANCER_STEM)
	{
		setup_chem_source(_source_CCL2, CHEM_CCL2, params.getVal(PARAM_CCL2_RELEASE));
		//std::cout << "Stem CCL2 source pointer: " << _source_CCL2 << std::endl;
		setup_chem_source(_source_TGFB, CHEM_TGFB, params.getVal(PARAM_CANCER_TGFB_RELEASE));
		setup_chem_source(_source_VEGFA, CHEM_VEGFA, params.getVal(PARAM_CANCER_STEM_VEGFA_RELEASE));
		//sink
		setup_chem_sink(_sink_O2, CHEM_O2, params.getVal(PARAM_O2_UPTAKE));
		setup_chem_sink(_sink_IFNg, CHEM_IFN, params.getVal(PARAM_CANCER_IFN_G_UPTAKE));
	}
	else if (_state == AgentStateEnum::CANCER_PROGENITOR)
	{
		//source
		setup_chem_source(_source_CCL2, CHEM_CCL2, params.getVal(PARAM_CCL2_RELEASE));
		//std::cout << "Progenitor CCL2 source pointer: " << _source_CCL2 << std::endl;
		setup_chem_source(_source_TGFB, CHEM_TGFB, params.getVal(PARAM_CANCER_PROGENITOR_TGFB_RELEASE));
		setup_chem_source(_source_VEGFA, CHEM_VEGFA, params.getVal(PARAM_CANCER_PRO_VEGFA_RELEASE));
		//sink
		setup_chem_sink(_sink_O2, CHEM_O2, params.getVal(PARAM_O2_UPTAKE));
		setup_chem_sink(_sink_IFNg, CHEM_IFN, params.getVal(PARAM_CANCER_IFN_G_UPTAKE));
	}
}

CancerCell::~CancerCell()
{
}

string CancerCell::toString()const{
	stringstream ss;
	ss << CellAgent::toString();
	ss << "division flag: " << _divideFlag << ", division cool down: " << _divideCD << std::endl;
	return ss.str();
}

bool CancerCell::agent_movement_step(double t, double dt, Coord& c){


	bool move = false;
	// if ECM density is high, cell is unlikely to move
	double ECM_density = get_tumor().get_coord_ECM(c);
	//double ECM_sat = (ECM_density - params.getVal(PARAM_FIB_ECM_BASELINE)) / (params.getVal(PARAM_FIB_ECM_SATURATION) - params.getVal(PARAM_FIB_ECM_BASELINE));
	//std::cout << "ECM_density: " << ECM_density << ", ECM_sat: " << ECM_sat << std::endl;
	double ECM_sat = ECM_density / (ECM_density + params.getVal(PARAM_FIB_ECM_MOT_EC50));
	if (rng.get_unif_01() < ECM_sat){
		return move;
	}

	// move
	int idx;
	const auto shape = getCellShape();
	if (_compartment->getOneOpenVoxel(shape->getMoveDestinationVoxels(),
		shape->getMoveDirectionAnchor(), getCoord(), getType(), idx, rng))
	{
		move = true;
		c = getCellShape()->getMoveDirectionAnchor()[idx] + getCoord();
	}

	return move;

}

bool CancerCell::agent_state_step(double t, double dt, Coord& c){

	bool divide = false;
	
	
	if (!isDead() && _state == AgentStateEnum::CANCER_SENESCENT)
	{
		_life -= 1;
		if (_life <= 0)
		{
			setDead();
		}
	}
	
	// senescence 
	//std::cout << _divideCountRemaining << std::endl;
	if (_state == AgentStateEnum::CANCER_PROGENITOR && _divideCountRemaining <= 0) {
		setSenescent();
	}

	//hypoxia
	if (_state == AgentStateEnum::CANCER_PROGENITOR) {
		double O2 = get_tumor().get_chem(getCoord(), CHEM_O2);
		if (O2 < params.getVal(PARAM_CANCER_HYPOXIA_TH)) {
			_divideFlag = false;
		} else {
			_divideFlag = true;
		}
	}
	const auto shape = getCellShape();

	// check IFNg
	Cell_Tumor::agent_state_step(t, dt, c);

	auto tumor = dynamic_cast<Tumor*>(_compartment);
	double nivo = tumor->get_Nivo();
	double ent = tumor->get_Ent();
	// decide if get killed by Tcyt
	//if (_state != AgentStateEnum::CANCER_STEM)
	{
		double cancer_count = double(get_tumor().get_cancer_counts());
		double cancer_count_max = double(_compartment->get_sizeX()) * double(_compartment->get_sizeY()) * double(_compartment->get_sizeZ());
		double cancer_count_ratio = cancer_count / cancer_count_max;
		//std::cout << "cancer cell count: " << cancer_count << ", cancer count ratio: " << cancer_count_ratio << std::endl;
		//std::cout << "cancer cell count: " << cancer_count << std::endl;
		if (_count_neighbor_Teff > 0 && cancer_count_ratio > 0.1) {
			// count neighborhood cancer cells
			int neighbor_cc = 0;
			_compartment->for_each_neighbor_ag(shape->getEnvironmentLocations(),
				getCoord(), [&](BaseAgent* ag){
				if (ag->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
					neighbor_cc += 1;
				}
				return true;
			});
			//std::cout << "PDL1 expression: " << _PDL1_syn << std::endl;
			double bond = TCell::get_PD1_PDL1(_PDL1_syn, nivo);
			double supp = TCell::get_PD1_supp(bond, params.getVal(PARAM_N_PD1_PDL1));
			double NO = get_tumor().get_chem(getCoord(), CHEM_NO);
			double ArgI = get_tumor().get_chem(getCoord(), CHEM_ARGI);
			double TGFB = get_tumor().get_chem(getCoord(), CHEM_TGFB);
			// double H_mdsc_c1 = 1 - ((1 - (NO / (NO + params.getVal(PARAM_MDSC_IC50_NO_CTL)))) * (1 - (ArgI / (ArgI + params.getVal(PARAM_MDSC_IC50_ArgI_CTL))))); //(NO entinostats) *(1-(ent/(ent+params.getVal(PARAM_IC50_ENT_ARGI)))));
			double H_mdsc_c1 = 1 - (1 - (ArgI / (ArgI + params.getVal(PARAM_MDSC_IC50_ArgI_CTL)))); //(NO entinostats) *(1-(ent/(ent+params.getVal(PARAM_IC50_ENT_ARGI)))));
			//(NO entinostats) *(1-(ent/(ent+params.getVal(PARAM_IC50_ENT_ARGI)))));
			double H_TGFB = (TGFB / (TGFB + params.getVal(PARAM_TEFF_TGFB_EC50)));
			double q = (1 - H_mdsc_c1) * double(_count_neighbor_Teff) / (_count_neighbor_Teff + neighbor_cc + params.getVal(PARAM_CELL)) * (1 - H_TGFB);

			//std::cout << "PARAM_ESCAPE_BASE : " << params.getVal(PARAM_ESCAPE_BASE) << ", q: " << q << ", supp: " << supp << std::endl;
			double ECM_density = get_tumor().get_coord_ECM(c);
			double ECM_sat = ECM_density / (ECM_density + params.getVal(PARAM_FIB_ECM_MOT_EC50));

			double cabo = get_tumor().get_Cabo();
			double R_cabo = get_tumor().get_R_Cabo();
			double killing_scaler = params.getVal(PARAM_ABM_KILL_SCALE) * (1 - supp);

			double p_kill = TCell::get_kill_prob(supp, q) * killing_scaler;
			

			/*
			if (rng.get_unif_01() < 0.05) {
				std::cout << "T cell killing:\n"
					<< "neighbors: " << _count_neighbor_Teff << ", neighbor_cc:" << neighbor_cc<< ", q: " << q << ", H_mdsc_c1: " << H_mdsc_c1 << std::endl;
				
				
			
				
				std::cout << "H_mdsc_c1 Variable: " << H_mdsc_c1
					<< ", NO: " << NO << ", ArgI: " << ArgI << "\n PARAM_IC50_NO_CTL: "
					<< params.getVal(PARAM_MDSC_IC50_NO_CTL) << ", PARAM_IC50_ARGI_CTL: " << params.getVal(PARAM_MDSC_IC50_ArgI_CTL)  << ", PARAM_IC50_ENT_ARGI: " << params.getVal(PARAM_IC50_ENT_ARGI) << std::endl;
				
				std::cout << "H_NO: " << (NO / (NO + params.getVal(PARAM_MDSC_IC50_NO_CTL)))
					<< ", H_ARGI: " << (ArgI / (ArgI + params.getVal(PARAM_MDSC_IC50_ArgI_CTL))) <<  std::endl;

				

				
				std::cout 
					<< "T cell PDL1: " << _PDL1_syn << ", pd1_PDL1_conc: " << bond 
					<< ", PD1_PDL1 Half max: " << params.getVal(PARAM_PD1_PDL1_HALF) << ", supp: " << supp << ", killing scaler: " << killing_scaler
					<< ", q: "<< q << ", TGFB: " << TGFB << ", PARAM_TEFF_TGFB_EC50: " << params.getVal(PARAM_TEFF_TGFB_EC50)<< ", H_TGFB: " << H_TGFB 
					<< ", pkill: " << p_kill 
					<< std::endl;
				
			}
			*/

			if (rng.get_unif_01() < p_kill){
				setDead();
				tumor->inc_abm_var_exchange(Tumor::TUMEX_CC_T_KILL);
				//total_kill += 1;
				/*
				if (_state = AgentStateEnum::CANCER_PROGENITOR)
				{
					std::cout << "total kill: " << total_kill << std::endl;
				}*/
				//std::cout << "Killed: " << _state 
				// << ": " << p_kill << std::endl;

			}
		}
		if (_count_neighbor_MacM1 > 0 && cancer_count_ratio > 0.1) {
			// count neighborhood cancer cells
			int neighbor_cc = 0;
			_compartment->for_each_neighbor_ag(shape->getEnvironmentLocations(),
				getCoord(), [&](BaseAgent* ag) {
					if (ag->getType() == AgentTypeEnum::CELL_TYPE_CANCER) {
						neighbor_cc += 1;
					}
					return true;
				});
			//std::cout << "PDL1 expression: " << _PDL1_syn << std::endl;
			double PD1_bond = Mac::get_PD1_PDL1(_PDL1_syn, nivo);
			double IL10 = get_tumor().get_chem(getCoord(), CHEM_IL_10);
			double TGFB = get_tumor().get_chem(getCoord(), CHEM_TGFB);

			//V_T.C1 -> V_T.C_x =  k_M1_phago*V_T.C1*V_T.Mac_M1/(V_T.Mac_M1+K_Mac_C*V_T.C1+cell)*(1-H_Mac_C)*(1-H_IL10_phago)
			//Calculate H_PD1 of macrophage as coefficient for cancer cell killing
			// params.getVal(PARAM_PD1_PDL1_HALF) is in bonds/ um^2. 
			// PD1_bond is in bonds Need to divided by params.getVal(PARAM_A_SYN)
			double H_PD1_M = Mac::get_H(PD1_bond / params.getVal(PARAM_A_SYN), params.getVal(PARAM_N_PD1_PDL1), params.getVal(PARAM_PD1_PDL1_HALF));

			double kd = (params.getVal(PARAM_KON_SIRPa_CD47) / params.getVal(PARAM_KOFF_SIRPa_CD47));
			double a = params.getVal(PARAM_C1_CD47_SYN);
			double b = params.getVal(PARAM_MAC_SIRPa_SYN);
			//Calculate H_SIRPa_CD47 of macrophage as coefficient for cancer cell killing
			double SIRPa_CD47_conc = ((a + b + 1 / kd) - std::sqrt((a + b + 1 / kd) * (a + b + 1 / kd) - 4 * a * b)) / 2;
			double SIRPa_CD47_k50 = params.getVal(PARAM_MAC_SIRPa_HALF);

			// H_SIRPa_CD47_M is in bonds / um^2, No Need to divided by params.getVal(PARAM_A_SYN)
			double H_SIRPa_CD47_M = Mac::get_H(SIRPa_CD47_conc, params.getVal(PARAM_N_SIRPa_CD47), SIRPa_CD47_k50);
			double H_Mac_C = 1 - (1 - H_SIRPa_CD47_M) * (1 - H_PD1_M);

			//Calculate H_IL10_phago of macrophage as coefficient for cancer cell killing
			double H_IL10_phago = (IL10 / (IL10 + params.getVal(PARAM_MAC_IL_10_HALF_PHAGO)));

			

			double q = double(_count_neighbor_MacM1) / (_count_neighbor_MacM1 + neighbor_cc + params.getVal(PARAM_CELL)) * (1 - H_Mac_C) * (1 - H_IL10_phago);
			//Get probability of cancer cell killed by macrophage
			double p_kill = Mac::get_kill_prob(params.getVal(PARAM_ESCAPE_MAC_BASE), q);
			
		


			/*
			std::cout
				<< "Macrophage, PDL1: " << _PDL1_syn << ", PD1_bond: " << PD1_bond 
				<< ", Cancer CD47_SYN: " << params.getVal(PARAM_C1_CD47_SYN)
				<< ", Macrophage SIRPa_SYN: " << params.getVal(PARAM_MAC_SIRPa_SYN)
				<< ", kon SIRPa_CD47: " << params.getVal(PARAM_KON_SIRPa_CD47)
				<< ", koff SIRPa_CD47: " << params.getVal(PARAM_KOFF_SIRPa_CD47)
				<< ", MAC_SIRPa_HALF: " << params.getVal(PARAM_MAC_SIRPa_HALF)
				<< ", SIRPa_CD47_bond concentration: " << SIRPa_CD47_conc
				<< ", H_SIRPa_CD47_M: " << H_SIRPa_CD47_M
				<< ", H_PD1_M: " << H_PD1_M << "\n"
				<< "pkill: " << p_kill << std::endl;
			*/
			if (rng.get_unif_01() < p_kill) {
				setDead();
				tumor->inc_abm_var_exchange(Tumor::TUMEX_CC_MAC_KILL);
				//total_kill += 1;
				/*
				if (_state = AgentStateEnum::CANCER_PROGENITOR)
				{
					std::cout << "total kill: " << total_kill << std::endl;
				}
				*/
				//std::cout << "Killed: " << _state 
				// << ": " << p_kill << std::endl;
			}
		}
	}

	if (isDead())
	{
		return divide;
	}

	/*
	if (_state == AgentStateEnum::CANCER_PROGENITOR) {
		double O2 = get_tumor().get_chem(getCoord(), CHEM_O2);
		//std::cout << getCoord() << ": " << O2 << std::endl;
		if (O2 < params.getVal(PARAM_CANCER_HYPOXIA_TH)) {
			update_chem_source(_source_VEGFA, params.getVal(PARAM_CANCER_PRO_VEGFA_RELEASE));
		}
		
	}
	*/
	//printCellInfo2(t, this ,"before divide/move");
	// divide
	if (_divideCD > 0)
	{
		_divideCD--;
	}
	
	//std::cout << "Cancer Divide CD: " << _divideCD << std::endl;

	if (_divideFlag && _divideCD <= 0)
	{
		// find location to divide
		// if ECM density is high, cell is unlikely to move
		//double ECM_density = get_tumor().get_coord_ECM(c);
		//double ECM_sat = (ECM_density - params.getVal(PARAM_FIB_ECM_BASELINE)) / (params.getVal(PARAM_FIB_ECM_SATURATION) - params.getVal(PARAM_FIB_ECM_BASELINE));
		//double ECM_coeff =  1 - ECM_sat;
		//if (rng.get_unif_01() < ECM_coeff){
		//	return divide;
		//}

		int idx;
		if (_compartment->getOneOpenVoxel(shape->getProlifDestinationVoxels(),
			shape->getProlifDestinationAnchor(), getCoord(), getType(), idx, rng))
		{
			divide = true;
			//divide = false;
			//cout << "idx: " << idx << ", " << getCellShape()->getProlif()[idx] << endl;
			c = getCellShape()->getProlifDestinationAnchor()[idx] + getCoord();

			//_divideFlag = true;
			if (_state == AgentStateEnum::CANCER_STEM)
			{
				double cabo = get_tumor().get_Cabo();
				double R_cabo = get_tumor().get_R_Cabo();
				double cabo_prolif_factor = 1 - (params.getVal(PARAM_LAMBDA_C_CABO) * cabo / (cabo + params.getVal(PARAM_IC50_MET))) * R_cabo;
				//_divideCD = int(params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE) / (1-H_ent_c1) + .5);
				_divideCD = int((params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE) / cabo_prolif_factor) + 0.5);

				//std::cout << "cabo: " << cabo << ", PARAM_IC50_MET: " << params.getVal(PARAM_IC50_MET) << ", PARAM_LAMBDA_C_CABO: " << params.getVal(PARAM_LAMBDA_C_CABO) 
				//<< ", R_cabo: " << R_cabo << ", cabo_prolif_factor: " << cabo_prolif_factor << std::endl;
				//_divideCD = int(params.getVal(PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) / (1-H_ent_c1) + .5);
				//std::cout << "Cancer Divide CD: " << _divideCD << std::endl;
				//_divideCD = int(params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE) + 0.5) ;
			}
			else if (_state == AgentStateEnum::CANCER_PROGENITOR) {
				_divideCountRemaining -= 1;
				_divideCD = int(params.getVal(PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) + 0.5);
				//double O2 = get_tumor().get_chem(getCoord(), CHEM_O2);
				//if (O2 < params.getVal(PARAM_CANCER_HYPOXIA_TH)) {
				//	_divideCD *= 2;
				//}
			}
		}
	}

	
	//printCellInfo2(t, this ,"after divide/move");

	// do other stuff

	_count_neighbor_Teff = 0;
	_count_neighbor_MacM1 = 0;

	return divide;
}

/*! change cell state of a CSC to progenitor.
	-# daughter cell from asymmetric division 
	-# during initialization (default state is stem)
*/
void CancerCell::setProgenitor() {
	_state = AgentStateEnum::CANCER_PROGENITOR;
	_divideCountRemaining = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
	
	_divideCD = int(params.getVal(PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) + 0.5);
	// only stem cell secrete TGFB 
	update_chem_source(_source_TGFB, params.getVal(PARAM_CANCER_PROGENITOR_TGFB_RELEASE));
	update_chem_source(_source_VEGFA, params.getVal(PARAM_CANCER_PRO_VEGFA_RELEASE));
}

/*! set cancer cell to senescent
*/
void CancerCell::setSenescent(){
	_state = AgentStateEnum::CANCER_SENESCENT;
	_divideCD = -1;
	_divideFlag = false;
	_life = getSenescentLife();

	// Senescent cells no longer secreting Cytokines
	update_chem_source(_source_CCL2, 0);
	update_chem_source(_source_TGFB, 0);
	update_chem_source(_source_VEGFA, 0);
	return;
}
/*! randomize division cooldown 
	Set to a number between [1, mean], uniform
*/
void CancerCell::randomize_div_cd(int mean){
	_divideCD = int(rng.get_unif_01()*mean) + 1;
	return;
}

/*! remaining slices for senescent cancer cells
	randomly drawn from exponential distribution.
*/
int CancerCell::getSenescentLife(void){
	double mean = params.getVal(PARAM_CANCER_SENESCENT_MEAN_LIFE);
	return int(rng.get_exponential(mean) + .5);
	//cout << "random cancer cell life: " << cLife << endl;
}

//! remarks to string
std::string CancerCell::getRemark() const{
	stringstream ss;
	ss << Cell_Tumor::getRemark()
		<< "|" << _stemID
		<< "|" << _divideCD
		<< "|" << _divideCountRemaining;
	return ss.str();
}

//! move sources (CCL2)
void CancerCell::move_all_source_sink(void)const{
	//std::cout << "moving sources: " << getCoord() << std::endl;
	move_source_sink(_source_CCL2);
	move_source_sink(_source_TGFB);
	move_source_sink(_source_VEGFA);
	move_source_sink(_sink_O2);
	move_source_sink(_sink_IFNg);
return;
}

//! remove sources (CCL2)
void CancerCell::remove_all_source_sink(void){

	remove_source_sink(_source_CCL2);
	remove_source_sink(_source_TGFB);
	remove_source_sink(_source_VEGFA);
	remove_source_sink(_sink_O2);
	remove_source_sink(_sink_IFNg);
	
	return;
}
};
};
