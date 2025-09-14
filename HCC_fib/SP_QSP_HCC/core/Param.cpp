#include "Param.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <math.h>
#include "../ode/Param.h"
#include "../ode/ODE_system.h"

// parameters from QSP module. 
extern CancerVCT::Param qsp_params;

#define QP(x) CancerVCT::ODE_system::get_class_param(x)
#define AVOGADROS 6.022140857E23 
#define PI 3.1415926525897932384626

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

namespace pt = boost::property_tree;

static int SEC_PER_DAY = 86400;
static int HOUR_PER_DAY = 24;

// must match the order in the ParamInt and ParamFloat enums
#define PARAM_DESCRIPTION_FIELD_COUNT 3
const char* _description[][PARAM_DESCRIPTION_FIELD_COUNT] =
{
	//{"fullpath", "desc", "constraint"}

	//------------------------ float --------------------------//
	/* QSP */
	{ "Param.QSP.simulation.weight_qsp", "", "prob" },
	{ "Param.QSP.simulation.t_steadystate", "days", "pos" },
	{ "Param.QSP.simulation.t_resection", "days", "pos" },
	{ "Param.QSP.simulation.presimulation_diam_frac", "", "pos" },
	/* ABM */
	//environmental
	{ "Param.ABM.Environment.SecPerSlice", "", "pos" },
	{ "Param.ABM.Environment.recSiteFactor", "number of adhesion site per port voxel", "pos" },
	//Adhesion molecule density Reference
	//Increased ICAM-1 Expression Causes Endothelial Cell Leakiness, Cytoskeletal Reorganizationand Junctional Alterations
	// 657 molecules per field --> 1 field = 2.06um * 2.06um --> 657/(2.06 * 2.06) um
	{ "Param.ABM.Environment.adhSiteDens", "total adhesion site density on tumor vasculature, molecule/um^3", "pos" },
	//pharmacokinetics	
	{ "Param.ABM.Pharmacokinetics.nivoDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.nivoDose", "mole/m^3", "pos" },
	{ "Param.ABM.Pharmacokinetics.ipiDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.ipiDose", "mole/m^3", "pos" },
	{ "Param.ABM.Pharmacokinetics.caboDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.caboDose", "mole/m^3", "pos" },
	//T cell
	{ "Param.ABM.TCell.lifespanSD", "days", "pos" },
	{ "Param.ABM.TCell.moveSteps", "", "pos" },
	{ "Param.ABM.TCell.IL2_release_time", "amount of time to release IL2 after stimulation, sec", "pos" },
	{ "Param.ABM.TCell.IL2_prolif_th", "accumulative IL2 exposure to proliferate, sec*ng/mL", "pos" },
	{ "Param.ABM.TCell.IFNg_recruit_Half", "Half-Maximal IFNg level for T cell recruitment, ng/mL", "pos" },
	// TCD4
	{ "Param.ABM.TCD4.TGFB_release_time", "amount of time to release TGFB after stimulation, sec", "pos"  },
	// MDSC
	{ "Param.ABM.MDSC.lifespanSD", "", "pos" },
	{ "Param.ABM.MDSC.moveSteps", "", "pos" },
	//Mac
	{ "Param.ABM.Mac.lifespanSD", "", "pos" },
	{ "Param.ABM.Mac.moveSteps", "", "pos" },
	//cancer cell
	{"Param.ABM.CancerCell.asymmetricDivProb", "", "pr"},
	{"Param.ABM.CancerCell.progGrowthRate", "per day", "pos"},
	{"Param.ABM.CancerCell.senescentDeathRate", "per day", "pos"},
	{"Param.ABM.CancerCell.moveSteps", "", "pos"},
	{"Param.ABM.CancerCell.moveSteps_csc", "", "pos"},
	{"Param.ABM.CancerCell.Tkill_scaler", "", "pos"},
	{"Param.ABM.CancerCell.mincc", "", "pos"},
	{"Param.ABM.CancerCell.C1_CD47", "CD47 expression on cancer cells, molecule/micrometer^2", "pos"},
	{"Param.ABM.CancerCell.IFNgUptake", "per sec", "pos"},
	{"Param.ABM.CancerCell.hypoxia_th", "ng/mL", "pos"},
	{"Param.ABM.CancerCell.density_csc", "cancer cell IC density in core", "prob"},
	//Vas
	{"Param.ABM.Vas.maxPerVoxel", "", "pos"},
	{"Param.ABM.Vas.vas_50", "half maximum of VEGFA of induction of angiogensis pg/ml", "pos"},
	{"Param.ABM.Vas.O2_conc", "", "pos"},
	{"Param.ABM.Vas.Rc", "radius of tumor capillary", "pos"},
	{"Param.ABM.Vas.sigma", "dimensionless ratio of intracapillary to extracapillary transport resistance from Sharan et al; D*alpha / k*Rc = 0.84", "pos"},
	//reference: Assessment of Blood Flow in Hepatocellular Carcinoma Correlations of Computed Tomography Perfusion Imagingan Circulating Angiogenic Factors
	//estimated: 20 / 0.62 mm^2
	{"Param.ABM.Vas.ref_vas_frac", "reference fraction of endotheilal cell per voxel", "pos"},
	{"Param.ABM.Vas.init_density", "inital density of endothelial cell", "pos"},
	{"Param.ABM.Vas.tumble", "tumbling rate due to chemotaxis", "pos"},
	{"Param.ABM.Vas.delta", "tumbling rate due to chemotaxis", "pos"},
	{"Param.ABM.Vas.branch_prob", "branching probability of phelanx cell", "pos"},
	{"Param.ABM.Vas.min_neighbor", "minimal distance to other vessel branch to generate a new tip cell", "pos"},
	// fibroblast
	{ "Param.ABM.Fib.moveSteps","", "pos" },
	//agent chemokine interaction
	{"Param.ABM.cell.PDL1_th", "percent of max PDL1_syn to be detectable", "prob"},
	{"Param.ABM.cell.IFNg_PDL1_half", "c of IFNg to induce PDL1 to half maximal level, ng/mL", "pos"},
	{"Param.ABM.cell.IFNg_PDL1_n", "hill coef for PDL1 expression", "pos"},
	/* molecular level */
	// diffusion grid
	{"Param.Molecular.biofvm.IFNg.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IFNg.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.IFNg.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.IFNg.molecularWeight","1/sec", "pos"},
	{"Param.Molecular.biofvm.IL_2.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IL_2.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.IL_2.uptake","1/s", "pos"},
	{"Param.Molecular.biofvm.IL_2.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.CCL2.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.CCL2.release","moles/day, one cell", "pos"},
	{"Param.Molecular.biofvm.CCL2.uptake","1/s", "pos"},
	{"Param.Molecular.biofvm.CCL2.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.CCL2.molecularWeight","kDa", "pos"},		
	{"Param.Molecular.biofvm.ArgI.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.ArgI.release","mU * Liter/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.ArgI.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.ArgI.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.NO.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.NO.release","mg/sec", "pos"},
	{"Param.Molecular.biofvm.NO.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.NO.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.TGFB.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.TGFB.release.Mac","moles/day, one cell", "pos"},
	{"Param.Molecular.biofvm.TGFB.release.Treg","moles/day, one cell", "pos"},
	{"Param.Molecular.biofvm.TGFB.release.CancerStem","mg/sec", "pos"},
    {"Param.Molecular.biofvm.TGFB.release.CancerProgenitor","mg/sec", "pos"},
	{"Param.Molecular.biofvm.TGFB.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.TGFB.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.IL10.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IL10.release.Treg","mg/sec", "pos"},
	{"Param.Molecular.biofvm.IL10.release.Mac","ng/sec", "pos"},
	{"Param.Molecular.biofvm.IL10.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.IL10.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.IL12.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IL12.release","ng/sec", "pos"},
	{"Param.Molecular.biofvm.IL12.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.IL12.molecularWeight","kDa", "pos"},
	{ "Param.Molecular.biofvm.VEGFA.diffusivity","cm^2/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.release.Mac","ng/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.release.CancerStem", "ng/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.release.CancerProgenitor", "ng/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.uptake", "1/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.decayRate","1/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.molecularWeight","kDa", "pos" },
	{ "Param.Molecular.biofvm.O2.diffusivity","cm^2/sec", "pos" },
	{ "Param.Molecular.biofvm.O2.uptake", "mg/sec", "pos" },
	{ "Param.Molecular.biofvm.O2.decayRate","1/sec", "pos" },
	{ "Param.Molecular.biofvm.O2.molecularWeight","kDa", "pos" },

	//------------------------ int ----------------------------//
	{"Param.ABM.Environment.Tumor.XSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.YSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.ZSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.VoxelSize", "voxel resolution, microns", "pos"},
	{"Param.ABM.Environment.Tumor.nr_T_voxel", "", "pos"},
	{"Param.ABM.Environment.Tumor.nr_T_voxel_C", "", "pos"},
	{"Param.ABM.Environment.Tumor.stem_mode", "", "pos" },
	{"Param.ABM.Environment.ShuffleInterval", "", "pos"},
	{"Param.ABM.Environment.gridshiftInterval", "", "pos"},
	{"Param.ABM.TCell.div_interval", "", "pos"},
	{"Param.ABM.TCell.div_limit", "", "pos"},
	{"Param.ABM.TCD4.div_interval", "", "pos"},
	{"Param.ABM.TCD4.div_limit", "", "pos"},		
	{"Param.ABM.CancerCell.progenitorDivMax", "", "pos"},
	{"Param.Molecular.stepPerSlice","", "pos"},

	// ---------------------- bool -----------------------------//
	{"Param.QSP.simulation.use_resection", "", "" },
	{"Param.Molecular.allMolecularOff", "", ""},
	{"Param.Molecular.diffusionOff", "", ""},
	{"Param.Molecular.cellOdeOff", "", ""},
	{"Param.ABM.Pharmacokinetics.nivoOn", "", "" },
	{"Param.ABM.Pharmacokinetics.ipiOn", "", "" },
	{"Param.ABM.Pharmacokinetics.entOn", "", "" },
	{"Param.ABM.Pharmacokinetics.caboOn", "", "" },

};

Param::Param()
	:ParamBase()
{
	setupParam();
}

/*! Setup parameter storage
	instantiation of pure virtual member of the base class.
	setup description vector;
	initialize parameter value vectors with 0/false, 
	with size determined by enums, 
	so that other base class members can access vector sizes
*/
void Param::setupParam(){

	size_t nrExternalParam = PARAM_FLOAT_COUNT + PARAM_INT_COUNT + PARAM_BOOL_COUNT;
	for (size_t i = 0; i < nrExternalParam; i++)
	{
		_paramDesc.push_back(std::vector<std::string>(_description[i], 
			_description[i]+ PARAM_DESCRIPTION_FIELD_COUNT));
	}
	_paramFloat = std::vector<double>(PARAM_FLOAT_COUNT, 0);
	_paramInt= std::vector<int>(PARAM_INT_COUNT, 0);
	_paramBool= std::vector<bool>(PARAM_BOOL_COUNT, false);
	_paramFloatInternal = std::vector<double>(PARAM_FLOAT_INTERNAL_COUNT, 0);
	_paramIntInternal = std::vector<int>(PARAM_INT_INTERNAL_COUNT, 0);
	_paramBoolInternal = std::vector<bool>(PARAM_BOOL_INTERNAL_COUNT, false);
}

/*! Calculate internal parameters
*/
void Param::processInternalParams(){
	
	_paramFloatInternal[PARAM_AVOGADROS] = AVOGADROS;

	//micrometer to cm
	_paramFloatInternal[PARAM_VOXEL_SIZE_CM] = _paramInt[PARAM_VOXEL_SIZE] / 1e4;

	/*
	_paramFloatInternal[PARAM_T_CELL_LIFE_MEAN_SLICE] = _paramFloat[PARAM_T_CELL_LIFE_MEAN]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;
	*/
	_paramFloatInternal[PARAM_T_CELL_LIFE_SD_SLICE] = _paramFloat[PARAM_T_CELL_LIFE_SD]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;

	_paramBoolInternal[PARAM_MOLECULAR_MODULES_ON]
		= !getVal(PARAM_ALL_MOLECULAR_OFF);

	_paramBoolInternal[PARAM_DIFFUSION_ON] 
		= getVal(PARAM_MOLECULAR_MODULES_ON) && !getVal(PARAM_DIFFUSION_OFF);
	
	_paramBoolInternal[PARAM_CELL_ODE_ON] 
		= getVal(PARAM_MOLECULAR_MODULES_ON) && !getVal(PARAM_ALL_CELL_ODE_OFF);

}

//! update from QSP parameters
void Param::update_from_qsp(void){

	// cell (need to be 1,  not in mole)
	_paramFloatInternal[PARAM_CELL] = QP(10) * AVOGADROS;

	// volume of cancer cells, vol_cell
	_paramFloatInternal[PARAM_V_CELL] = QP(12);

	// volume of T cells, vol_Tcell
	_paramFloatInternal[PARAM_V_TCELL] = QP(13);

	// surface area of T cells, area_Tcell
	_paramFloatInternal[PARAM_A_TCELL] = QP(94);

	// initial tumor volume calculated by initial tumor diameter.
	_paramFloatInternal[PARAM_INIT_TUM_VOL] = (PI * std::pow(QP(17), 3)) / 6;

	_paramFloatInternal[PARAM_INIT_TUM_DIAM] = QP(17);

	std::cout << "cell: " << _paramFloatInternal[PARAM_CELL]
		<< ", minimum tumor size : " << _paramFloatInternal[PARAM_VT_MIN]
		<< ", cancer cell size: " << _paramFloatInternal[PARAM_V_CELL]
		<< ", T cell size: " << _paramFloatInternal[PARAM_V_TCELL] << std::endl;
//  ==========================================================
//|| Treg and MDSC does not have maximum number this version  ||
//  ===========================================================
	// maximum concentration of Tregs per volume in the tumor
	//_paramFloatInternal[PARAM_TREGMAX] = QP(178);		

	// maximum concentration of MDSC per volume
	//_paramFloatInternal[PARAM_MDSCMAX] = QP(177);


	//==============================================
	//  MDSC module are only in ABM  in this version
	//*================================================
	 
	// MDSC base recruitment
	_paramFloatInternal[PARAM_MDSC_RECRUIT_K] = QP(205);	
	//unit: nanomolarity (1e-9 mole/L) -> ng/ml; conversion factor: 1e9 (mole/L to mole/m^3), 1e6 (m^3 to ml)
	_paramFloatInternal[PARAM_MDSC_EC50_CCL2_REC] = QP(204) * _paramFloat[PARAM_CCL2_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;

	// half maximal inhibitory concentration of NO on inhibition of CD8+ T cell cytotoxic activity (ng/ml)
	_paramFloatInternal[PARAM_MDSC_IC50_NO_CTL] = QP(212) * _paramFloat[PARAM_NO_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;

	// half maximal inhibitory concentration of Arg I on inhibition of CD8+ T cell cytotoxic activity (ng/ml) 
	_paramFloatInternal[PARAM_MDSC_IC50_ArgI_CTL] = QP(211) * _paramFloat[PARAM_ARGI_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;

	// half maximal effective concentration of arginase I on Treg expansion (ng/ml) (molecular weights in kDa
	_paramFloatInternal[PARAM_MDSC_EC50_ArgI_Treg] = QP(213)  * _paramFloat[PARAM_ARGI_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
	
	std::cout << ", PARAM_MDSC_IC50_NO_CTL: " << _paramFloatInternal[PARAM_MDSC_IC50_NO_CTL]
		<< ", PARAM_MDSC_IC50_ArgI_CTL: " << _paramFloatInternal[PARAM_MDSC_IC50_ArgI_CTL]
		<< ", PARAM_MDSC_EC50_ArgI_Treg: " << _paramFloatInternal[PARAM_MDSC_EC50_ArgI_Treg] << std::endl;

	// half maximal inhibitory concentration of entinostat on anti-proliferation of tumor cell (mol/m^3)
	//_paramFloatInternal[PARAM_IC50_ENT_C] = QP(164);	

	// half maximal inhibitory concentration of entinostat on inhibition of CCL2 production (mol/m^3)
	//_paramFloatInternal[PARAM_IC50_ENT_CCL2] = QP(179);	

	// half maximal inhibitory concentration of entinostat on inhibition of NO production (mol/m^3)
	//_paramFloatInternal[PARAM_IC50_ENT_NO] = QP(171);

	// half maximal inhibitory concentration of entinostat on inhibition of Arg I production (mol/m^3)
	//_paramFloatInternal[PARAM_IC50_ENT_ARGI] = QP(181);	
	
	// Area of synapse
	_paramFloatInternal[PARAM_A_SYN] = QP(93);
	// Area of Cancer cell
	_paramFloatInternal[PARAM_A_CELL] = QP(95);
	// Area of T cell
	_paramFloatInternal[PARAM_A_TCELL] = QP(94);
	// Area of Macrophages
	_paramFloatInternal[PARAM_A_MAC] = QP(258);
	// number of PD1/PDL1 binding for half maximal inhibition, PD1_50
	_paramFloatInternal[PARAM_PD1_PDL1_HALF] = QP(175) * _paramFloatInternal[PARAM_A_CELL] * 3;

	// total number of PD1 per synapse on T cell = T_PD1_total*A_syn/A_Tcell; T_PD1_total in density (molecule per micrometer^2)
	_paramFloatInternal[PARAM_PD1_SYN] = QP(179);  /// _paramFloatInternal[PARAM_A_TCELL];
	// total number of PD1 per synapse on Macrophages = M_PD1_total*A_syn; T_PD1_total in density (molecule per micrometer^2)
	// _paramFloat[PARAM_MAC_PD1_AREA] in molecule / micrometer^2, 1m^2 = 1e12 micrometer^2
	_paramFloatInternal[PARAM_MAC_PD1_SYN] = QP(256); ///	 _paramFloatInternal[PARAM_A_MAC];

	// _paramFloat[PARAM_C1_CD47_SYN] in molecule / micrometer^2, 1m^2 = 1e12 micrometer^2
	_paramFloatInternal[PARAM_C1_CD47_SYN] = QP(255);

	// _paramFloat[PARAM_MAC_SIRPa] in molecule / micrometer^2, mole/m^2 1m^2 = 1e12 micrometer^2
	_paramFloatInternal[PARAM_MAC_SIRPa_SYN] = QP(257);

	// Number of PDL1 per synapse = C1_PDL1_total; C1_PDL1_total in density (molecule per micrometer^2)
	_paramFloatInternal[PARAM_PDL1_SYN_MAX] = QP(183); 
	
	// total number of PDL1 per Treg cell = T_PDL1_total
	// _paramFloatInternal[PARAM_PDL1_CELL] = QP(182);

	// total number of PDL1 per cancer cell = C1_PDL1_total
	_paramFloatInternal[PARAM_PDL1_CELL] = QP(183);

	// k1 for PDL1-PD1 calculation =  kon_PD1_PDL1 / (koff_PD1_PDL1* A_CELL)
	_paramFloatInternal[PARAM_PDL1_K1] = QP(147) / (QP(161) * _paramFloatInternal[PARAM_A_CELL]);

	// k2 for PDL1-aPD1 calculation = 2* kon_PD1_aPD1 / (koff_PD1_aPD1 * gamma_T_nivo)
	_paramFloatInternal[PARAM_PDL1_K2] = 2 * QP(152)  / (QP(163) * QP(125));

	// k3 for PDL1-PD1 calculation = (Chi_PD1 * kon_PD1_aPD1) / (2 * koff_PD1_aPD1)
	_paramFloatInternal[PARAM_PDL1_K3] = QP(172) * QP(152) / (2 * QP(163));
	// hill coefficient
	_paramFloatInternal[PARAM_N_PD1_PDL1] = QP(176);

	
	std::cout << "" 
		<< "k1: " << _paramFloatInternal[PARAM_PDL1_K1]
		<< ", k2: " << _paramFloatInternal[PARAM_PDL1_K2]
		<< ", k3: " << _paramFloatInternal[PARAM_PDL1_K3]
		<< ", PD1 per synapse on T cell: " << _paramFloatInternal[PARAM_PD1_SYN]
		<< ", PD1 per synapse on Mac: " << _paramFloatInternal[PARAM_MAC_PD1_SYN]
		<< ", PDL1 per synapse on Cancer cell: " << _paramFloatInternal[PARAM_PDL1_SYN_MAX]
		<< std::endl;
	std::cout << "k50: " << _paramFloatInternal[PARAM_PD1_PDL1_HALF]<< std::endl;
	
	// Binding rate between kon_CTLA4 and ipi
	//_paramFloatInternal[PARAM_KON_CTLA4_IPI] = QP(24);

	// Unbinding rate between CTLA4 and ipi
	//_paramFloatInternal[PARAM_KOFF_CTLA4_IPI] = QP(150);
	_paramFloatInternal[PARAM_KOFF_CTLA4_IPI] = 6.96e-06;

	// Volume fraction available to ipi in tumor compartment
	//_paramFloatInternal[PARAM_GAMMA_T_IPI] = QP(121);
	_paramFloatInternal[PARAM_GAMMA_T_IPI] = 0.718;

	// Antibody cross-arm binding efficiency  that also includes the conversion of kon from 3D to 2D (estimated)
	//_paramFloatInternal[PARAM_CHI_CTLA4_IPI] = QP(151);
	_paramFloatInternal[PARAM_CHI_CTLA4_IPI] = 3;

	// CTLA4 occupancy for half-maximal Treg inactivation by macrophages (estimated)  = CTLA_50 * A_treg
	_paramFloatInternal[PARAM_TREG_CTLA4_50] = QP(163) * QP(104);

	// total number of CTLA4 per synapse = Treg_CTLA4_total * A_treg; T_CTLA_total in density (molecule per micrometer^2)
	_paramFloatInternal[PARAM_CTLA4_TREG] = QP(140) * QP(104);

	// Anti-CTLA4 ADCC (antibody-dependent cellular cytotoxicity) rate of Treg (Richards 2008, PMID: 18723496)
	//_paramFloatInternal[PARAM_K_ADCC] = QP(159) * _paramFloat[PARAM_SEC_PER_TIME_SLICE];
	_paramFloatInternal[PARAM_K_ADCC] = 0.1 * _paramFloat[PARAM_SEC_PER_TIME_SLICE];
	/*
	std::cout << "kon_CTLA_ipi: " << _paramFloatInternal[PARAM_KON_CTLA4_IPI]
		<< ", koff_CTLA_ipi: " << _paramFloatInternal[PARAM_KOFF_CTLA4_IPI]
		<< ", gamma_T_ipi: " << _paramFloatInternal[PARAM_GAMMA_T_IPI]
		<< ", chi_CTLA_ipi: " << _paramFloatInternal[PARAM_CHI_CTLA4_IPI]
		<< ", CTLA4 total in Treg: " << _paramFloatInternal[PARAM_CTLA4_TREG] << std::endl;
	*/
	// Parameters calculated from QSP parameter values
	double t_step_sec = _paramFloat[PARAM_SEC_PER_TIME_SLICE];

	// Update cabozantinib module from QSP model (values are temperary)
	// IC50	for receptors inhibited by cabozantinib
	_paramFloatInternal[PARAM_IC50_AXL] = QP(231);
	_paramFloatInternal[PARAM_IC50_VEGFR2] = QP(232);
	_paramFloatInternal[PARAM_IC50_MET] = QP(229);
	_paramFloatInternal[PARAM_IC50_RET] = QP(230);

	// theraputic effects parameters elicited by cabozantinib
	_paramFloatInternal[PARAM_LAMBDA_C_CABO] = QP(234);

	// poisson process is calculated in TCD4.cpp, here is the rate of the process
	// T cell killing of Cancer cell
	// QP(48): k_C_death_by_T (day^-1, sec^-1 internal)
	// Becuase the ABM contain modules that QSP, which contains extra immunesuppresive modules
	_paramFloatInternal[PARAM_ESCAPE_BASE] = std::exp(-t_step_sec * QP(70));
	// Macrophage killing of Cancer cell
	//_paramFloat[PARAM_MAC_K_M1_PHAGO] k_C_death_by_Macrophage (day^-1, sec^-1 internal)
	// Becuase the ABM contain modules that QSP, which contains extra immunesuppresive modules
	// 5 is a arbitary coefficient to compensate immuno-suppresive module not present in QSP model
	_paramFloatInternal[PARAM_ESCAPE_MAC_BASE] = std::exp(-t_step_sec * QP(249));

	// T cell exhaustion from PDL1; 
	_paramFloatInternal[PARAM_EXHUAST_BASE_PDL1] = std::exp(-t_step_sec * QP(69));

	// T cell exhaustion from Treg inhibition, k_Treg;
	_paramFloatInternal[PARAM_EXHUAST_BASE_TREG] = std::exp(-t_step_sec * QP(51));

	// rate of Th to Treg transformation , units: 1/days -> 1/timestep
	// poisson process is calculated in TCD4.cpp, here is the rate of the process
	_paramFloatInternal[PARAM_K_TH_TREG] = t_step_sec * QP(192);
	// Macrophage module related parameter covnersion
	// Rate of M1 to M2 macrophage polarization, units: 1/days -> 1/timestep

	_paramFloatInternal[PARAM_MAC_M2_POL] = t_step_sec * QP(244);
	// Rate of M2 to M1 macrophage polarization, units: 1/days -> 1/timestep
	_paramFloatInternal[PARAM_MAC_M1_POL] = t_step_sec * QP(245);
	//unit: nanomolarity (1e-9 mole/L) -> ng/ml; conversion factor: 1e9 (mole/L to mole/m^3), 1e6 (m^3 to ml)
	_paramFloatInternal[PARAM_MAC_EC50_CCL2_REC] = QP(204) * _paramFloat[PARAM_CCL2_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
	// TGFb_50 Half-Maximal TGFb level for Th-to-Treg differentiation / chemoresistance development / M1-to-M2 polarization
    //unit: nanomolarity (1e-9 mole/L) -> ng/ml; conversion factor: 1e9 (mole/L to mole/m^3), 1e6 (m^3 to ml)
	_paramFloatInternal[PARAM_MAC_TGFB_EC50] = QP(195) * _paramFloat[PARAM_TGFB_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
	//unit: nanomolarity (1e-9 mole/L) -> ng/ml
	_paramFloatInternal[PARAM_TEFF_TGFB_EC50] = QP(196) * _paramFloat[PARAM_TGFB_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
	//unit: picomolarity (1e-12 mole/L) -> ng/ml
	_paramFloatInternal[PARAM_MAC_IL_10_EC50] = QP(246) * _paramFloat[PARAM_IL_10_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
	//unit: picomolarity (1e-12 mole/L) -> ng/ml
	_paramFloatInternal[PARAM_MAC_IL_10_HALF_PHAGO] = QP(259) * _paramFloat[PARAM_IL_10_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
	//unit: picomolarity (1e-12 mole/L) -> ng/ml
	_paramFloatInternal[PARAM_MAC_IFN_G_EC50] = QP(248) * _paramFloat[PARAM_IFN_G_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
    //unit: picomolarity (1e-12 mole/L) -> ng/ml
	_paramFloatInternal[PARAM_MAC_IL_12_EC50] = QP(247) * _paramFloat[PARAM_IL_12_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
	//unit: moles/m^2
	_paramFloatInternal[PARAM_MAC_SIRPa_HALF] = QP(253);
	// 1/(micromolarity*minute*nanometer) ->   0.001 (micromolarity to mole/m^3), 60 (minute to second), 1e-9 (nanometer to meter)
	_paramFloatInternal[PARAM_KON_SIRPa_CD47] = QP(251);
	// 1/minute
	_paramFloatInternal[PARAM_KOFF_SIRPa_CD47] = QP(252);
	// hill coefficient
	_paramFloatInternal[PARAM_N_SIRPa_CD47] = QP(254);
	std::cout << "M1 to M2 Mac polarization: " << _paramFloatInternal[PARAM_MAC_M2_POL]
		<< ", M2 to M1 Mac polarization: " << _paramFloatInternal[PARAM_MAC_M1_POL]
		<< ", EC50 TGFB: " << _paramFloatInternal[PARAM_MAC_TGFB_EC50]
		<< ", PARAM_TEFF_TGFB_EC50: " << _paramFloatInternal[PARAM_TEFF_TGFB_EC50]
		<< ", PARAM_MAC_IL_10_EC50: " << _paramFloatInternal[PARAM_MAC_IL_10_EC50]
		<< ", PARAM_MAC_IFN_G_EC50: " << _paramFloatInternal[PARAM_MAC_IFN_G_EC50] << std::endl;

	// ECM secretion per time step for both fibroblast and caf (parameter inputs are ng/s) 6970071747.68519 is the conversion factor to nanomole/cell/day
	_paramFloatInternal[PARAM_FIB_ECM_RELEASE_FIB] = QP(264)  / 6970071747.68519;
	_paramFloatInternal[PARAM_FIB_ECM_RELEASE_CAF] = QP(265)  / 6970071747.68519;

	//fibroblast activation rate due to tgfb (1/day)
	_paramFloatInternal[PARAM_FIB_CAF_ACTIVATION] = QP(263) * t_step_sec;
	_paramFloatInternal[PARAM_FIB_CAF_EC50] = QP(195) * _paramFloat[PARAM_TGFB_MOLECULAR_WEIGHT] * 1e3 * 1e9 / 1e6;
	_paramFloatInternal[PARAM_FIB_ECM_TGFB_FACTOR] = 2;

	_paramFloatInternal[PARAM_FIB_ECM_DECAY_RATE] = QP(266);
	// unit: mole/m^3 -> nmole/ml
	_paramFloatInternal[PARAM_FIB_ECM_BASELINE] = QP(267) * 1e3;
	_paramFloatInternal[PARAM_FIB_ECM_SATURATION] = QP(268) * 1e3;
	_paramFloatInternal[PARAM_FIB_ECM_MOT_EC50] = QP(274) * 5e3;

	std::cout << "normal fib ECM release rate: " << _paramFloatInternal[PARAM_FIB_ECM_RELEASE_FIB]
		<< ", normal caf ECM release rate: " << _paramFloatInternal[PARAM_FIB_ECM_RELEASE_CAF]
		<< ", Firboblast activation rate: " << _paramFloatInternal[PARAM_FIB_CAF_ACTIVATION]
		<< ", PARAM_CAF_TGFB_EC50: " << _paramFloatInternal[PARAM_FIB_CAF_EC50]
		<< ", PARAM ECM tgfb factor: " << _paramFloatInternal[PARAM_FIB_ECM_TGFB_FACTOR] 
		<< ", PARAM_FIB_ECM_DECAY_RATE: " << _paramFloatInternal[PARAM_FIB_ECM_DECAY_RATE]
		<< ", PARAM_FIB_ECM_BASELINE: " << _paramFloatInternal[PARAM_FIB_ECM_BASELINE]
		<< ", PARAM_FIB_ECM_SATURATION: " << _paramFloatInternal[PARAM_FIB_ECM_SATURATION]
		<< ", PARAM_FIB_ECM_RELEASE_FIB: " << _paramFloatInternal[PARAM_FIB_ECM_RELEASE_FIB]
		<< ", PARAM_FIB_ECM_RELEASE_CAF: " << _paramFloatInternal[PARAM_FIB_ECM_RELEASE_CAF]
		<< ", PARAM_FIB_ECM_MOT_EC50: " << _paramFloatInternal[PARAM_FIB_ECM_MOT_EC50]
		<< std::endl;

	// time for resection
	_paramFloatInternal[PARAM_RESECT_TIME_STEP] = _paramFloat[PARAM_QSP_T_RESECTION] * SEC_PER_DAY / t_step_sec;


	// Recruitment
	// T effector recruitment
	/* for each mole of adhesion site, the amount of T cell recruited (in unit of mole):
	dt * q_T1_T_in * Cent.T * V_T
	# Weighted QSP:
	# Units:
	All parameters come in SI units.
	Cent.T should use SI units --- in mole

	The result is mole recruited per mole site, or number per site, so no conversion needed.
	When calculating recruitment probability:
	p = 1 / (s * m^3) * cell  * (m^3 / cell) * dt
	*/

	//The number of adhesion site per voxel is:
	double site_per_voxel = _paramFloat[PARAM_ADH_SITE_DENSITY] * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]), 3);
	//number of adhesion sites needed to recruit a single cell (becoming a recruitment port)
	double site_per_port = _paramFloat[PARAM_REC_SITE_FACTOR];
	//how many port per voxel have
	_paramFloatInternal[PARAM_REC_PORT] = site_per_voxel / site_per_port;
	std::cout << "PARAM_REC_PORT_PROB : " << _paramFloatInternal[PARAM_REC_PORT] << std::endl;
	/*When calculating recruitment probability:
	*/
	double  w = _paramFloat[PARAM_WEIGHT_QSP];
	// Teff -> k (1/mol) // p = k (1/mol) * Cent.T (mol), q_T1_T_in
	//_paramFloatInternal[PARAM_TEFF_RECRUIT_K] = QP(36) * site_per_port * t_step_sec  / w / _paramFloat[PARAM_ADH_SITE_DENSITY];

	_paramFloatInternal[PARAM_TEFF_RECRUIT_K] = QP(62) *  t_step_sec * AVOGADROS * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]) / 1e6, 3) * _paramFloatInternal[PARAM_REC_PORT];
	// TCD4 -> k (1/mol) // p = k (1/mol) * Cent.T (mol), q_T0_T_in
	//_paramFloatInternal[PARAM_TCD4_RECRUIT_K] = QP(59) * site_per_port * t_step_sec  / w / _paramFloat[PARAM_ADH_SITE_DENSITY];
	_paramFloatInternal[PARAM_TREG_RECRUIT_K] = QP(35) * t_step_sec * AVOGADROS * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]) / 1e6, 3) * _paramFloatInternal[PARAM_REC_PORT];
	_paramFloatInternal[PARAM_TH_RECRUIT_K] = QP(35) * t_step_sec * AVOGADROS * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]) / 1e6, 3) * _paramFloatInternal[PARAM_REC_PORT];
	std::cout << "keff: " << _paramFloatInternal[PARAM_TEFF_RECRUIT_K] << std::endl;
	std::cout << "kreg: " << _paramFloatInternal[PARAM_TREG_RECRUIT_K] << std::endl;
	std::cout << "kth: " << _paramFloatInternal[PARAM_TH_RECRUIT_K] << std::endl;

	// The recruitment of MDSC does not depend on central compartment MDSC (central compartment does not have MDSC)
	// so the recruitment equation is different from Teff and Treg
	// 11000 * cells/(ml*day) -> cells / (m^3 * s)
	// 1 m^3 = 1e6 ml
	_paramFloatInternal[PARAM_MDSC_RECRUIT_K] = QP(205) * t_step_sec * AVOGADROS * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]) / 1e6, 3);;
	_paramFloatInternal[PARAM_MAC_RECRUIT_K] = QP(235) * t_step_sec * AVOGADROS * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]) / 1e6, 3);

	std::cout << "MDSC_RECRUIT_K: " << _paramFloatInternal[PARAM_MDSC_RECRUIT_K] << ", MAC_RECRUIT_K: " << _paramFloatInternal[PARAM_MAC_RECRUIT_K] << std::endl;
	// APC -> k (m^3/mol) // p = k (m^3/mol) * (APC0_T*V_T-V_T.APC)
	_paramFloatInternal[PARAM_APC_RECRUIT_K] = 0;
	
	// APC density in the tumour
	_paramFloatInternal[PARAM_APC0_T] = QP(79);

	// APC transmigration rate from tumor to lymph node
	_paramFloatInternal[PARAM_APC_TRANSMIG] = QP(76);

	// Number of T0 cell Clonality to lymph node
	_paramFloatInternal[PARAM_T0_CLONE] = QP(26);

	// Number of T1 cell Clonality to lymph node
	_paramFloatInternal[PARAM_T1_CLONE] = QP(53);

	// antigen concentration in cancer cell
	_paramFloatInternal[PARAM_ANTIGEN_PER_CELL] = QP(110);

	// antigen uptake rate by mAPC
	_paramFloatInternal[PARAM_ANTIGEN_UP] = QP(84);

	// antigen degrdation rate in the tumor
	_paramFloatInternal[PARAM_K_xP_DEG] = QP(85);
	
	// mean life of Tcell, unit: time step, 
	// which is different from the QSP model.
	_paramFloatInternal[PARAM_T_CELL_LIFE_MEAN_SLICE] = 1 / QP(59) / t_step_sec / 5;
	// mean life of TCD4, unit: time step
	_paramFloatInternal[PARAM_TCD4_LIFE_MEAN] = 1 / QP(32) / t_step_sec / 5;
	std::cout << "T cell life mean: " << _paramFloatInternal[PARAM_T_CELL_LIFE_MEAN_SLICE]
		<< ", TCD4_LIFE_MEAN: " << _paramFloatInternal[PARAM_TCD4_LIFE_MEAN] << std::endl;
	// mean life of MDSC, unit: time step
	_paramFloatInternal[PARAM_MDSC_LIFE_MEAN] = 1 / QP(206) / t_step_sec;
	// mean life of MAC, unit: time step
	_paramFloatInternal[PARAM_MAC_LIFE_MEAN] = 1 / QP(236) / t_step_sec;
	// mean life of APC, unit: time step
	_paramFloatInternal[PARAM_APC_LIFE_MEAN] = 1 / QP(77) / t_step_sec;
	std::cout << "PARAM_APC_LIFE_MEAN: " << QP(77)  <<  ", "  << getVal(PARAM_APC_LIFE_MEAN) << std::endl;
	
	// Maximum rate of APC maturation
	_paramFloatInternal[PARAM_K_APC_MAT] = QP(75);
	std::cout
	<< "PARAM_ESCAPE_BASE, " << _paramFloatInternal[PARAM_ESCAPE_BASE] << "\n"
	<< "PARAM_ESCAPE_MAC_BASE, " << _paramFloatInternal[PARAM_ESCAPE_MAC_BASE] << "\n"
	<< "PARAM_K_TH_TREG, " << _paramFloatInternal[PARAM_K_TH_TREG] << "\n"
	<< "PARAM_EXHUAST_BASE_PDL1, " << _paramFloatInternal[PARAM_EXHUAST_BASE_PDL1] << "\n"
	<< "PARAM_EXHUAST_BASE_TREG, " << _paramFloatInternal[PARAM_EXHUAST_BASE_TREG] << "\n"
	<< std::endl;
	

	/*Cancer cell dynamics parameters*/

	// stem cell division rate is calculated from QSP parameter
	// unit: s^-1, k_C1_growth
	double rs = QP(15) / (1 - _paramFloat[PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB]);
	// unit: day^-1
	_paramFloatInternal[PARAM_CSC_GROWTH_RATE] = rs * SEC_PER_DAY;

	_paramFloatInternal[PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE]
		= std::log(2)/rs / getVal(PARAM_SEC_PER_TIME_SLICE);

	_paramFloatInternal[PARAM_CANCER_SENESCENT_MEAN_LIFE] =
		1 / _paramFloat[PARAM_CANCER_SENESCENT_DEATH_RATE]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;
		
	
	_paramFloatInternal[PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE]
		= std::log(2)/ getVal(PARAM_CANCER_PROG_GROWTH_RATE) 
		* SEC_PER_DAY / getVal(PARAM_SEC_PER_TIME_SLICE) + .5;
	
	//_paramFloatInternal[PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE]
		//= std::log(2) / rs / getVal(PARAM_SEC_PER_TIME_SLICE);
	
	std::cout << "PARAM_CSC_GROWTH_RATE: "<< getVal(PARAM_CSC_GROWTH_RATE) << ", k_C1_growth: " << QP(15)
		<< ", PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE: " << _paramFloatInternal[PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE] << std::endl;
	/*
	std::cout << getVal(PARAM_INT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE)
		<< ", "<< getVal(PARAM_CANCER_SENESCENT_MEAN_LIFE)
		<< ", "<< getVal(PARAM_INT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE)
		<< std::endl;
	*/

	return;
}

};
};