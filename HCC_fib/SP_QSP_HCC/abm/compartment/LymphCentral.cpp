#include "LymphCentral.h"
#include "Tumor.h"
#include "../../core/GlobalUtilities.h"

// shorthands
// get raw value (original units)
#define GET_PARAM_RAW(x) _QSP_model.getSystem()->getParameterVal(x, true)
#define SET_PARAM_RAW(x, y) _QSP_model.getSystem()->setParameterVal(x, y, true)
#define GET_VAR_RAW(x) _QSP_model.getSystem()->getSpeciesVar(x, true)
#define SET_VAR_RAW(x, y) _QSP_model.getSystem()->setSpeciesVar(x, y, true)
// get value (SI units)
#define GET_PARAM(x) _QSP_model.getSystem()->getParameterVal(x, false)
#define SET_PARAM(x, y) _QSP_model.getSystem()->setParameterVal(x, y, false)
#define GET_VAR(x) _QSP_model.getSystem()->getSpeciesVar(x, false)
#define SET_VAR(x, y) _QSP_model.getSystem()->setSpeciesVar(x, y, false)
// parameter (SI units)
#define QSP_CONST(x) LymphBloodQSP::get_class_param(x)

// indices of parameter/variables in their vectors
// y
#define QSP_ID_TUM_C1 23
#define QSP_ID_TUM_C2 26
#define QSP_ID_CENT_TREG 1
#define QSP_ID_CENT_TEFF 3
#define QSP_ID_CENT_TH 7
#define QSP_ID_TUM_TREG 27
#define QSP_ID_TUM_TEFF 28
#define QSP_ID_TUM_TH 37 
#define QSP_ID_TUM_MDSC 39
#define QSP_ID_TUM_M1 44
#define QSP_ID_TUM_M2 45
#define QSP_ID_CENT_NIVO 4
#define QSP_ID_TUM_NIVO 34
#define QSP_ID_SITE1_CABO 9
#define QSP_ID_SITE2_CABO 10
#define QSP_ID_TUM_CABO 43
#define QSP_ID_P0 32
#define QSP_ID_P1 33
#define QSP_ID_TUM_CX 20
#define QSP_ID_TUM_TEXH 21
#define QSP_ID_TUM_TH_EXH 22
#define QSP_C_MAX 24
#define QSP_ID_TUM_FIB 48
#define QSP_ID_TUM_CAF 49
#define QSP_ID_TUM_ECM 50
// class_param
#define QSP_V_C 0
#define QSP_CELL 10
#define QSP_P1_C1 91
#define QSP_P0_C1 110
#define QSP_N_T0_CLONES 26
#define QSP_N_T1_CLONES 53
#define QSP_INIT_TUM_DIAM 17
#define QSP_VOL_CELL 12
#define QSP_VOL_TCELL 13
#define QSP_VOL_MCELL 250
#define QSP_VOL_FIBCELL 275
#define QSP_VOL_CAFCELL 276
#define QSP_ECM_MW 269
#define QSP_ECM_DENSITY 270
#define QSP_VE_T 14

// constants
#define AVOGADROS 6.022140857E23 
#define SEC_PER_DAY 86400
#define PI 3.1415926525897932384626

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

LymphCentral::LymphCentral()
	: _QSP_model()
	, _var_qsp_to_abm()
{
	// Cent.Teff, Cent.TCD4, Cent.Nivo
	_var_qsp_to_abm = std::vector<double>(QSPEX_VAR_NUM, 0);
}

LymphCentral::~LymphCentral()
{
}

void LymphCentral::setup_param(LymphBloodParam& p){

	bool steadystate = true;
	unsigned int n = _QSP_model.getSystem()->get_num_variables();
	std::vector<double> ss_val(n, 0);

	if (steadystate)
	{
		QSP ss;
		LymphBloodQSP::_QSP_weight = 1;
		LymphBloodQSP::use_steady_state = true;
		LymphBloodQSP::use_resection = false;
		LymphBloodQSP::setup_class_parameters(p);
		ss.getSystem()->setup_instance_tolerance(p);
		ss.getSystem()->setup_instance_varaibles(p);
		ss.getSystem()->eval_init_assignment();
		// run to steady state or until volume condition is met
		double tss = params.getVal(PARAM_QSP_STEADYSTATE) * SEC_PER_DAY;
		double tt = 0;
		double deltatt = params.getVal(PARAM_SEC_PER_TIME_SLICE);

		
		// First define the auxiliary variables needed for tumor volume calculation
		double AUX_VAR_C_total = ss_val[QSP_ID_TUM_C1] + ss_val[QSP_ID_TUM_C2];  // V_T.C1 + V_T.C2

		double AUX_VAR_T_total = ss_val[QSP_ID_TUM_TREG] + ss_val[QSP_ID_TUM_TEFF] + ss_val[QSP_ID_TUM_TH];  // V_T.T0 + V_T.T1 + V_T.Th

		double AUX_VAR_M_total = ss_val[QSP_ID_TUM_M1] + ss_val[QSP_ID_TUM_M2];  // V_T.Mac_M1 + V_T.Mac_M2

		// Then calculate tumor volume using these auxiliary variables
		double tumor_volume = ((ss_val[20] + AUX_VAR_C_total) / AVOGADROS * QSP_CONST(QSP_VOL_CELL) + 
							(ss_val[21] + ss_val[22] + AUX_VAR_T_total) / AVOGADROS * QSP_CONST(QSP_VOL_TCELL)) / QSP_CONST(QSP_VE_T) + 
							AUX_VAR_M_total / AVOGADROS * QSP_CONST(QSP_VOL_MCELL) / QSP_CONST(QSP_VE_T) + 
							ss_val[48] / AVOGADROS* QSP_CONST(QSP_VOL_FIBCELL) / QSP_CONST(QSP_VE_T) + 
							ss_val[49] / AVOGADROS* QSP_CONST(QSP_VOL_CAFCELL) / QSP_CONST(QSP_VE_T) + 
							ss_val[50] / 1e9 * QSP_CONST(QSP_ECM_MW) / QSP_CONST(QSP_ECM_DENSITY); // ECM in nanomoles, convert to moles

		double presimulation_radius = params.getVal(PARAM_QSP_PRE_SIMULATION_DIAM_FRAC) * QSP_CONST(QSP_INIT_TUM_DIAM);
		std::cout << "presimulation radius: " << QSP_CONST(QSP_INIT_TUM_DIAM) << std::endl;
		std::cout << "central volume: " << QSP_CONST(QSP_V_C) << std::endl;
		double tumor_volume_ref_presimulation = (PI * std::pow(presimulation_radius, 3)) / 6;
		//ss_val[QSP_ID_TUM_C1] < cancer_count_ref
		//double const cancer_count_ref = params.getVal(PARAM_QSP_PRE_SIMULATION_CANCER);
		while (tt < tss && tumor_volume < tumor_volume_ref_presimulation)
		{
			//std::cout << "current tumor volume: " << tumor_volume << ", target tumor volume : " << tumor_volume_ref_presimulation << std::endl;

			ss.solve(tt, deltatt);
			for (size_t i = 0; i < n; i++)
			{
				ss_val[i] = ss.getSystem()->getSpeciesVar(i);
			}

			
			AUX_VAR_C_total = ss_val[QSP_ID_TUM_C1] + ss_val[QSP_ID_TUM_C2];  // V_T.C1 + V_T.C2

			AUX_VAR_T_total = ss_val[QSP_ID_TUM_TREG] + ss_val[QSP_ID_TUM_TEFF] + ss_val[QSP_ID_TUM_TH];  // V_T.T0 + V_T.T1 + V_T.Th

			AUX_VAR_M_total = ss_val[QSP_ID_TUM_M1] + ss_val[QSP_ID_TUM_M2];  // V_T.Mac_M1 + V_T.Mac_M2

			tumor_volume = ((ss_val[20] + AUX_VAR_C_total) / AVOGADROS * QSP_CONST(QSP_VOL_CELL) + 
							(ss_val[21] + ss_val[22] + AUX_VAR_T_total) / AVOGADROS * QSP_CONST(QSP_VOL_TCELL)) / QSP_CONST(QSP_VE_T) + 
							AUX_VAR_M_total / AVOGADROS * QSP_CONST(QSP_VOL_MCELL) / QSP_CONST(QSP_VE_T) + 
							ss_val[48] / AVOGADROS* QSP_CONST(QSP_VOL_FIBCELL) / QSP_CONST(QSP_VE_T) + 
							ss_val[49] / AVOGADROS* QSP_CONST(QSP_VOL_CAFCELL) / QSP_CONST(QSP_VE_T) + 
							ss_val[50] / 1e9 * QSP_CONST(QSP_ECM_MW) / QSP_CONST(QSP_ECM_DENSITY); // ECM in nanomoles, convert to moles
			

			tt += deltatt;
		}
		if (tumor_volume < tumor_volume_ref_presimulation)
		{		
			std::cout<<"QSP: tumor volume condition is not met"<<std::endl;
			exit(0);
		}	
	}

	// setup

	LymphBloodQSP::_QSP_weight = params.getVal(PARAM_WEIGHT_QSP);
	LymphBloodQSP::use_steady_state = false;
	LymphBloodQSP::setup_class_parameters(p);
	_QSP_model.getSystem()->setup_instance_tolerance(p);
	_QSP_model.getSystem()->setup_instance_varaibles(p);

	// load steady state
	if (steadystate)
	{
		for (size_t i = 0; i < n; i++)
		{
			_QSP_model.getSystem()->setSpeciesVar(i, ss_val[i]);
		}
	}

	_QSP_model.getSystem()->eval_init_assignment();
	_QSP_model.getSystem()->updateVar();
}

/*! solve QSP from t to t + dt during presimulation
*/
void LymphCentral::time_step_preSimulation(double t, double dt) {
	
	//No nivolumab is being added (nor any treatment) at this time
	// solve QSP for dt
	_QSP_model.solve(t, dt);

	return;
}

/*! solve QSP from t to t + dt
*/
void LymphCentral::time_step(double t, double dt){

	double final_dose_week = 8;
	double day_per_week = 7;
	if (t / (SEC_PER_DAY * day_per_week) < final_dose_week) {
		// Pharmacokinetics
		if (params.getVal(PARAM_NIVO_ON) != 0)
		{
			double week = t / (SEC_PER_DAY * params.getVal(PARAM_NIVO_DOSE_INTERVAL_TIME));
			int week_int = floor(week);
			double nivo_dose = params.getVal(PARAM_NIVO_DOSE);
			double cent_nivo = GET_VAR(QSP_ID_CENT_NIVO);

			if (week == week_int)
			{
				cent_nivo += nivo_dose;
				SET_VAR(QSP_ID_CENT_NIVO, cent_nivo);
			}
		}

		if (params.getVal(PARAM_IPI_ON) != 0)
		{
			double week = t / (SEC_PER_DAY * params.getVal(PARAM_IPI_DOSE_INTERVAL_TIME));
			int week_int = floor(week);
			double ipi_dose = params.getVal(PARAM_IPI_DOSE);
			//double cent_ipi = GET_VAR(QSP_ID_CENT_IPI);
			double cent_ipi = 0;
			if (week == week_int)
			{
				cent_ipi += ipi_dose;
				//SET_VAR(QSP_ID_CENT_IPI, cent_ipi);
			}
		}

		if (params.getVal(PARAM_CABO_ON) != 0)
		{
			double cabo_interval = SEC_PER_DAY * params.getVal(PARAM_CABO_DOSE_INTERVAL_TIME);
			double f1 = 0.675; //fraction of dose in the absorption site
			double cabo_dose_site1 = f1 * params.getVal(PARAM_CABO_DOSE) / QSP_CONST(QSP_V_C);
			double cabo_dose_site2 = (1 - f1) * params.getVal(PARAM_CABO_DOSE) / QSP_CONST(QSP_V_C);
			double site1_cabo = GET_VAR(QSP_ID_SITE1_CABO);
			double site2_cabo = GET_VAR(QSP_ID_SITE2_CABO);

			if (int(t) % int(cabo_interval) == 0)
			{
				site1_cabo += cabo_dose_site1;
				site2_cabo += cabo_dose_site2;
				SET_VAR(QSP_ID_SITE1_CABO, site1_cabo);
				SET_VAR(QSP_ID_SITE2_CABO, site2_cabo);
			}
		}
	}
	// solve QSP for dt
	_QSP_model.solve(t, dt);

	return;
}

/*! Get QSP variables for ABM.
	
	# Tum.C1 (unit: cell)
    # Cent.Teff (unit: convert from cell to mole)
	# Cent.TCD4 (unit: convert from cell to mole)
	# Tum.MDSC (unit: convert from cell to mole)
	# Tum.Nivo (unit: convert from 1e-6 mole/m^3 to mole/m^3)
	# Tum.ENT (unit: convert from 1e-6 mole/m^3 to mole/m^3)
    # Tum.Cx (unit: convert from cell to mole)
	# Tum.Texh (unit: convert from cell to mole)
	# Tum.CCL2 (unit: convert from 1e-6 mole/m^3 to mole/m^3)	
*/
const std::vector<double>& LymphCentral::get_var_exchange(void){

	// need to be cell count for calculating ABM scalor
	_var_qsp_to_abm[QSPEX_TUM_C] = GET_VAR_RAW(QSP_ID_TUM_C1);
	// internal SI unit for calculating rates and probabilities.
	_var_qsp_to_abm[QSPEX_CENT_TEFF] = GET_VAR(QSP_ID_CENT_TEFF);
	_var_qsp_to_abm[QSPEX_CENT_TREG] = GET_VAR(QSP_ID_CENT_TREG);
	_var_qsp_to_abm[QSPEX_CENT_TH] = GET_VAR(QSP_ID_CENT_TH);
	_var_qsp_to_abm[QSPEX_TUM_TEFF] = GET_VAR(QSP_ID_TUM_TEFF);
	_var_qsp_to_abm[QSPEX_TUM_TREG] = GET_VAR(QSP_ID_TUM_TREG);
	_var_qsp_to_abm[QSPEX_TUM_TH] = GET_VAR(QSP_ID_TUM_TH);
	_var_qsp_to_abm[QSPEX_TUM_NIVO] = GET_VAR(QSP_ID_TUM_NIVO);
	_var_qsp_to_abm[QSPEX_TUM_IPI] = 0;
	_var_qsp_to_abm[QSPEX_TUM_CABO] = GET_VAR(QSP_ID_TUM_CABO);
	_var_qsp_to_abm[QSPEX_TUM_CMAX] = GET_VAR(QSP_C_MAX);
	_var_qsp_to_abm[QSPEX_TUM_CX] = GET_VAR(QSP_ID_TUM_CX);
	_var_qsp_to_abm[QSPEX_TUM_TEXH] = GET_VAR(QSP_ID_TUM_TEXH);	
	_var_qsp_to_abm[QSPEX_TUM_CAF] = GET_VAR(QSP_ID_TUM_CAF);
	_var_qsp_to_abm[QSPEX_TUM_M1] = GET_VAR(QSP_ID_TUM_M1);
	_var_qsp_to_abm[QSPEX_TUM_M2] = GET_VAR(QSP_ID_TUM_M2);
	_var_qsp_to_abm[QSPEX_TUM_MDSC] = GET_VAR(QSP_ID_TUM_MDSC);
	//! tumor volume, calculated in QSP, updated every time step
	double AUX_VAR_C_total = GET_VAR(QSP_ID_TUM_C1) + GET_VAR(QSP_ID_TUM_C2);  // V_T.C1 + V_T.C2

	double AUX_VAR_T_total = GET_VAR(QSP_ID_TUM_TREG) + GET_VAR(QSP_ID_TUM_TH) + GET_VAR(QSP_ID_TUM_TEFF);  // V_T.T0 + V_T.T1 + V_T.Th

	double AUX_VAR_M_total = GET_VAR(QSP_ID_TUM_M1) + GET_VAR(QSP_ID_TUM_M2) + GET_VAR(QSP_ID_TUM_MDSC);  // V_T.Mac_M1 + V_T.Mac_M2 + V_T.Mac_MDSC

	double tumor_volume = ((GET_VAR(QSP_ID_TUM_CX) + AUX_VAR_C_total) * QSP_CONST(QSP_VOL_CELL) + 
					(GET_VAR(QSP_ID_TUM_TEXH) + GET_VAR(QSP_ID_TUM_TH_EXH) + AUX_VAR_T_total) * QSP_CONST(QSP_VOL_TCELL)) / QSP_CONST(QSP_VE_T) + 
					AUX_VAR_M_total * QSP_CONST(QSP_VOL_MCELL) / QSP_CONST(QSP_VE_T) + 
					GET_VAR(QSP_ID_TUM_FIB)  * QSP_CONST(QSP_VOL_FIBCELL) / QSP_CONST(QSP_VE_T) + 
					GET_VAR(QSP_ID_TUM_CAF)  * QSP_CONST(QSP_VOL_CAFCELL) / QSP_CONST(QSP_VE_T) + 
					GET_VAR(QSP_ID_TUM_ECM)  * QSP_CONST(QSP_ECM_MW) / QSP_CONST(QSP_ECM_DENSITY);
	_var_qsp_to_abm[QSPEX_TUM_VOL] = tumor_volume;

	return _var_qsp_to_abm;
}

/*! update QSP module with output from ABM
	unit convert: item to mole
    # cancer cell death (total)
	# cancer cell death (Teff kill)
	# Teff recruitment
	# TCD4 recruitment
*/
void LymphCentral::update_qsp_var(const std::vector<double>& var_abm){

	// convert item to internal units
	double scalar = 1 / AVOGADROS;

	// CC death total, CC death Teff, Teff recruit, TCD4 recruit
	double cc_death_total = var_abm[Tumor::TUMEX_CC_DEATH] * scalar;
	double cc_death_Teff = var_abm[Tumor::TUMEX_CC_T_KILL] * scalar;
	double Teff_recruit = var_abm[Tumor::TUMEX_TEFF_REC] * scalar;
	double TREG_recruit = var_abm[Tumor::TUMEX_TREG_REC] * scalar;
	double TH_recruit = var_abm[Tumor::TUMEX_TH_REC] * scalar;
	std::cout << "DEAD CANCER CELL: " << var_abm[Tumor::TUMEX_CC_DEATH]<<std::endl;

	double p0 = GET_VAR(QSP_ID_P0);
	double p1 = GET_VAR(QSP_ID_P1);
	double factor_p0 = QSP_CONST(QSP_N_T0_CLONES) * QSP_CONST(QSP_P0_C1);
	double factor_p1 = QSP_CONST(QSP_N_T1_CLONES) * QSP_CONST(QSP_P1_C1);
	p0 += cc_death_total * factor_p0; 
	p1 += cc_death_total *  factor_p1; 
	SET_VAR(QSP_ID_P0, p0);
	SET_VAR(QSP_ID_P1, p1);
	double cent_t_eff = GET_VAR(QSP_ID_CENT_TEFF);
	double cent_t_reg = GET_VAR(QSP_ID_CENT_TREG);
	double cent_t_h = GET_VAR(QSP_ID_CENT_TH);

	cent_t_eff -= Teff_recruit;
	cent_t_reg -= TREG_recruit;
	cent_t_h -= TH_recruit;
	cent_t_eff = std::max(0.0, cent_t_eff);
	cent_t_reg = std::max(0.0, cent_t_reg);
	cent_t_h = std::max(0.0, cent_t_h);

	SET_VAR(QSP_ID_CENT_TEFF, cent_t_eff);
	SET_VAR(QSP_ID_CENT_TREG, cent_t_reg);
	SET_VAR(QSP_ID_CENT_TH, cent_t_h);
	return;
}

};
};
