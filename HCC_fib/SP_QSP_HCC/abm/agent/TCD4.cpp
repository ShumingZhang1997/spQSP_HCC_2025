//#include <boost/serialization/export.hpp>
#include "TCD4.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(TCD4)

#include <iostream>
#include <sstream>

#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"
#include "TCell.h"

namespace SP_QSP_IO {
namespace SP_QSP_HCC {

	TCD4::TCD4(SpatialCompartment* c)
		:Cell_Tumor(c)
		, _divide_cd_TCD4_exp(params.getVal(PARAM_TCD4_DIV_INTERVAL))
		, _divide_limit_TCD4_exp(params.getVal(PARAM_TCD4_DIV_LIMIT))
		, _divide_flag(false)
		, _IL2_exposure(0)
		, _TGFB_release_remain(params.getVal(PARAM_TGFB_RELEASE_TIME))
		//, _source_IFNg(NULL)
		, _source_IL_2(NULL)
		, _source_IL_10(NULL)
		, _source_TGFB(NULL)
		, _CTLA4(0)
	{
		_state = AgentStateEnum::TCD4_Th;
		_life = getTCD4Life();
	}

	TCD4::TCD4(const TCD4& c)
		:Cell_Tumor(c)
		, _divide_cd_TCD4_exp(params.getVal(PARAM_TCD4_DIV_INTERVAL))
		, _divide_limit_TCD4_exp(c._divide_limit_TCD4_exp)
		, _divide_flag(c._divide_flag)
		, _IL2_exposure(c._IL2_exposure)
		, _TGFB_release_remain(params.getVal(PARAM_TGFB_RELEASE_TIME))
		//, _source_IFNg(NULL)
		, _source_IL_2(NULL)
		, _source_IL_10(NULL)
		, _source_TGFB(NULL)
		, _CTLA4(c._CTLA4)
		
	{
		_state = c._state;
		_life = getTCD4Life();
		if (_state == TCD4_Th) {
			//setup_chem_source(_source_IFNg, CHEM_IFN, 0);
			setup_chem_source(_source_IL_2, CHEM_IL_2, 0);
			setup_chem_source(_source_IL_10, CHEM_IL_10, 0);
			setup_chem_source(_source_TGFB, CHEM_TGFB, 0);
		}
		if (_state == AgentStateEnum::TCD4_TREG)
		{
			//setup_chem_source(_source_IFNg, CHEM_IFN, 0);
			setup_chem_source(_source_IL_2, CHEM_IL_2, 0);
			setup_chem_source(_source_IL_10, CHEM_IL_10, params.getVal(PARAM_TREG_IL_10_RELEASE));
			setup_chem_source(_source_TGFB, CHEM_TGFB, params.getVal(PARAM_TREG_TGFB_RELEASE));
		}
	}

	TCD4::~TCD4()
	{
	}

	std::string TCD4::toString()const {
		std::stringstream ss;
		ss << Cell_Tumor::toString();
		return ss.str();
	}

	/*void TCD4::setDead(void)
	{
		//remove_source_sink(_source_IL_10);
		CellAgent::setDead();
	}*/

	bool TCD4::agent_movement_step(double t, double dt, Coord& c)  {

	bool move = false;
	// if ECM density is high, cell is unlikely to move
	double ECM_density = get_tumor().get_coord_ECM(c);
	double ECM_sat = ECM_density / (ECM_density + params.getVal(PARAM_FIB_ECM_MOT_EC50));
	if (rng.get_unif_01() < ECM_sat){
		return move;
	}

	if (_tumble == false) {
		//calculate the probability of tumble according to tgfb

		double v_x = _moveDirection.x / params.getVal(PARAM_SEC_PER_TIME_SLICE);
		double v_y = _moveDirection.y / params.getVal(PARAM_SEC_PER_TIME_SLICE);
		double v_z = _moveDirection.z / params.getVal(PARAM_SEC_PER_TIME_SLICE);

		
		std::vector< std::vector<double> > gradient = get_tumor().get_chem_grid().get_gradient(getCoord());

		/*! 
			get gradient of(x,y,z); with given dimension d; x, y, z  are integer voxel coordinates
			d=0: x axis
			d=1: y axis
			d=2: z axis
			substrate i
		*/
		//! calculate the dot product and norm of the gradient
		double dot_product = 0;
		double norm_gradient = 0;
		dot_product = v_x * gradient[CHEM_IFN][0] + v_y * gradient[CHEM_IFN][1] + v_z * gradient[CHEM_IFN][2];
		norm_gradient = std::sqrt(gradient[CHEM_IFN][0] * gradient[CHEM_IFN][0] + gradient[CHEM_IFN][1] * gradient[CHEM_IFN][1] + gradient[CHEM_IFN][2] * gradient[CHEM_IFN][2]);
		// T helper follows IL12 gradient
		//if (_state == AgentStateEnum::TCD4_Th) {
		//	dot_product = v_x * gradient[CHEM_IFN][0] + v_y * gradient[CHEM_IFN][1] + v_z * gradient[CHEM_IFN][2];
		//	norm_gradient = std::sqrt(gradient[CHEM_IFN][0] * gradient[CHEM_IFN][0] + gradient[CHEM_IFN][1] * gradient[CHEM_IFN][1] + gradient[CHEM_IFN][2] * gradient[CHEM_IFN][2]);
		//}
		// Treg follows IL10 gradient
		//else if (_state == AgentStateEnum::TCD4_TREG) {
		//	dot_product = v_x * gradient[CHEM_ARGI][0] + v_y * gradient[CHEM_ARGI][1] + v_z * gradient[CHEM_ARGI][2];
		//	norm_gradient = std::sqrt(gradient[CHEM_ARGI][0] * gradient[CHEM_ARGI][0] + gradient[CHEM_ARGI][1] * gradient[CHEM_ARGI][1] + gradient[CHEM_ARGI][2] * gradient[CHEM_ARGI][2]);
		//}
		//double delta = 5e-6;
		double lambda = 0.0000168;
		double norm_v = std::sqrt(v_x * v_x + v_y * v_y + v_z * v_z);

		double cos_theta = dot_product / (norm_v * norm_gradient);

		double EC50_grad = 1;
		//Hill function of cell sensing the chemokine gradient, higher gradient makes cell easier to move (less likely to tumble)
		double H_grad = norm_gradient / (norm_gradient + EC50_grad);
		// if movement direction is opposite with tgfb gradient
		if (cos_theta < 0) {
			H_grad = -H_grad;
		}
		double Tumble_rate = (lambda/2) * (1 - cos_theta) * (1 - H_grad) * params.getVal(PARAM_SEC_PER_TIME_SLICE);
		//double Tumble_rate = 0.3;
		// poisson process (tumble times N >= 1)
		double pTumble = 1 - std::exp(-Tumble_rate);
		/*
		if (_state == AgentStateEnum::TCD4_Th) {
			std::cout << "IFN concentration: " << get_tumor().get_chem(getCoord(), CHEM_IFN) << ", x grad: " << gradient[CHEM_IFN][0] << ", y grad: " << gradient[CHEM_IFN][1] << ", z grad: " << gradient[CHEM_IFN][2] << std::endl;
			std::cout << "v_x: " << v_x << ", v_y: " << v_y << ", v_z: " << v_z << ", dot_product: " << dot_product 
			<< ", norm_v: " << norm_v << ", norm_gradient: " << norm_gradient << ", cos_theta: " << cos_theta 
			<< ", EC50_grad: " << EC50_grad << ", H_grad: " << H_grad << ", Tumble_rate: " << Tumble_rate << ", pTumble: " << pTumble << std::endl;
		}
		if (_state == AgentStateEnum::TCD4_TREG) {
			std::cout << "ARGI concentration: " << get_tumor().get_chem(getCoord(), CHEM_ARGI) << ", x grad: " << gradient[CHEM_ARGI][0] << ", y grad: " << gradient[CHEM_ARGI][1] << ", z grad: " << gradient[CHEM_ARGI][2] << std::endl;
			std::cout << "v_x: " << v_x << ", v_y: " << v_y << ", v_z: " << v_z << ", dot_product: " << dot_product 
			<< ", norm_v: " << norm_v << ", norm_gradient: " << norm_gradient << ", cos_theta: " << cos_theta 
			<< ", EC50_grad: " << EC50_grad << ", H_grad: " << H_grad << ", Tumble_rate: " << Tumble_rate << ", pTumble: " << pTumble << std::endl;
		}
		*/
		if (rng.get_unif_01() < pTumble)
		{
			_tumble = true;
			//std::cout << "cell tumbled, no move" << std::endl;
			return move;
		}
	}
	
	// when the cell decides if it moves it tracks the tumble status from last step
	if (_tumble == true)
	{
		// move
		int idx;
		const auto shape = getCellShape();
		if (_compartment->getOneOpenVoxel(shape->getMoveDestinationVoxels(),
			shape->getMoveDirectionAnchor(), getCoord(), getType(), idx, rng))
		{
			move = true;
			c = getCellShape()->getMoveDirectionAnchor()[idx] + getCoord();
			_moveDirection = getCellShape()->getMoveDirectionAnchor()[idx];
			//std::cout << "Tumble: Previous location: " << getCoord() << ", move direction: " << _moveDirection << ", new location: " << c << std::endl;
			_tumble = false;
		}
	}
	else {
		Coord3D target_voxel = _moveDirection + getCoord();
		if (_compartment->voxelIsOpen(target_voxel, getType())) {
			move = true;
			c = target_voxel;
			//std::cout << "Run: Previous location: " << getCoord() << ", move direction: " << _moveDirection << ", new location: " << target_voxel << std::endl;
		}
	}
	
	return move;
}

	bool TCD4::agent_state_step(double t, double dt, Coord& c) {
		bool divide = false;
		if (!isDead())
		{
			_life--;
			if (_life <= 0)
			{
				setDead();
				// remove source when cell die
				return divide;
			}
		}

		const auto shape = getCellShape();
		Cell_Tumor::agent_state_step(t, dt, c);

		auto tumor = dynamic_cast<Tumor*>(_compartment);
		double tumvol = tumor->get_Tum_Vol();
		double Treg_count = tumor->get_Treg();
		double Th_count = tumor->get_Th();
		double tum_ipi = tumor->get_Ipi();
		// ARGI and TGFB induced Treg transformation
		double ArgI = get_tumor().get_chem(getCoord(), CHEM_ARGI);
		double TGFB = get_tumor().get_chem(getCoord(), CHEM_TGFB);

		bool cancer_found = _compartment->hasTypeStateInTarget(shape->getEnvironmentLocations(), 
			getCoord(), AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_PROGENITOR);

		if (!isDead() && _state == AgentStateEnum::TCD4_Th)
		{
			//std::cout << "TCD4 Th divide interval: " << _divide_cd_TCD4_exp << std::endl;
			if (cancer_found)
			{
				//std::cout << "TCD4 Th found cancer" << std::endl;
				//update_chem_source(_source_IFNg, params.getVal(PARAM_IFN_G_RELEASE));
				update_chem_source(_source_IL_2, params.getVal(PARAM_IL_2_RELEASE));
				double IL2 = get_tumor().get_chem(getCoord(), CHEM_IL_2);
				//_divide_flag = true;
			}
			else {
				//update_chem_source(_source_IFNg, 0.0);
				update_chem_source(_source_IL_2, 0.0);
				//_divide_flag = false;
			}

			_divide_cd_TCD4_exp--;
			// same parameter value of TGFB_50 for macrophage polarization and Treg transformation
			double alpha = params.getVal(PARAM_K_TH_TREG) * (1 + TGFB / (TGFB + params.getVal(PARAM_MAC_TGFB_EC50)));

			//! calculate probability of Th to Treg transition based on alpha
			double p_th_treg = 1 - std::exp(-alpha);

			
			// std::cout << "H_TGFB: " << (TGFB / (TGFB + params.getVal(PARAM_MAC_TGFB_EC50))) << ", T helper to Treg constant alpha: " << alpha << ", probability p_th_treg: " << p_th_treg << std::endl;
			
			if (rng.get_unif_01() < p_th_treg)
			{
				setTreg();
				// std::cout << "Treg transformation by TGFB polarization" << std::endl;
				return divide;
			}
		}

		
		if (!isDead() && _state == AgentStateEnum::TCD4_TREG)
		{
			// TGFB release, if cancer is found and TGFB release time is not over
			bool cancer_found = _compartment->hasTypeStateInTarget(shape->getEnvironmentLocations(), 
			getCoord(), AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_PROGENITOR);

			if (_TGFB_release_remain >= 0 && cancer_found)
			{
				update_chem_source(_source_TGFB, params.getVal(PARAM_TREG_TGFB_RELEASE));
				// TGFB release time limit
				_TGFB_release_remain -= params.getVal(PARAM_SEC_PER_TIME_SLICE);
			}
			else 
			{
				update_chem_source(_source_TGFB, 0.0);
			}

			double p_ADCC_death = get_CTLA4_Ipi(tum_ipi);
			if (!isDead())
			{
				if (rng.get_unif_01() < p_ADCC_death)
				{
					setDead();
					// remove source when cell die
					// std::cout << "TCD4 dead due to ADCC" << std::endl;
					return divide;
				}
			}

			//ArgI induced Treg expansion
			double H_ArgI = ArgI / (ArgI + params.getVal(PARAM_MDSC_EC50_ArgI_Treg));
			//std::cout << "H_ArgI: " << ArgI / (ArgI + params.getVal(PARAM_MDSC_EC50_ArgI_Treg)) << ", ArgI: " << ArgI <<  std::endl;

			//ArgI induced Treg 
			if (rng.get_unif_01() < H_ArgI && _divide_limit_TCD4_exp > 0)
			{
				_divide_cd_TCD4_exp = 0;
			}

			_divide_cd_TCD4_exp--;
		}

		if (_divide_limit_TCD4_exp > 0 && _divide_cd_TCD4_exp <= 0)
		{
			int idx = 0;
			if (_compartment->getOneOpenVoxel(shape->getProlifDestinationVoxels(),
				shape->getProlifDestinationAnchor(), getCoord(), getType(), idx, rng))
			{
				divide = true;
				//std::cout << "TCD4 divide" << std::endl;
				//std::cout << "idx: " << idx << ", " << getCellShape()->getProlif()[idx] << std::endl;
				c = getCellShape()->getProlifDestinationAnchor()[idx] + getCoord();

				_divide_limit_TCD4_exp -= 1;
				//_divide_cd_TCD4_exp = int(params.getVal(PARAM_TCD4_EXP_INTERVAL_SLICE) / (ArgI / (ArgI + params.getVal(PARAM_EC50_ARGI_TCD4)) * (1 - ent / (ent + params.getVal(PARAM_IC50_ENT_ARGI)))) + .5);
				_divide_cd_TCD4_exp = params.getVal(PARAM_TCD4_DIV_INTERVAL);


			}
		}
			
		return divide;
	}

	void TCD4::move_all_source_sink(void) const
	{
		//move_source_sink(_source_IFNg);
		move_source_sink(_source_IL_2);
		move_source_sink(_source_TGFB);
		move_source_sink(_source_IL_10);
	}

	void TCD4::remove_all_source_sink(void) {
		//remove_source_sink(_source_IFNg);
		remove_source_sink(_source_IL_2);
		remove_source_sink(_source_IL_10);
		remove_source_sink(_source_TGFB);
		return;
	}

	void TCD4::setTh()
	{
		_state = AgentStateEnum::TCD4_Th;
		_CTLA4 = 0;
		//reset the division time cool down
		_divide_cd_TCD4_exp = params.getVal(PARAM_TCD4_LIFE_MEAN);
		//update_chem_source(_source_IFNg, params.getVal(PARAM_IFN_G_RELEASE));
		update_chem_source(_source_IL_2, params.getVal(PARAM_IL_2_RELEASE));
		update_chem_source(_source_TGFB, 0);
		update_chem_source(_source_IL_10, 0);
		return;
	}

	void TCD4::setTreg()
	{
		//std::cout << "Transforming Th to Treg" << std::endl;
		_state = AgentStateEnum::TCD4_TREG;
		_CTLA4 = params.getVal(PARAM_CTLA4_TREG);
		_divide_cd_TCD4_exp = params.getVal(PARAM_TCD4_DIV_INTERVAL);
		//std::cout << "ARGI conc: "<< ArgI << ", cd4 treg divide limit: " << params.getVal(PARAM_TREG_EXP_INTERVAL_SLICE) / (ArgI / (ArgI + params.getVal(PARAM_MDSC_EC50_ArgI_Treg))) + .5 << std::endl;
		//update_chem_source(_source_IFNg, 0);
		update_chem_source(_source_IL_2, 0);
		update_chem_source(_source_IL_10, params.getVal(PARAM_TREG_IL_10_RELEASE));
		return;
	}

	
	int TCD4::getTCD4Life() {

		double lifeMean = params.getVal(PARAM_TCD4_LIFE_MEAN);
		double lifeSd = params.getVal(PARAM_T_CELL_LIFE_SD_SLICE);
		double tLifeD = lifeMean + rng.get_norm_std() * lifeSd;

		int tLife = int(tLifeD + 0.5);
		tLife = tLife > 0 ? tLife : 0;
		//std::cout << "random TCD4 life: " << tLife << std::endl;
		return tLife;
	}

	double TCD4::get_CTLA4_Ipi(double Ipi) {
		//determine the number of bonds between CTLA4 and Ipi using following equations at equilibrium
		// 
		//kon_CTLA4_ipi* (V_T.TCD4_CTLA4 * V_T.ipi / gamma_C_ipi) - koff_CTLA4_ipi * V_T.TCD4_CTLA4_ipi
		//Chi_CTLA4_ipi * kon_CTLA4_ipi * (V_T.TCD4_CTLA4 * V_T.TCD4_CTLA4_ipi) / A_Tcell - koff_CTLA4_ipi * V_T.TCD4_CTLA4_ipi_CTLA4
		//V_T.TCD4_CTLA4 + V_T.TCD4_CTLA4_ipi + V_T.TCD4_CTLA4_ipi_CTLA4 = N

		double k_on = params.getVal(PARAM_KON_CTLA4_IPI);
		double k_off = params.getVal(PARAM_KOFF_CTLA4_IPI);
		double gamma_T_ipi = params.getVal(PARAM_GAMMA_T_IPI);
		double chi_CTLA4 = params.getVal(PARAM_CHI_CTLA4_IPI);
		double a_Tcell = params.getVal(PARAM_A_TCELL);
		double n_CTLA4_TCD4 = params.getVal(PARAM_CTLA4_TREG);

		double a = k_on * chi_CTLA4 / a_Tcell / k_off;
		double b = k_on * Ipi / k_off / gamma_T_ipi;
		//total number of CTLA4 molecule on TCD4 
		double c = n_CTLA4_TCD4;
		double d = -1;

		//std::cout << "a: " << a << ", b: " << b << ", c: " << c << ", d: " << d << std::endl;
		//Newton_Raphson_root
		int max_iter = 20;
		double tol_rel = 1E-5;
		double root = 0;
		double res, root_new, f, f1;
		int i = 0;
		while (i < max_iter) {


			f = 2 * a * b * std::pow(root, 2) + (b + 1) * root - c;
			f1 = 4 * a * b * root + (b + 1);

			root_new = root - f / f1;
			res = std::abs(root_new - root) / root_new;
			if (res > tol_rel) {
				i++;
				root = root_new;
			}
			else {
				break;
			}
		}

		double free_CTLA4 = root;
		double CTLA4_Ipi = b * free_CTLA4;
		double CTLA4_Ipi_CTLA4 = a * b * free_CTLA4 * free_CTLA4;

		/*
		std::cout << "fraction of free CTLA4: " << root
			<< ", free CTLA4: " << free_CTLA4
			<< ", CTLA4_Ipi: " << CTLA4_Ipi
			<< ", CTLA4_Ipi_CTLA4: " << CTLA4_Ipi_CTLA4 <<std::endl;
		*/

		// H_TCD4_T = ((V_T.TCD4_CTLA4_ipi+2*V_T.TCD4_CTLA4_ipi_CTLA4)/TCD4_CTLA4_50)^n_TCD4_CTLA4/(((V_T.TCD4_CTLA4_ipi+2*V_T.TCD4_CTLA4_ipi_CTLA4)/TCD4_CTLA4_50)^n_TCD4_CTLA4 + 1)
		double CTLA_50 = params.getVal(PARAM_TREG_CTLA4_50);
		double H_TCD4 = std::pow(((CTLA4_Ipi + 2 * CTLA4_Ipi_CTLA4) / CTLA_50), 2) / (std::pow(((CTLA4_Ipi + 2 * CTLA4_Ipi_CTLA4) / CTLA_50), 2) + 1);
		//std::cout << "H_TCD4: " << H_TCD4 << std::endl;

		double rate = H_TCD4 * params.getVal(PARAM_K_ADCC);
		double p = 1 - exp(-rate);
		//std::cout << "p_death: " << p << std::endl;
		return p;

	}


};
};