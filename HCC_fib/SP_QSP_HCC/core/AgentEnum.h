#pragma once

#include <utility>

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

//! Agent type enum
/*!
  Each agent class has one unique type.
  Thus no need to use member variable to store.
  Just implement getType() function to access its type.
*/
enum AgentTypeEnum 
{
	AGENT_DUMMY,
	CELL_TYPE_CANCER,
	CELL_TYPE_T,
	CELL_TYPE_TCD4,
	CELL_TYPE_MDSC,
	CELL_TYPE_MAC,
	CELL_TYPE_FIB,
	CELL_TYPE_APC,
	CELL_TYPE_VAS
};

//! Agent state emum
/*! 
  Each Agent class may have multiple states.
  Default state for classes with only one state and will
  not change during a simulation.
*/
enum AgentStateEnum
{
	// default
	DEFAULT_STATE,

	// cancer cell states: PDL1
	CANCER_PDL1_POS,
	CANCER_PDL1_NEG,

	// t cell states
	T_CELL_EFF,
	T_CELL_CYT,
	T_CELL_SUPP,

	// cancer cell stem/progenitor/senescence
	CANCER_STEM,
	CANCER_PROGENITOR,
	CANCER_SENESCENT,

	// cd4 t cell states
	TCD4_Th,
	TCD4_TREG,

	// APC
	APC_NATIVE,
	APC_MATURE,

	//MAC
	MAC_M1,
	MAC_M2,

	//fibroblast
	FIB_NORMAL,
	FIB_CAF,

	//vasculature
	VAS_TIP,
	VAS_STALK,
	VAS_PHALANX

//typedef std::pair<AgentType, AgentState> AgentTypeState;

};

};
};
