#include "ODE_system.h"
    
#define SPVAR(x) NV_DATA_S(y)[x]
#define NSPVAR(x) ptrOde->_nonspecies_var[x]
#define PARAM(x) _class_parameter[x]
#define PFILE(x) param.getVal(x)

namespace CancerVCT{
#define QSP_W ODE_system::_QSP_weight

bool ODE_system::use_steady_state = false;
bool ODE_system::use_resection = false;
double ODE_system::_QSP_weight = 1.0;

ODE_system::ODE_system()
:CVODEBase()
{
    setupVariables();
    setupEvents();
    setupCVODE();
    update_y_other();
}

ODE_system::ODE_system(const ODE_system& c)
{
    setupCVODE();
}

ODE_system::~ODE_system()
{
}

void ODE_system::initSolver(realtype t){

    restore_y();
    int flag;

    flag = CVodeInit(_cvode_mem, f, t, _y);
    check_flag(&flag, "CVodeInit", 1);

    /* Call CVodeRootInit to specify the root function g */
    flag = CVodeRootInit(_cvode_mem, _nroot, g);
    check_flag(&flag, "CVodeRootInit", 1);
    
    	/*Do not do this. Event only trigger when turn from false to true.
	  If this is reset before trigger evaluation at the beginning of simulation,
	  t=0 events might be missed.*/
    //updateTriggerComponentConditionsOnValue(t);
    //resetEventTriggers();

    return;
} 

state_type ODE_system::_class_parameter = state_type(277, 0);

void ODE_system::setup_class_parameters(Param& param){
    //V_C, mw447d336b_d627_4aca_a678_c2ecc527d6c5, index: 0
    //Unit: metre^(3)
    _class_parameter[0] = PFILE(5) * 0.0010000000000000007;
    //V_P, mw4cdba6b2_c1a7_4911_b534_53467a7678b7, index: 1
    //Unit: metre^(3)
    _class_parameter[1] = PFILE(6) * 0.0010000000000000007;
    //V_LN, mw8e8d62d4_92de_423b_86c8_3735d614371a, index: 2
    //Unit: metre^(3)
    _class_parameter[2] = PFILE(8) * 1.0000000000000013e-09;
    //V_e, mw5159e97a_33f5_4a97_845a_3234fab4a1aa, index: 3
    //Unit: metre^(3)
    _class_parameter[3] = PFILE(9) * 0.0010000000000000002;
    //A_e, mw6a5d2c2d_dfeb_4b1b_ad8c_4cdc07fd2574, index: 4
    //Unit: metre^(2)
    _class_parameter[4] = PFILE(10) * 1e-12;
    //A_s, mw28485be3_4442_4825_a9b6_9b0851e97e33, index: 5
    //Unit: metre^(2)
    _class_parameter[5] = PFILE(11) * 1e-12;
    //syn_T_C1, mw49251f78_1263_4781_99d8_3c68ed5cd8f8, index: 6
    //Unit: metre^(2)
    _class_parameter[6] = PFILE(12) * 1e-12;
    //syn_T_APC, mwd43b5439_6fe4_46fc_ad28_645d36ac6de6, index: 7
    //Unit: metre^(2)
    _class_parameter[7] = PFILE(13) * 1e-12;
    //syn_M_C, mwe23344c0_0523_42d0_98ff_5710bee3e60d, index: 8
    //Unit: metre^(2)
    _class_parameter[8] = PFILE(14) * 1e-12;
    //k_cell_clear, mwb5df1c3f_1cb2_486a_a2d6_94674ea00fa8, index: 9
    //Unit: second^(-1)
    _class_parameter[9] = PFILE(170) * 1.15740740740741e-05;
    //cell, mw3babb9cd_c324_485e_a4df_c931523274e1, index: 10
    //Unit: mole^(1)
    _class_parameter[10] = PFILE(171) * 1.66053872801495e-24;
    //day, mwde11fb8e_e77b_4377_919b_6574159fd6e8, index: 11
    //Unit: second^(1)
    _class_parameter[11] = PFILE(172) * 86400.0;
    //vol_cell, mw67cf8754_a07e_4861_a28e_4126bb7179ed, index: 12
    //Unit: metre^(3)mole^(-1)
    _class_parameter[12] = PFILE(173) * 602214.1989999996;
    //vol_Tcell, mw084a7f49_ed7d_4fc0_987b_0a72c5cdd9a3, index: 13
    //Unit: metre^(3)mole^(-1)
    _class_parameter[13] = PFILE(174) * 602214.1989999996;
    //Ve_T, mwfae255b5_cd84_4fd8_b417_c3274e1724d0, index: 14
    //Unit: dimensionless^(1)
    _class_parameter[14] = PFILE(175) * 1.0;
    //k_C1_growth, mw3841a431_2b43_46cc_a876_f49c41238fea, index: 15
    //Unit: second^(-1)
    _class_parameter[15] = PFILE(188) * 1.15740740740741e-05;
    //k_C1_death, mwcf3b6450_696c_4bc2_8198_05b113491c33, index: 16
    //Unit: second^(-1)
    _class_parameter[16] = PFILE(190) * 1.15740740740741e-05;
    //initial_tumour_diameter, mw5601ebb1_16fb_4a1b_bbb0_40fbfc5044cc, index: 17
    //Unit: metre^(1)
    _class_parameter[17] = PFILE(192) * 0.01;
    //k_K_g, mw22eda974_659e_47ff_8aa2_e5a3c938b36f, index: 18
    //Unit: second^(-1)
    _class_parameter[18] = PFILE(193) * 1.15740740740741e-05;
    //k_K_d, mw2653c4b8_c546_45fe_b396_9a2a7da9e775, index: 19
    //Unit: second^(-1)
    _class_parameter[19] = PFILE(194) * 1.15740740740741e-05;
    //k_vas_Csec, mwbff9cc9e_aa67_4ebb_9cf6_7819ce781ebd, index: 20
    //Unit: kilogram^(1)mole^(-1)second^(-1)
    _class_parameter[20] = PFILE(195) * 6970.0717476851905;
    //k_vas_deg, mw76d5b688_a703_4655_98c4_1ce6beb26354, index: 21
    //Unit: second^(-1)
    _class_parameter[21] = PFILE(196) * 1.15740740740741e-05;
    //c_vas_50, mw5bb2785d_dba5_4f53_b6d8_8d84cba5afc0, index: 22
    //Unit: kilogram^(1)metre^(-3)
    _class_parameter[22] = PFILE(197) * 1e-09;
    //k_C2_growth, mw901ac016_5b1b_4bfb_9d78_8517a96f7eca, index: 23
    //Unit: second^(-1)
    _class_parameter[23] = PFILE(198) * 1.15740740740741e-05;
    //k_C2_death, mw254f3149_a13e_49f2_883c_5bb5c4d350de, index: 24
    //Unit: second^(-1)
    _class_parameter[24] = PFILE(199) * 1.15740740740741e-05;
    //div_T0, mwea8799ed_78b9_4201_b612_671faf1b43a6, index: 25
    //Unit: dimensionless^(1)
    _class_parameter[25] = PFILE(201) * 1.0;
    //n_T0_clones, mw8e422a79_821e_43b1_bf65_4bbdaa8c1db8, index: 26
    //Unit: dimensionless^(1)
    _class_parameter[26] = PFILE(202) * 1.0;
    //q_nT0_LN_in, mwe67b44ff_72b2_4607_966d_c6789085a72c, index: 27
    //Unit: second^(-1)
    _class_parameter[27] = PFILE(203) * 1.15740740740741e-05;
    //q_T0_LN_out, mwf19ee9d6_d288_4ecd_8260_6b68740bbb50, index: 28
    //Unit: second^(-1)
    _class_parameter[28] = PFILE(204) * 1.15740740740741e-05;
    //q_nT0_LN_out, mwa00e6f2d_66eb_4101_b27a_1e888a5fbc36, index: 29
    //Unit: second^(-1)
    _class_parameter[29] = PFILE(205) * 1.15740740740741e-05;
    //k_T0_act, mw5086ae8b_8756_43aa_ad0c_95c78cbc56ed, index: 30
    //Unit: second^(-1)
    _class_parameter[30] = PFILE(206) * 1.15740740740741e-05;
    //k_T0_pro, mwc64f8d14_d998_42cc_bc05_110bb7ce655c, index: 31
    //Unit: second^(-1)
    _class_parameter[31] = PFILE(207) * 1.15740740740741e-05;
    //k_T0_death, mw1b702bf9_2357_482b_891b_b562ce07dca8, index: 32
    //Unit: second^(-1)
    _class_parameter[32] = PFILE(208) * 1.15740740740741e-05;
    //q_T0_P_in, mw7dfd9b2d_ba04_47c1_8ec4_159c8010cd2f, index: 33
    //Unit: second^(-1)
    _class_parameter[33] = PFILE(209) * 0.0166666666666667;
    //q_T0_P_out, mwd9b56c91_14f4_4b58_8350_08df79d96d2d, index: 34
    //Unit: second^(-1)
    _class_parameter[34] = PFILE(210) * 1.15740740740741e-05;
    //q_T0_T_in, mwc7e52afa_8e27_43c4_ae5c_7229a4a03d64, index: 35
    //Unit: metre^(-3)second^(-1)
    _class_parameter[35] = PFILE(211) * 16666.666666666693;
    //q_nT0_P_in, mwac162bb7_568e_48f4_8e64_94f123e45fe3, index: 36
    //Unit: second^(-1)
    _class_parameter[36] = PFILE(212) * 0.0166666666666667;
    //q_nT0_P_out, mw659be865_72cc_49e8_981a_c40cf68340b1, index: 37
    //Unit: second^(-1)
    _class_parameter[37] = PFILE(213) * 1.15740740740741e-05;
    //Q_nT0_thym, mw182210e3_de8e_44e6_8d97_69ff0b831095, index: 38
    //Unit: mole^(1)second^(-1)
    _class_parameter[38] = PFILE(214) * 1.92191982409137e-29;
    //k_nT0_pro, mw74d6a9f8_d8e6_4f09_ab0e_affb3ba1099c, index: 39
    //Unit: mole^(1)second^(-1)
    _class_parameter[39] = PFILE(215) * 1.92191982409137e-29;
    //K_nT0_pro, mw6ca4b231_768f_4fd7_a120_8b8eddd1ea73, index: 40
    //Unit: mole^(1)
    _class_parameter[40] = PFILE(216) * 1.66053872801495e-24;
    //k_nT0_death, mw0ec7b3b6_6be1_499d_af54_a65b0531d6db, index: 41
    //Unit: second^(-1)
    _class_parameter[41] = PFILE(217) * 1.15740740740741e-05;
    //k_IL2_deg, mwc17396b0_c8c2_43c6_aaaa_d1eff3b88b22, index: 42
    //Unit: second^(-1)
    _class_parameter[42] = PFILE(218) * 0.0166666666666667;
    //k_IL2_cons, mw7186370a_109c_4a6f_9c63_7b762041f59a, index: 43
    //Unit: second^(-1)
    _class_parameter[43] = PFILE(219) * 167281721944.444;
    //k_IL2_sec, mwc7b0f143_de1b_4228_bffe_bd4e09ef72a5, index: 44
    //Unit: second^(-1)
    _class_parameter[44] = PFILE(220) * 167281721944.444;
    //IL2_50, mwabc21c57_532a_4331_90d9_1c8e82972b32, index: 45
    //Unit: metre^(-3)mole^(1)
    _class_parameter[45] = PFILE(221) * 1.0000000000000008e-06;
    //IL2_50_Treg, mw9e7c3585_b24a_4507_b708_d18c23c7782d, index: 46
    //Unit: metre^(-3)mole^(1)
    _class_parameter[46] = PFILE(222) * 1.0000000000000008e-06;
    //N0, mwa69714be_76fe_4667_8af9_371f27d0c85b, index: 47
    //Unit: dimensionless^(1)
    _class_parameter[47] = PFILE(223) * 1.0;
    //N_costim, mwb86d17a1_97b2_4b5e_8616_defe966da1ca, index: 48
    //Unit: dimensionless^(1)
    _class_parameter[48] = PFILE(224) * 1.0;
    //N_IL2_CD8, mwdfd09d93_2eb2_400d_8c5a_7bd2ace86ca0, index: 49
    //Unit: dimensionless^(1)
    _class_parameter[49] = PFILE(225) * 1.0;
    //N_IL2_CD4, mwa823a51d_a79a_4e75_a879_83ed31e32d01, index: 50
    //Unit: dimensionless^(1)
    _class_parameter[50] = PFILE(226) * 1.0;
    //k_Treg, mw4de985be_fbe9_451c_a9d6_0b6e852bbb47, index: 51
    //Unit: second^(-1)
    _class_parameter[51] = PFILE(227) * 1.15740740740741e-05;
    //div_T1, mw5f5c10e2_90de_4395_89e2_d6d38537c3db, index: 52
    //Unit: dimensionless^(1)
    _class_parameter[52] = PFILE(231) * 1.0;
    //n_T1_clones, mw0bff1abc_c490_42de_b1d6_cb7bddefe464, index: 53
    //Unit: dimensionless^(1)
    _class_parameter[53] = PFILE(232) * 1.0;
    //q_nT1_LN_in, mwbcfbfcd8_8d76_4525_8682_82c068323dff, index: 54
    //Unit: second^(-1)
    _class_parameter[54] = PFILE(233) * 1.15740740740741e-05;
    //q_T1_LN_out, mw60a12da8_a62a_417a_9a55_2c178a89c415, index: 55
    //Unit: second^(-1)
    _class_parameter[55] = PFILE(234) * 1.15740740740741e-05;
    //q_nT1_LN_out, mwb4f2070c_3a32_438e_ad77_cdea203c08f0, index: 56
    //Unit: second^(-1)
    _class_parameter[56] = PFILE(235) * 1.15740740740741e-05;
    //k_T1_act, mwc79c750f_8601_49e2_98bf_70ef7d623929, index: 57
    //Unit: second^(-1)
    _class_parameter[57] = PFILE(236) * 1.15740740740741e-05;
    //k_T1_pro, mwc36e507e_ae63_4c44_989e_c2926bc82841, index: 58
    //Unit: second^(-1)
    _class_parameter[58] = PFILE(237) * 1.15740740740741e-05;
    //k_T1_death, mwe6596f33_c5dc_4d80_ae06_991ad6cc926c, index: 59
    //Unit: second^(-1)
    _class_parameter[59] = PFILE(238) * 1.15740740740741e-05;
    //q_T1_P_in, mw2cb333cc_8491_42ae_a1d6_d4daf5c14238, index: 60
    //Unit: second^(-1)
    _class_parameter[60] = PFILE(239) * 0.0166666666666667;
    //q_T1_P_out, mw177f669e_95a3_4982_806b_465aafd4945a, index: 61
    //Unit: second^(-1)
    _class_parameter[61] = PFILE(240) * 1.15740740740741e-05;
    //q_T1_T_in, mw9c4acc69_b282_44d5_803c_fc17066938fe, index: 62
    //Unit: metre^(-3)second^(-1)
    _class_parameter[62] = PFILE(241) * 16666.666666666693;
    //q_nT1_P_in, mwfa8f5b8f_7e76_42e7_8b26_8ac45a8375ad, index: 63
    //Unit: second^(-1)
    _class_parameter[63] = PFILE(242) * 0.0166666666666667;
    //q_nT1_P_out, mwa3d13cf1_f8f5_4b77_a34c_5d2fcf77793c, index: 64
    //Unit: second^(-1)
    _class_parameter[64] = PFILE(243) * 1.15740740740741e-05;
    //Q_nT1_thym, mwa02288b9_1aaa_4999_8a16_8791dcf99865, index: 65
    //Unit: mole^(1)second^(-1)
    _class_parameter[65] = PFILE(244) * 1.92191982409137e-29;
    //k_nT1_pro, mw2807d248_3dde_49b0_933f_0cd7020e3d1c, index: 66
    //Unit: mole^(1)second^(-1)
    _class_parameter[66] = PFILE(245) * 1.92191982409137e-29;
    //K_nT1_pro, mw5ad2fe4f_5be0_4b18_90ca_58824f4230be, index: 67
    //Unit: mole^(1)
    _class_parameter[67] = PFILE(246) * 1.66053872801495e-24;
    //k_nT1_death, mwde78f44b_f5a3_41b8_aba6_27fb6a3ca871, index: 68
    //Unit: second^(-1)
    _class_parameter[68] = PFILE(247) * 1.15740740740741e-05;
    //k_T1, mwdf82bfd6_aa7e_4c22_8940_8ec5f85fd018, index: 69
    //Unit: second^(-1)
    _class_parameter[69] = PFILE(248) * 1.15740740740741e-05;
    //k_C_T1, mw22f858f1_d23f_4890_a866_204c774877ef, index: 70
    //Unit: second^(-1)
    _class_parameter[70] = PFILE(249) * 1.15740740740741e-05;
    //k_Tcell_ECM, mwf40ea669_c057_48d1_ad8f_28aa2201857a, index: 71
    //Unit: second^(-1)
    _class_parameter[71] = PFILE(250) * 1.15740740740741e-05;
    //k_IFNg_Tsec, mw466ce451_0d7f_42c6_84f7_7607dfb099c3, index: 72
    //Unit: second^(-1)
    _class_parameter[72] = PFILE(253) * 6970071747.68519;
    //K_T_C, mw649eed5d_37cc_470a_a2d8_ab3eeb0d8427, index: 73
    //Unit: dimensionless^(1)
    _class_parameter[73] = PFILE(255) * 1.0;
    //K_T_Treg, mwc3246427_b5b2_4429_a769_243936385711, index: 74
    //Unit: dimensionless^(1)
    _class_parameter[74] = PFILE(256) * 1.0;
    //k_APC_mat, mw69d0baa1_624d_41ce_8ae9_4f1e2f68bebe, index: 75
    //Unit: second^(-1)
    _class_parameter[75] = PFILE(257) * 1.15740740740741e-05;
    //k_APC_mig, mwb8514ae8_739f_443c_9c47_4d327fb47e32, index: 76
    //Unit: second^(-1)
    _class_parameter[76] = PFILE(258) * 1.15740740740741e-05;
    //k_APC_death, mw5f435670_8cd8_468c_970b_269b499c65c1, index: 77
    //Unit: second^(-1)
    _class_parameter[77] = PFILE(259) * 1.15740740740741e-05;
    //k_mAPC_death, mw0d7c7656_34c7_45a5_a323_4618b69267c3, index: 78
    //Unit: second^(-1)
    _class_parameter[78] = PFILE(260) * 1.15740740740741e-05;
    //APC0_T, mw4c56d425_9419_43d5_9b1f_6d9eb37906ba, index: 79
    //Unit: metre^(-3)mole^(1)
    _class_parameter[79] = PFILE(261) * 1.6605387280149534e-18;
    //APC0_LN, mw9ede0315_0ae0_4c1d_ad4d_53febbc661f3, index: 80
    //Unit: metre^(-3)mole^(1)
    _class_parameter[80] = PFILE(262) * 1.6605387280149534e-18;
    //n_sites_APC, mw50afb61e_0e53_4fcc_9b60_1118ad985397, index: 81
    //Unit: dimensionless^(1)
    _class_parameter[81] = PFILE(263) * 1.0;
    //kin, mw0b05ecc6_6c68_4b5f_a3e4_639e7a0c7bdb, index: 82
    //Unit: second^(-1)
    _class_parameter[82] = PFILE(264) * 1.15740740740741e-05;
    //kout, mwaa59a431_d069_475c_82d9_ef959cce3df6, index: 83
    //Unit: second^(-1)
    _class_parameter[83] = PFILE(265) * 1.15740740740741e-05;
    //k_P0_up, mwcbcee841_070f_4b93_b551_3ef11dec034b, index: 84
    //Unit: mole^(-1)second^(-1)
    _class_parameter[84] = PFILE(266) * 6.97007174768519e+18;
    //k_xP0_deg, mw216250d6_b9b1_4422_aacb_8a869a52a8cb, index: 85
    //Unit: second^(-1)
    _class_parameter[85] = PFILE(267) * 1.15740740740741e-05;
    //k_P0_deg, mwe40aa6e0_346e_4b7d_8c92_f49ce2d8221e, index: 86
    //Unit: second^(-1)
    _class_parameter[86] = PFILE(268) * 1.15740740740741e-05;
    //k_p0_deg, mwa7eb6936_2824_4b2a_ac26_d0375b29e409, index: 87
    //Unit: second^(-1)
    _class_parameter[87] = PFILE(269) * 1.15740740740741e-05;
    //k_P0_on, mwfbb81ef5_b4fa_4f23_884e_29aa06b53610, index: 88
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[88] = PFILE(270) * 1.1574074074074112e-08;
    //k_P0_d1, mw67cae976_ef30_4c17_92d0_a640c4bd72b5, index: 89
    //Unit: metre^(-3)mole^(1)
    _class_parameter[89] = PFILE(271) * 999.9999999999994;
    //p0_50, mw40aef2ed_2d1f_4126_84a3_fb3bf301ae24, index: 90
    //Unit: metre^(-2)mole^(1)
    _class_parameter[90] = PFILE(272) * 1.66053872801495e-12;
    //P0_C1, mw53713495_b6fd_400e_b5be_7f024958a44a, index: 91
    //Unit: dimensionless^(1)
    _class_parameter[91] = PFILE(273) * 6.02214199e+23;
    //P0_C2, mw93d016c7_1c17_40bb_b9d3_576f78dd9a50, index: 92
    //Unit: dimensionless^(1)
    _class_parameter[92] = PFILE(274) * 6.02214199e+23;
    //A_syn, mw3458b4bc_36d7_486d_a973_f2322a4806a5, index: 93
    //Unit: metre^(2)
    _class_parameter[93] = PFILE(275) * 1e-12;
    //A_Tcell, mwa53bd64b_15fb_4ed2_926d_69cf8fe42320, index: 94
    //Unit: metre^(2)
    _class_parameter[94] = PFILE(276) * 1e-12;
    //A_cell, mw52f710d2_a8a0_42a5_b327_99283f547473, index: 95
    //Unit: metre^(2)
    _class_parameter[95] = PFILE(277) * 1e-12;
    //A_APC, mwa21d667a_6a85_4ae2_bc87_c7697d55781e, index: 96
    //Unit: metre^(2)
    _class_parameter[96] = PFILE(278) * 1e-12;
    //k_M1p0_TCR_on, mw6f110553_db6e_4780_bcfc_efadc54f4328, index: 97
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[97] = PFILE(279) * 602214199000.0;
    //k_M1p0_TCR_off, mw1e699e4e_46dc_4669_bce0_b4de73989cc7, index: 98
    //Unit: second^(-1)
    _class_parameter[98] = PFILE(280) * 1.0;
    //k_M1p0_TCR_p, mw5b165628_9a99_45d9_ac7a_74b3ff589f26, index: 99
    //Unit: second^(-1)
    _class_parameter[99] = PFILE(281) * 1.0;
    //phi_M1p0_TCR, mwca7f29b4_c91d_4303_a4bf_2dce64444360, index: 100
    //Unit: second^(-1)
    _class_parameter[100] = PFILE(282) * 1.0;
    //N_M1p0_TCR, mwc25cfe94_3e5f_4296_a5e6_d6a03cafde0b, index: 101
    //Unit: dimensionless^(1)
    _class_parameter[101] = PFILE(283) * 1.0;
    //TCR_p0_tot, mw6133f052_3f62_4069_b8e2_cf2b483aab2e, index: 102
    //Unit: metre^(-2)mole^(1)
    _class_parameter[102] = PFILE(284) * 1.66053872801495e-12;
    //k_P1_up, mwc0e65e38_0938_471e_b747_f34f6a9e22f6, index: 103
    //Unit: mole^(-1)second^(-1)
    _class_parameter[103] = PFILE(286) * 6.97007174768519e+18;
    //k_xP1_deg, mw4792e234_80c6_4194_b2fc_be3206990758, index: 104
    //Unit: second^(-1)
    _class_parameter[104] = PFILE(287) * 1.15740740740741e-05;
    //k_P1_deg, mw2dedbcdb_e234_46a7_988e_a0dbf01f1697, index: 105
    //Unit: second^(-1)
    _class_parameter[105] = PFILE(288) * 1.15740740740741e-05;
    //k_p1_deg, mwe80e7b2f_a122_42f5_827f_518260093198, index: 106
    //Unit: second^(-1)
    _class_parameter[106] = PFILE(289) * 1.15740740740741e-05;
    //k_P1_on, mwdf8d0dce_cf49_4772_84e1_2bf86c9ab84c, index: 107
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[107] = PFILE(290) * 1.1574074074074112e-08;
    //k_P1_d1, mw1f8a78c1_7249_4604_8575_56c7d0fc3e93, index: 108
    //Unit: metre^(-3)mole^(1)
    _class_parameter[108] = PFILE(291) * 999.9999999999994;
    //p1_50, mw4de0d283_e3eb_4beb_a466_6bc54437ed0a, index: 109
    //Unit: metre^(-2)mole^(1)
    _class_parameter[109] = PFILE(292) * 1.66053872801495e-12;
    //P1_C1, mwd21a702e_802a_464b_8dc1_26b1596516d6, index: 110
    //Unit: dimensionless^(1)
    _class_parameter[110] = PFILE(293) * 6.02214199e+23;
    //P1_C2, mw14953708_f27a_4aab_8c93_95dcf45f3625, index: 111
    //Unit: dimensionless^(1)
    _class_parameter[111] = PFILE(294) * 6.02214199e+23;
    //k_M1p1_TCR_on, mwe13c623c_8941_413f_b9ea_4ed29aa74c2b, index: 112
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[112] = PFILE(295) * 602214199000.0;
    //k_M1p1_TCR_off, mwf3689e73_e8ac_4c5f_b812_dba72320c4d4, index: 113
    //Unit: second^(-1)
    _class_parameter[113] = PFILE(296) * 1.0;
    //k_M1p1_TCR_p, mwef3e4a1b_d1d3_4c0f_aea3_0c368f8fdc4a, index: 114
    //Unit: second^(-1)
    _class_parameter[114] = PFILE(297) * 1.0;
    //phi_M1p1_TCR, mwc0bb5137_49c6_4247_af7f_01f15ded629f, index: 115
    //Unit: second^(-1)
    _class_parameter[115] = PFILE(298) * 1.0;
    //N_M1p1_TCR, mw867eb046_0db4_434a_8a74_e89542db4f6a, index: 116
    //Unit: dimensionless^(1)
    _class_parameter[116] = PFILE(299) * 1.0;
    //TCR_p1_tot, mw0579ff47_9e57_4ecb_b005_f350b26f56d1, index: 117
    //Unit: metre^(-2)mole^(1)
    _class_parameter[117] = PFILE(300) * 1.66053872801495e-12;
    //q_P_aPD1, mw0d71978a_5d9e_4b91_9ce1_f801800b5634, index: 118
    //Unit: metre^(3)second^(-1)
    _class_parameter[118] = PFILE(302) * 0.0010000000000000007;
    //q_T_aPD1, mw247a2d62_966c_4345_ae0d_f15ff8cf1171, index: 119
    //Unit: metre^(3)second^(-1)
    _class_parameter[119] = PFILE(303) * 1.0000000000000006e-06;
    //q_LN_aPD1, mwd043eb23_1acc_467a_a4e8_d6c1d7cc2e6c, index: 120
    //Unit: metre^(3)second^(-1)
    _class_parameter[120] = PFILE(304) * 1.0000000000000006e-06;
    //q_LD_aPD1, mwb56ec0c0_6589_4b61_b243_d7263dfc0cb5, index: 121
    //Unit: second^(-1)
    _class_parameter[121] = PFILE(305) * 0.0166666666666667;
    //k_cl_aPD1, mw07352214_68b1_415c_8f0b_738acc5eb5d0, index: 122
    //Unit: metre^(3)second^(-1)
    _class_parameter[122] = PFILE(306) * 1.1574074074074112e-08;
    //gamma_C_aPD1, mwf325e503_79e9_4e88_9953_8a423e251038, index: 123
    //Unit: dimensionless^(1)
    _class_parameter[123] = PFILE(307) * 1.0;
    //gamma_P_aPD1, mwb57b4cc2_4b34_4037_96ac_17388e4ce29f, index: 124
    //Unit: dimensionless^(1)
    _class_parameter[124] = PFILE(308) * 1.0;
    //gamma_T_aPD1, mw74dedd98_91ec_47e1_b0d7_1d9d572842b2, index: 125
    //Unit: dimensionless^(1)
    _class_parameter[125] = PFILE(309) * 1.0;
    //gamma_LN_aPD1, mweb796bb9_5bfd_46f2_b1ba_086c622d84dd, index: 126
    //Unit: dimensionless^(1)
    _class_parameter[126] = PFILE(310) * 1.0;
    //q_P_aPDL1, mwd90bc5a4_45c2_4d21_a9b0_4934d36094ae, index: 127
    //Unit: metre^(3)second^(-1)
    _class_parameter[127] = PFILE(311) * 0.0010000000000000007;
    //q_T_aPDL1, mw63c9a357_b78f_4667_93e3_93d606c2b068, index: 128
    //Unit: metre^(3)second^(-1)
    _class_parameter[128] = PFILE(312) * 1.0000000000000006e-06;
    //q_LN_aPDL1, mw8fb3b848_03a3_4bbd_ad2c_3b556e941560, index: 129
    //Unit: metre^(3)second^(-1)
    _class_parameter[129] = PFILE(313) * 1.0000000000000006e-06;
    //q_LD_aPDL1, mw2bf884b3_637d_4222_9802_f44c33ad24c2, index: 130
    //Unit: second^(-1)
    _class_parameter[130] = PFILE(314) * 0.0166666666666667;
    //k_cl_aPDL1, mw57496210_dbbd_472b_96eb_feb7147bc271, index: 131
    //Unit: metre^(3)second^(-1)
    _class_parameter[131] = PFILE(315) * 1.1574074074074112e-08;
    //gamma_C_aPDL1, mw1ef08a8c_5e3f_438c_b11a_93034c414798, index: 132
    //Unit: dimensionless^(1)
    _class_parameter[132] = PFILE(316) * 1.0;
    //gamma_P_aPDL1, mwfd2f8ec6_15b1_4996_9406_19934897c6a2, index: 133
    //Unit: dimensionless^(1)
    _class_parameter[133] = PFILE(317) * 1.0;
    //gamma_T_aPDL1, mw51dac9a1_c586_44cb_9200_b32d4621bac1, index: 134
    //Unit: dimensionless^(1)
    _class_parameter[134] = PFILE(318) * 1.0;
    //gamma_LN_aPDL1, mw61de9dda_5fb9_431b_a9ef_cdc97a38e39b, index: 135
    //Unit: dimensionless^(1)
    _class_parameter[135] = PFILE(319) * 1.0;
    //k_cln_aPDL1, mw2fe353be_4375_4bf4_ad7e_066996d9ee7f, index: 136
    //Unit: mole^(1)second^(-1)
    _class_parameter[136] = PFILE(320) * 1.15740740740741e-14;
    //Kc_aPDL1, mw4e30351e_4108_456c_920b_8ce65088a703, index: 137
    //Unit: metre^(-3)mole^(1)
    _class_parameter[137] = PFILE(321) * 1.0000000000000008e-06;
    //q_P_aCTLA4, mw93b499a3_020e_4635_9757_5da00c6eeffd, index: 138
    //Unit: metre^(3)second^(-1)
    _class_parameter[138] = PFILE(322) * 0.0010000000000000007;
    //q_T_aCTLA4, mwfaf46650_ede9_433a_8614_c73124d770aa, index: 139
    //Unit: metre^(3)second^(-1)
    _class_parameter[139] = PFILE(323) * 1.0000000000000006e-06;
    //q_LN_aCTLA4, mwbfa771a1_bdae_47dc_a37e_ad142f88dc56, index: 140
    //Unit: metre^(3)second^(-1)
    _class_parameter[140] = PFILE(324) * 1.0000000000000006e-06;
    //q_LD_aCTLA4, mw0a5a6120_c89e_4900_8849_886fccb3162b, index: 141
    //Unit: second^(-1)
    _class_parameter[141] = PFILE(325) * 0.0166666666666667;
    //k_cl_aCTLA4, mw03195f74_105c_40b3_a8de_1d1d35456094, index: 142
    //Unit: metre^(3)second^(-1)
    _class_parameter[142] = PFILE(326) * 1.1574074074074112e-08;
    //gamma_C_aCTLA4, mw50130d83_ba85_4b69_aa3e_9b3e1bfe44cf, index: 143
    //Unit: dimensionless^(1)
    _class_parameter[143] = PFILE(327) * 1.0;
    //gamma_P_aCTLA4, mw85de870d_ada0_4c2f_9b00_b43cae7b7141, index: 144
    //Unit: dimensionless^(1)
    _class_parameter[144] = PFILE(328) * 1.0;
    //gamma_T_aCTLA4, mwa5aab5f0_306d_4b35_93d1_9a89c295b927, index: 145
    //Unit: dimensionless^(1)
    _class_parameter[145] = PFILE(329) * 1.0;
    //gamma_LN_aCTLA4, mwe66d6119_cd9b_4cd0_80dd_c7f92c96e2d6, index: 146
    //Unit: dimensionless^(1)
    _class_parameter[146] = PFILE(330) * 1.0;
    //kon_PD1_PDL1, mwff66f420_3847_4e07_8078_2b7bc4358873, index: 147
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[147] = PFILE(331) * 1000000000000.0;
    //k_out_PDL1, mw334268ec_dcd0_49ae_b851_ee3cd478d44b, index: 148
    //Unit: mole^(1)second^(-1)
    _class_parameter[148] = PFILE(332) * 1.92191982409137e-29;
    //k_in_PDL1, mwfec4fb28_ee8d_4bab_adf7_7f834e4c6b75, index: 149
    //Unit: second^(-1)
    _class_parameter[149] = PFILE(333) * 1.15740740740741e-05;
    //r_PDL1_IFNg, mw56d0a446_cc2e_4720_9b27_7980bd7a43bd, index: 150
    //Unit: dimensionless^(1)
    _class_parameter[150] = PFILE(334) * 1.0;
    //kon_PD1_PDL2, mw17a44964_e36f_4fcf_80bb_3f78967f899a, index: 151
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[151] = PFILE(335) * 1000000000000.0;
    //kon_PD1_aPD1, mw76081873_87f8_4895_aafe_8078e3c9f067, index: 152
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[152] = PFILE(336) * 0.0010000000000000007;
    //kon_PDL1_aPDL1, mw18562cee_4661_4a65_b5e1_a1daa8c2c42e, index: 153
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[153] = PFILE(337) * 0.0010000000000000007;
    //kon_CD28_CD80, mwfac3709c_96be_429f_bca4_d1fc835858f0, index: 154
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[154] = PFILE(338) * 1000000000000.0;
    //kon_CD28_CD86, mw4e8d3473_cdde_4c35_9b76_ce2b0f67b8e3, index: 155
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[155] = PFILE(339) * 1000000000000.0;
    //kon_CTLA4_CD80, mw310ad742_1114_44f8_8c9c_4937df550c34, index: 156
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[156] = PFILE(340) * 1000000000000.0;
    //kon_CTLA4_CD86, mw9880ebd5_c4ae_44df_8fdc_0a225bace74e, index: 157
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[157] = PFILE(341) * 1000000000000.0;
    //kon_CD80_PDL1, mwd5915e36_4e17_4915_ba39_796a7e07ee18, index: 158
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[158] = PFILE(342) * 1000000000000.0;
    //kon_CTLA4_aCTLA4, mwdd058bcb_f321_455a_b6bc_e1272b299026, index: 159
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[159] = PFILE(343) * 0.0010000000000000007;
    //kon_CD80_CD80, mw01ad60af_72b5_4ce5_8526_d778ab7f29fe, index: 160
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[160] = PFILE(344) * 1000000000000.0;
    //koff_PD1_PDL1, mwabd26b6e_42fb_48b9_9a61_c0683c2ded46, index: 161
    //Unit: second^(-1)
    _class_parameter[161] = PFILE(345) * 1.0;
    //koff_PD1_PDL2, mw123de0ae_b4c3_442d_8396_58bf6ccadaff, index: 162
    //Unit: second^(-1)
    _class_parameter[162] = PFILE(346) * 1.0;
    //koff_PD1_aPD1, mw537f9c03_e8de_440a_a2fb_d2e26f56dac7, index: 163
    //Unit: second^(-1)
    _class_parameter[163] = PFILE(347) * 1.0;
    //koff_PDL1_aPDL1, mwf0cfb962_7f15_4afc_8f23_58e4a79ed194, index: 164
    //Unit: second^(-1)
    _class_parameter[164] = PFILE(348) * 1.0;
    //koff_CD28_CD80, mw0a038a15_7972_468e_9fce_38d474af85da, index: 165
    //Unit: second^(-1)
    _class_parameter[165] = PFILE(349) * 1.0;
    //koff_CD28_CD86, mw01e837b1_34f8_4466_8a75_e46e48143333, index: 166
    //Unit: second^(-1)
    _class_parameter[166] = PFILE(350) * 1.0;
    //koff_CTLA4_CD80, mw55b4e3df_1364_4b0d_bdb9_9400e44b45dc, index: 167
    //Unit: second^(-1)
    _class_parameter[167] = PFILE(351) * 1.0;
    //koff_CTLA4_CD86, mw2553cfa1_c025_4a8a_916f_c34b247591f5, index: 168
    //Unit: second^(-1)
    _class_parameter[168] = PFILE(352) * 1.0;
    //koff_CD80_PDL1, mw689ceae3_5ae0_4258_ad3e_e245b1d24522, index: 169
    //Unit: second^(-1)
    _class_parameter[169] = PFILE(353) * 1.0;
    //koff_CTLA4_aCTLA4, mw7ae90c13_f671_4f6d_b3a6_a05bb2bbeec3, index: 170
    //Unit: second^(-1)
    _class_parameter[170] = PFILE(354) * 1.0;
    //koff_CD80_CD80, mw8447f3bf_9fe0_42ed_88a3_bf6bf5c548c5, index: 171
    //Unit: second^(-1)
    _class_parameter[171] = PFILE(355) * 1.0;
    //Chi_PD1_aPD1, mw172d818c_96e9_4e4c_8fab_29818288527d, index: 172
    //Unit: metre^(-1)
    _class_parameter[172] = PFILE(356) * 999999999.9999999;
    //Chi_PDL1_aPDL1, mw20993cea_bc5f_48c4_8b69_823c1cbf76aa, index: 173
    //Unit: metre^(-1)
    _class_parameter[173] = PFILE(357) * 999999999.9999999;
    //Chi_CTLA4_aCTLA4, mw8aad30b3_82c6_4170_805e_0a483f3f6e54, index: 174
    //Unit: metre^(-1)
    _class_parameter[174] = PFILE(358) * 999999999.9999999;
    //PD1_50, mw97fdab53_85c9_4c20_92f4_998d3eeec3ce, index: 175
    //Unit: metre^(-2)mole^(1)
    _class_parameter[175] = PFILE(359) * 1.66053872801495e-12;
    //n_PD1, mw82e76e3a_475a_4c9f_b336_bf520eb1fac2, index: 176
    //Unit: dimensionless^(1)
    _class_parameter[176] = PFILE(360) * 1.0;
    //CD28_CD8X_50, mw76d6e8d1_0e20_4458_b6b5_f59932c8104d, index: 177
    //Unit: metre^(-2)mole^(1)
    _class_parameter[177] = PFILE(361) * 1.66053872801495e-12;
    //n_CD28_CD8X, mw30014f58_209b_41ab_be10_79dd7b10b444, index: 178
    //Unit: dimensionless^(1)
    _class_parameter[178] = PFILE(362) * 1.0;
    //T_PD1_total, mwfd8cd76f_c493_49eb_854b_baf010851b09, index: 179
    //Unit: mole^(1)
    _class_parameter[179] = PFILE(363) * 1.66053872801495e-24;
    //T_CD28_total, mw8163b229_f823_40f4_b19d_aedcafaba33b, index: 180
    //Unit: mole^(1)
    _class_parameter[180] = PFILE(364) * 1.66053872801495e-24;
    //T_CTLA4_syn, mw9b85daec_15ff_45d7_b024_fca433c7cfdd, index: 181
    //Unit: mole^(1)
    _class_parameter[181] = PFILE(365) * 1.66053872801495e-24;
    //T_PDL1_total, mwfe730607_a906_4687_b6d2_9f6e0ada65ae, index: 182
    //Unit: mole^(1)
    _class_parameter[182] = PFILE(366) * 1.66053872801495e-24;
    //C1_PDL1_base, mw379cda8b_7d4d_491b_8d9c_f903e9b1e1e6, index: 183
    //Unit: mole^(1)
    _class_parameter[183] = PFILE(367) * 1.66053872801495e-24;
    //r_PDL2C1, mw71bdf5d9_25d9_4986_81f0_ff22aea2f66e, index: 184
    //Unit: dimensionless^(1)
    _class_parameter[184] = PFILE(368) * 1.0;
    //C1_CD80_total, mw9e2b686d_b0be_4d64_8d03_0c652c9b722f, index: 185
    //Unit: mole^(1)
    _class_parameter[185] = PFILE(369) * 1.66053872801495e-24;
    //C1_CD86_total, mw18176e4f_0aeb_437f_baf6_869e63079e9d, index: 186
    //Unit: mole^(1)
    _class_parameter[186] = PFILE(370) * 1.66053872801495e-24;
    //APC_PDL1_base, mwa6775827_434f_415a_9e1d_b3fd1778b65d, index: 187
    //Unit: mole^(1)
    _class_parameter[187] = PFILE(373) * 1.66053872801495e-24;
    //r_PDL2APC, mw6d5be92a_2eb3_40a2_aae9_b479c79272e7, index: 188
    //Unit: dimensionless^(1)
    _class_parameter[188] = PFILE(374) * 1.0;
    //APC_CD80_total, mwd45f9369_d639_4b74_9e4e_a751744167a2, index: 189
    //Unit: mole^(1)
    _class_parameter[189] = PFILE(375) * 1.66053872801495e-24;
    //APC_CD86_total, mw14a57954_3422_4637_808d_f25fab8090f2, index: 190
    //Unit: mole^(1)
    _class_parameter[190] = PFILE(376) * 1.66053872801495e-24;
    //k_Th_act, mw2ba47ae2_b227_40c5_acea_6688d2ea96be, index: 191
    //Unit: second^(-1)
    _class_parameter[191] = PFILE(379) * 1.15740740740741e-05;
    //k_Th_Treg, mwe6cccfc9_8812_4e06_bec3_e1a033bb8e9f, index: 192
    //Unit: second^(-1)
    _class_parameter[192] = PFILE(380) * 1.15740740740741e-05;
    //k_TGFb_Tsec, mwf9c1e4d1_0abc_4a71_8c12_8351382ddeaf, index: 193
    //Unit: second^(-1)
    _class_parameter[193] = PFILE(381) * 6970071747.68519;
    //k_TGFb_deg, mw308d09f9_04ae_406a_8703_c58d23d478d7, index: 194
    //Unit: second^(-1)
    _class_parameter[194] = PFILE(382) * 1.15740740740741e-05;
    //TGFb_50, mw11ecc560_923a_4d65_9dca_e14c806950c6, index: 195
    //Unit: metre^(-3)mole^(1)
    _class_parameter[195] = PFILE(383) * 1.0000000000000008e-06;
    //TGFb_50_Teff, mwccd63c25_5a26_4d79_bd2f_96b5f220c1d2, index: 196
    //Unit: metre^(-3)mole^(1)
    _class_parameter[196] = PFILE(384) * 1.0000000000000008e-06;
    //Kc_rec, mwe0a7df87_114e_4668_a65f_300a12842847, index: 197
    //Unit: mole^(2)
    _class_parameter[197] = PFILE(385) * 2.757388867237499e-48;
    //TGFbase, mw59028808_9688_420c_a590_848806223bf1, index: 198
    //Unit: metre^(-3)mole^(1)
    _class_parameter[198] = PFILE(386) * 1.0000000000000008e-06;
    //k_IFNg_Thsec, mwb5708e4b_0b95_4985_b1e2_28acef5c0ce0, index: 199
    //Unit: second^(-1)
    _class_parameter[199] = PFILE(387) * 1.15740740740741e-05;
    //k_IFNg_deg, mw229ec9c0_7aaa_4ff2_9347_86a3ddaa786c, index: 200
    //Unit: second^(-1)
    _class_parameter[200] = PFILE(388) * 1.15740740740741e-05;
    //IFNg_50_ind, mw5af13b97_66ac_496c_93f8_9cd5842908b9, index: 201
    //Unit: metre^(-3)mole^(1)
    _class_parameter[201] = PFILE(389) * 1.0000000000000013e-09;
    //k_CCL2_sec, mw7c8f4495_3016_47e8_8768_b1cd1205d56a, index: 202
    //Unit: second^(-1)
    _class_parameter[202] = PFILE(393) * 6970071747.68519;
    //k_CCL2_deg, mwd5bc6f41_5176_4646_949d_0d8a31ae65fe, index: 203
    //Unit: second^(-1)
    _class_parameter[203] = PFILE(394) * 0.000277777777777778;
    //CCL2_50, mw93c5f5d4_638f_4a78_93b8_ac04ac87ede8, index: 204
    //Unit: metre^(-3)mole^(1)
    _class_parameter[204] = PFILE(395) * 1.0000000000000008e-06;
    //k_MDSC_rec, mwa98ef6ca_6465_400d_8647_508cd775812c, index: 205
    //Unit: metre^(-3)mole^(1)second^(-1)
    _class_parameter[205] = PFILE(396) * 1.9219198240913757e-23;
    //k_MDSC_death, mw5b261931_758c_4f2f_b8ba_da94752e20dc, index: 206
    //Unit: second^(-1)
    _class_parameter[206] = PFILE(397) * 1.15740740740741e-05;
    //k_NO_deg, mw11903951_c66d_4bdf_98dd_b7fd761019b2, index: 207
    //Unit: second^(-1)
    _class_parameter[207] = PFILE(398) * 1.15740740740741e-05;
    //k_ArgI_deg, mw70fa49b8_8bc0_4cf9_9403_363daae48887, index: 208
    //Unit: second^(-1)
    _class_parameter[208] = PFILE(399) * 1.15740740740741e-05;
    //k_NO_sec, mw6d0d7c1b_7a87_4152_aa0e_5b0a5bc60384, index: 209
    //Unit: second^(-1)
    _class_parameter[209] = PFILE(400) * 6970071747.68519;
    //k_ArgI_sec, mwc7175d78_d1fd_41f9_a12b_9983ca949f50, index: 210
    //Unit: second^(-1)
    _class_parameter[210] = PFILE(401) * 6970071747.68519;
    //ArgI_50_Teff, mw825a6198_050f_4137_9e73_3a35c4496a84, index: 211
    //Unit: metre^(-3)mole^(1)
    _class_parameter[211] = PFILE(402) * 1.0000000000000008e-06;
    //NO_50_Teff, mwa669730d_b8a1_46b4_837a_b4bddda42102, index: 212
    //Unit: metre^(-3)mole^(1)
    _class_parameter[212] = PFILE(403) * 1.0000000000000008e-06;
    //ArgI_50_Treg, mw0b1f35ad_7c9a_4b55_a619_5e5988f8b03e, index: 213
    //Unit: metre^(-3)mole^(1)
    _class_parameter[213] = PFILE(404) * 1.0000000000000008e-06;
    //k_a1_cabozantinib, mw0fdca9f0_033f_4ddf_a908_09f09bf91bdb, index: 214
    //Unit: second^(-1)
    _class_parameter[214] = PFILE(409) * 0.000277777777777778;
    //k_a2_cabozantinib, mw7bef680a_da0d_4f86_83ed_e863a8ff5143, index: 215
    //Unit: second^(-1)
    _class_parameter[215] = PFILE(410) * 0.000277777777777778;
    //k_cln_cabozantinib, mw77f8fe1d_68b0_4a58_b99c_e2fac8d3c5a8, index: 216
    //Unit: metre^(-3)mole^(1)second^(-1)
    _class_parameter[216] = PFILE(411) * 0.277777777777778;
    //Kc_cabozantinib, mw14c96337_3255_45cb_9532_fec809d13715, index: 217
    //Unit: metre^(-3)mole^(1)
    _class_parameter[217] = PFILE(412) * 999.9999999999994;
    //lagP1_cabozantinib, mwfbc76909_9521_4f71_b029_48c761b39cd0, index: 218
    //Unit: second^(1)
    _class_parameter[218] = PFILE(413) * 3600.0;
    //lagP2_cabozantinib, mwd72e0587_93e5_40eb_b048_c732d6eb6a61, index: 219
    //Unit: second^(1)
    _class_parameter[219] = PFILE(414) * 3600.0;
    //F_cabozantinib, mw4624d04f_eef7_4e9c_8731_c76b0b62f034, index: 220
    //Unit: dimensionless^(1)
    _class_parameter[220] = PFILE(415) * 1.0;
    //q_P_cabozantinib, mwd623a2b1_7e33_403c_89e0_59160aeb0177, index: 221
    //Unit: second^(-1)
    _class_parameter[221] = PFILE(416) * 1.0;
    //q_T_cabozantinib, mw196e8df6_69ce_451c_9dc9_bab9cb6601f2, index: 222
    //Unit: second^(-1)
    _class_parameter[222] = PFILE(417) * 1.0;
    //q_LN_cabozantinib, mw808483ee_c4c9_4421_bbea_844befc834a6, index: 223
    //Unit: second^(-1)
    _class_parameter[223] = PFILE(418) * 1.0;
    //q_LD_cabozantinib, mw305d9a04_5de6_4e25_9b61_1dfb59ae0c51, index: 224
    //Unit: second^(-1)
    _class_parameter[224] = PFILE(419) * 1.0;
    //gamma_C_cabozantinib, mw46a08d93_7d41_4c9c_a7a8_204944ab5eaa, index: 225
    //Unit: dimensionless^(1)
    _class_parameter[225] = PFILE(420) * 1.0;
    //gamma_P_cabozantinib, mwa6aecb21_6944_4637_8c8e_269f7b2bf270, index: 226
    //Unit: dimensionless^(1)
    _class_parameter[226] = PFILE(421) * 1.0;
    //gamma_T_cabozantinib, mwd4c74c65_50b8_49c9_9753_d571e441b65d, index: 227
    //Unit: dimensionless^(1)
    _class_parameter[227] = PFILE(422) * 1.0;
    //gamma_LN_cabozantinib, mwe9850eb4_e8d5_4e3c_9d68_3da471a5d387, index: 228
    //Unit: dimensionless^(1)
    _class_parameter[228] = PFILE(423) * 1.0;
    //IC50_MET, mwc34b3d62_4989_4558_a923_bff87f4a3c68, index: 229
    //Unit: metre^(-3)mole^(1)
    _class_parameter[229] = PFILE(424) * 1.0000000000000008e-06;
    //IC50_RET, mwfb83812f_9f6d_424d_a488_a8c81ac0d308, index: 230
    //Unit: metre^(-3)mole^(1)
    _class_parameter[230] = PFILE(425) * 1.0000000000000008e-06;
    //IC50_AXL, mwa4eeb876_99a0_4eab_847e_de5b43f7ff37, index: 231
    //Unit: metre^(-3)mole^(1)
    _class_parameter[231] = PFILE(426) * 1.0000000000000008e-06;
    //IC50_VEGFR2, mwcc105608_4b1a_48e9_902b_5119a489f2ca, index: 232
    //Unit: metre^(-3)mole^(1)
    _class_parameter[232] = PFILE(427) * 1.0000000000000008e-06;
    //k_C_resist, mw42e5ca97_bc06_4c44_9416_3df300b8e815, index: 233
    //Unit: second^(-1)
    _class_parameter[233] = PFILE(428) * 1.15740740740741e-05;
    //k_K_cabo, mw1b364385_3d67_4b90_8aac_17a2ffccd4b0, index: 234
    //Unit: dimensionless^(1)
    _class_parameter[234] = PFILE(429) * 1.0;
    //k_Mac_rec, mwe3c99b9f_f5eb_42d1_8f55_ef5e2e19ba82, index: 235
    //Unit: metre^(-3)mole^(1)second^(-1)
    _class_parameter[235] = PFILE(432) * 1.9219198240913757e-23;
    //k_Mac_death, mw99c3bef5_362c_49d2_990e_f96a624c2d1d, index: 236
    //Unit: second^(-1)
    _class_parameter[236] = PFILE(433) * 1.15740740740741e-05;
    //k_TGFb_Msec, mw93eb49f4_cc26_4b79_987c_af4fcc8e01ad, index: 237
    //Unit: second^(-1)
    _class_parameter[237] = PFILE(434) * 6970071747.68519;
    //k_vas_Msec, mwc39bef2c_131b_4115_b8b0_46e7f837341d, index: 238
    //Unit: kilogram^(1)mole^(-1)second^(-1)
    _class_parameter[238] = PFILE(435) * 6970.0717476851905;
    //k_IL12_sec, mwd3b75f6f_5dab_4844_b9ed_d6f6eac05b92, index: 239
    //Unit: second^(-1)
    _class_parameter[239] = PFILE(436) * 6970071747.68519;
    //k_IL12_Msec, mwa5d91037_c49a_4a70_be6f_cbebb432e82a, index: 240
    //Unit: second^(-1)
    _class_parameter[240] = PFILE(437) * 6970071747.68519;
    //k_IL12_deg, mw85c991a8_1f5f_49b0_bba3_474c7536f09b, index: 241
    //Unit: second^(-1)
    _class_parameter[241] = PFILE(438) * 0.000277777777777778;
    //k_IL10_sec, mw3d50483d_0c91_4ca2_a070_7fbd0449d0ea, index: 242
    //Unit: second^(-1)
    _class_parameter[242] = PFILE(439) * 6970071747.68519;
    //k_IL10_deg, mw824b99a1_90d1_41ff_9cb8_4ec86f1b7bc3, index: 243
    //Unit: second^(-1)
    _class_parameter[243] = PFILE(440) * 1.15740740740741e-05;
    //k_M2_pol, mwffe3614c_2f62_44a3_aad1_4d1059c58641, index: 244
    //Unit: second^(-1)
    _class_parameter[244] = PFILE(441) * 1.15740740740741e-05;
    //k_M1_pol, mw20615208_8448_40a3_8951_b950a4ab0205, index: 245
    //Unit: second^(-1)
    _class_parameter[245] = PFILE(442) * 1.15740740740741e-05;
    //IL10_50, mw3096917c_12fc_4e48_916f_104071254a3b, index: 246
    //Unit: metre^(-3)mole^(1)
    _class_parameter[246] = PFILE(443) * 1.0000000000000013e-09;
    //IL12_50, mwae0e24f3_e335_43b8_989d_3a5ebe65e7af, index: 247
    //Unit: metre^(-3)mole^(1)
    _class_parameter[247] = PFILE(444) * 1.0000000000000013e-09;
    //IFNg_50, mw0943861c_347e_45d6_92ea_69b74b18edea, index: 248
    //Unit: metre^(-3)mole^(1)
    _class_parameter[248] = PFILE(445) * 1.0000000000000013e-09;
    //k_M1_phago, mw8cc92bd5_61f7_488e_94c4_8ae080576ea9, index: 249
    //Unit: second^(-1)
    _class_parameter[249] = PFILE(446) * 1.15740740740741e-05;
    //vol_Mcell, mw8da8eca5_bf10_4b48_aa02_cf70d507a74f, index: 250
    //Unit: metre^(3)mole^(-1)
    _class_parameter[250] = PFILE(447) * 602214.1989999996;
    //kon_CD47_SIRPa, mw10c7a7aa_b286_4413_991a_642a31a5e29c, index: 251
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[251] = PFILE(448) * 16666666666.6667;
    //koff_CD47_SIRPa, mw5af6e8a6_10b1_4816_a7c0_940f9a30bc9a, index: 252
    //Unit: second^(-1)
    _class_parameter[252] = PFILE(449) * 0.0166666666666667;
    //SIRPa_50, mwb5f0663a_04fb_415c_ac0e_540277b94486, index: 253
    //Unit: metre^(-2)mole^(1)
    _class_parameter[253] = PFILE(450) * 1.66053872801495e-12;
    //n_SIRPa, mw6b436961_8d90_4c37_a466_4816fafd02e3, index: 254
    //Unit: dimensionless^(1)
    _class_parameter[254] = PFILE(451) * 1.0;
    //C_CD47, mwc8dc1861_9e27_49c0_80c3_e7a7fbd2dac5, index: 255
    //Unit: metre^(-2)mole^(1)
    _class_parameter[255] = PFILE(455) * 1.66053872801495e-12;
    //M_PD1_total, mwa778f883_d302_47ef_8955_91a8f86af530, index: 256
    //Unit: mole^(1)
    _class_parameter[256] = PFILE(456) * 1.66053872801495e-24;
    //M_SIRPa, mw7f451eee_0516_4f28_aa46_bb46c3b59326, index: 257
    //Unit: metre^(-2)mole^(1)
    _class_parameter[257] = PFILE(457) * 1.66053872801495e-12;
    //A_Mcell, mw5a3e1746_8a76_4777_9b3b_56cba45b82fa, index: 258
    //Unit: metre^(2)
    _class_parameter[258] = PFILE(458) * 1e-12;
    //IL10_50_phago, mw8b19f0da_183d_47cc_9cac_3a410ae5d37c, index: 259
    //Unit: metre^(-3)mole^(1)
    _class_parameter[259] = PFILE(459) * 1.0000000000000013e-09;
    //K_Mac_C, mw2d96c898_143b_4523_b4b0_4b4e460436b3, index: 260
    //Unit: dimensionless^(1)
    _class_parameter[260] = PFILE(460) * 1.0;
    //k_fib_rec, mw51f28e7c_0a49_417a_9f91_965b4ad8cf74, index: 261
    //Unit: metre^(-3)mole^(1)second^(-1)
    _class_parameter[261] = PFILE(465) * 1.9219198240913757e-23;
    //k_fib_const, mwba602a18_1e5f_42d0_a857_38494db62ad8, index: 262
    //Unit: metre^(-3)mole^(1)
    _class_parameter[262] = PFILE(466) * 1.6605387280149534e-18;
    //k_caf_tran, mw96386e4f_5247_4936_b4ff_417b2335cd4b, index: 263
    //Unit: second^(-1)
    _class_parameter[263] = PFILE(467) * 1.15740740740741e-05;
    //k_ECM_fib_sec, mwb994b081_1a77_4fee_b815_e7cbf71eab99, index: 264
    //Unit: second^(-1)
    _class_parameter[264] = PFILE(468) * 6970071747.68519;
    //k_ECM_CAF_sec, mwc676eb61_706c_4b4f_b50c_c3e73c0a4976, index: 265
    //Unit: second^(-1)
    _class_parameter[265] = PFILE(469) * 6970071747.68519;
    //k_ECM_deg, mwe122b976_daf5_4f12_b1d3_bd5bdc880541, index: 266
    //Unit: second^(-1)
    _class_parameter[266] = PFILE(470) * 1.15740740740741e-05;
    //ECM_base, mw3c8aa3ec_7775_40e0_a184_7ca813637151, index: 267
    //Unit: metre^(-3)mole^(1)
    _class_parameter[267] = PFILE(471) * 1.0000000000000008e-06;
    //ECM_max, mwef227d48_a64b_4a09_a0ea_ad1b46c900ad, index: 268
    //Unit: metre^(-3)mole^(1)
    _class_parameter[268] = PFILE(472) * 1.0000000000000008e-06;
    //ECM_MW, mw2ff6cf30_72e7_4ba9_b8d5_b8c5dbfe2c23, index: 269
    //Unit: kilogram^(1)mole^(-1)
    _class_parameter[269] = PFILE(474) * 0.001;
    //ECM_density, mw099a3e0b_c07a_4410_a51a_f549518d1304, index: 270
    //Unit: kilogram^(1)metre^(-3)
    _class_parameter[270] = PFILE(475) * 1000.0;
    //k_fib_death, mwc6b50f08_639b_4031_8979_7f36162ce302, index: 271
    //Unit: second^(-1)
    _class_parameter[271] = PFILE(476) * 1.15740740740741e-05;
    //k_CAF_death, mw22574ee9_f67b_4954_bc67_2551b0d09a2a, index: 272
    //Unit: second^(-1)
    _class_parameter[272] = PFILE(477) * 1.15740740740741e-05;
    //ECM_50_T_exh, mw808ab991_b45b_4680_9511_bcf1bdcde036, index: 273
    //Unit: metre^(-3)mole^(1)
    _class_parameter[273] = PFILE(478) * 1.0000000000000008e-06;
    //ECM_50_T_mot, mw44abfe67_22a8_4c34_ae30_ffb1fc746613, index: 274
    //Unit: metre^(-3)mole^(1)
    _class_parameter[274] = PFILE(479) * 1.0000000000000008e-06;
    //vol_Fibcell, mw46e575ad_032b_4207_807c_36cdfb715f71, index: 275
    //Unit: metre^(3)mole^(-1)
    _class_parameter[275] = PFILE(480) * 602214.1989999996;
    //vol_CAFcell, mw94db9b8c_8464_4249_a359_b0e36f7ee8a2, index: 276
    //Unit: metre^(3)mole^(-1)
    _class_parameter[276] = PFILE(481) * 602214.1989999996;
}

void ODE_system::setupVariables(void){

    _species_var = std::vector<realtype>(153, 0);
    _nonspecies_var = std::vector<realtype>(0, 0);
    //species not part of ode left-hand side
    _species_other =  std::vector<realtype>(2, 0);
    
    return;
}


void ODE_system::setup_instance_varaibles(Param& param){

    //V_C.nT0, mw990b1a98_ff85_43ca_8bde_a72464820812, index: 0
    //Unit: mole^(1)
    _species_var[0] = PFILE(15) * 1.66053872801495e-24;
    //V_C.T0, mwbe4c5d3b_8194_4fa2_a6a0_b206aee50ca8, index: 1
    //Unit: mole^(1)
    _species_var[1] = PFILE(16) * 1.66053872801495e-24;
    //V_C.nT1, mw4b248b6e_21bd_47c2_bcff_8ee7efb6466d, index: 2
    //Unit: mole^(1)
    _species_var[2] = PFILE(17) * 1.66053872801495e-24;
    //V_C.T1, mwec02e6cb_81df_48b7_9ba6_e3811850cc36, index: 3
    //Unit: mole^(1)
    _species_var[3] = PFILE(18) * 1.66053872801495e-24;
    //V_C.aPD1, mw599ee6e0_5e8b_48e0_85ba_da151dee0b44, index: 4
    //Unit: mole^(1)metre^(-3)
    _species_var[4] = PFILE(19) * 1.0000000000000002e-06;
    //V_C.aPDL1, mwe6b103e9_401f_4766_a3e4_9e2441b8976b, index: 5
    //Unit: mole^(1)metre^(-3)
    _species_var[5] = PFILE(20) * 1.0000000000000002e-06;
    //V_C.aCTLA4, mw29d146a7_81c3_4114_bfcf_25cf9bcbdc09, index: 6
    //Unit: mole^(1)metre^(-3)
    _species_var[6] = PFILE(21) * 1.0000000000000002e-06;
    //V_C.Th, mw54775622_03d2_4992_bbb6_37bbe035c2af, index: 7
    //Unit: mole^(1)
    _species_var[7] = PFILE(22) * 1.66053872801495e-24;
    //V_C.cabozantinib, mw8e960f92_3a8a_4334_b8f6_b67b3ab9fc5b, index: 8
    //Unit: mole^(1)metre^(-3)
    _species_var[8] = PFILE(23) * 1.0000000000000002e-06;
    //V_C.A_site1, mw6c9f8cfa_a6a4_4c2b_9826_8fc7826a1746, index: 9
    //Unit: mole^(1)metre^(-3)
    _species_var[9] = PFILE(24) * 1.0000000000000002e-06;
    //V_C.A_site2, mwd03a4cc8_9c0a_486f_93b2_d4ac04d97db6, index: 10
    //Unit: mole^(1)metre^(-3)
    _species_var[10] = PFILE(25) * 1.0000000000000002e-06;
    //V_P.nT0, mw8ff04f75_bfdc_4e9b_b0c1_853c7b92d240, index: 11
    //Unit: mole^(1)
    _species_var[11] = PFILE(26) * 1.66053872801495e-24;
    //V_P.T0, mw446f2088_0470_49b0_93fd_dbb280284c6b, index: 12
    //Unit: mole^(1)
    _species_var[12] = PFILE(27) * 1.66053872801495e-24;
    //V_P.nT1, mwd9b14030_9855_4784_a319_c35483f01010, index: 13
    //Unit: mole^(1)
    _species_var[13] = PFILE(28) * 1.66053872801495e-24;
    //V_P.T1, mwe7795e50_c215_40e9_b89b_a8af06aeb52e, index: 14
    //Unit: mole^(1)
    _species_var[14] = PFILE(29) * 1.66053872801495e-24;
    //V_P.aPD1, mw02a3d31f_5669_4e64_90cc_3be8aa2ed3b2, index: 15
    //Unit: mole^(1)metre^(-3)
    _species_var[15] = PFILE(30) * 1.0000000000000002e-06;
    //V_P.aPDL1, mw7399b483_69eb_4c81_9f41_60682a3bfe21, index: 16
    //Unit: mole^(1)metre^(-3)
    _species_var[16] = PFILE(31) * 1.0000000000000002e-06;
    //V_P.aCTLA4, mwb44350c6_7e40_459c_90ae_bea9822bf60e, index: 17
    //Unit: mole^(1)metre^(-3)
    _species_var[17] = PFILE(32) * 1.0000000000000002e-06;
    //V_P.Th, mw8e792130_345b_49e6_9baa_9a85bbcafafc, index: 18
    //Unit: mole^(1)
    _species_var[18] = PFILE(33) * 1.66053872801495e-24;
    //V_P.cabozantinib, mw7628cb05_4b44_497b_bfc1_97ff16470bda, index: 19
    //Unit: mole^(1)metre^(-3)
    _species_var[19] = PFILE(34) * 1.0000000000000002e-06;
    //V_T.C_x, mw695eb91b_0911_4b47_87cc_8ba71a0bb820, index: 20
    //Unit: mole^(1)
    _species_var[20] = PFILE(35) * 1.66053872801495e-24;
    //V_T.T1_exh, mw55c274cd_d175_4d06_a3c6_efc583b199e3, index: 21
    //Unit: mole^(1)
    _species_var[21] = PFILE(36) * 1.66053872801495e-24;
    //V_T.Th_exh, mwf215a94c_08a8_4769_8445_c1ed1d71337f, index: 22
    //Unit: mole^(1)
    _species_var[22] = PFILE(37) * 1.66053872801495e-24;
    //V_T.C1, mwf899106c_3e91_4c2f_a4b3_655b9b697e18, index: 23
    //Unit: mole^(1)
    _species_var[23] = PFILE(38) * 1.66053872801495e-24;
    //V_T.K, mwc62b3ce9_b2c2_46b8_a6f0_6035bedf2552, index: 24
    //Unit: mole^(1)
    _species_var[24] = PFILE(39) * 1.66053872801495e-24;
    //V_T.c_vas, mwd782a63e_61fc_4d79_93a2_6c9297f0e219, index: 25
    //Unit: kilogram^(1)metre^(-3)
    _species_var[25] = PFILE(40) * 1e-09;
    //V_T.C2, mw3a41c2e9_391e_4655_a5c4_c334f8153e9d, index: 26
    //Unit: mole^(1)
    _species_var[26] = PFILE(41) * 1.66053872801495e-24;
    //V_T.T0, mwbe41d0f3_b878_4cde_aeb1_572adddfabd9, index: 27
    //Unit: mole^(1)
    _species_var[27] = PFILE(42) * 1.66053872801495e-24;
    //V_T.T1, mwf0e5bcf5_4425_40a0_a101_f06ce25dce06, index: 28
    //Unit: mole^(1)
    _species_var[28] = PFILE(43) * 1.66053872801495e-24;
    //V_T.IFNg, mwab7f54d6_776c_4130_8fa2_17f63a7fd139, index: 29
    //Unit: mole^(1)metre^(-3)
    _species_var[29] = PFILE(44) * 1e-06;
    //V_T.APC, mw3fe06188_fa1b_4478_91b8_89abae923b83, index: 30
    //Unit: mole^(1)
    _species_var[30] = PFILE(45) * 1.66053872801495e-24;
    //V_T.mAPC, mwd88254ee_7e6c_4da5_8438_5045eb70d862, index: 31
    //Unit: mole^(1)
    _species_var[31] = PFILE(46) * 1.66053872801495e-24;
    //V_T.P0, mwab0c6c04_55e8_4031_b96e_9269f3d305ae, index: 32
    //Unit: mole^(1)metre^(-3)
    _species_var[32] = PFILE(47) * 1000.0;
    //V_T.P1, mw09f3e0e9_76dc_49b6_a4ac_1b6cedca9d50, index: 33
    //Unit: mole^(1)metre^(-3)
    _species_var[33] = PFILE(48) * 1000.0;
    //V_T.aPD1, mw3b639500_a4fd_499e_8149_c35f297ecc9b, index: 34
    //Unit: mole^(1)metre^(-3)
    _species_var[34] = PFILE(49) * 1e-06;
    //V_T.aPDL1, mw9d536365_fda9_4c2f_9699_f708f43d955a, index: 35
    //Unit: mole^(1)metre^(-3)
    _species_var[35] = PFILE(50) * 1e-06;
    //V_T.aCTLA4, mwfa43acfd_66da_4452_99f7_5949011f60cd, index: 36
    //Unit: mole^(1)metre^(-3)
    _species_var[36] = PFILE(51) * 1e-06;
    //V_T.Th, mw971dfc71_6580_45df_9115_2770b9b21ebb, index: 37
    //Unit: mole^(1)
    _species_var[37] = PFILE(52) * 1.66053872801495e-24;
    //V_T.TGFb, mw2892056a_adb1_481e_8c7e_9a8d3f2f6730, index: 38
    //Unit: mole^(1)metre^(-3)
    _species_var[38] = PFILE(53) * 1e-06;
    //V_T.MDSC, mw41b70159_ec33_4917_9561_ba5ea1c408d7, index: 39
    //Unit: mole^(1)
    _species_var[39] = PFILE(54) * 1.66053872801495e-24;
    //V_T.NO, mw1fbdd8b5_4cfc_4525_bdb2_8462e62ad940, index: 40
    //Unit: mole^(1)metre^(-3)
    _species_var[40] = PFILE(55) * 1e-06;
    //V_T.ArgI, mwe75acd14_780e_46cb_b03b_566c62e54368, index: 41
    //Unit: mole^(1)metre^(-3)
    _species_var[41] = PFILE(56) * 1e-06;
    //V_T.CCL2, mw139135ac_0efa_46e0_a7a3_29bbd82d0924, index: 42
    //Unit: mole^(1)metre^(-3)
    _species_var[42] = PFILE(57) * 1e-06;
    //V_T.cabozantinib, mw2f1f6254_2774_4a30_9717_0056e9f235ac, index: 43
    //Unit: mole^(1)metre^(-3)
    _species_var[43] = PFILE(58) * 1e-06;
    //V_T.Mac_M1, mwa8ae5143_2284_4149_a6dc_4dcd2d9457c1, index: 44
    //Unit: mole^(1)
    _species_var[44] = PFILE(59) * 1.66053872801495e-24;
    //V_T.Mac_M2, mw89614aa3_c485_4488_8a41_cde493e86019, index: 45
    //Unit: mole^(1)
    _species_var[45] = PFILE(60) * 1.66053872801495e-24;
    //V_T.IL12, mw16e02140_05a6_4836_a2da_31e70d4a4b98, index: 46
    //Unit: mole^(1)metre^(-3)
    _species_var[46] = PFILE(61) * 1e-06;
    //V_T.IL10, mw650ad75e_c9f7_45ff_83e0_9d3e5eac3adc, index: 47
    //Unit: mole^(1)metre^(-3)
    _species_var[47] = PFILE(62) * 1e-06;
    //V_T.Fib, mw7d80bee7_6394_4935_90c3_79c7a52de3b9, index: 48
    //Unit: mole^(1)
    _species_var[48] = PFILE(63) * 1.66053872801495e-24;
    //V_T.CAF, mw67a0c834_e9c7_49d7_87d5_f26e954e33f8, index: 49
    //Unit: mole^(1)
    _species_var[49] = PFILE(64) * 1.66053872801495e-24;
    //V_T.ECM, mwb457314d_897b_4dbb_8136_ea3bd4d83674, index: 50
    //Unit: mole^(1)
    _species_var[50] = PFILE(65) * 1e-09;
    //V_LN.nT0, mw1346e2df_84a0_452b_876b_fbd03b119264, index: 51
    //Unit: mole^(1)
    _species_var[51] = PFILE(66) * 1.66053872801495e-24;
    //V_LN.aT0, mwe1ffcb62_da78_4c05_97ff_db7d696e5baa, index: 52
    //Unit: mole^(1)
    _species_var[52] = PFILE(67) * 1.66053872801495e-24;
    //V_LN.T0, mwc740cedc_cd5e_47a1_9bb5_3b000a801389, index: 53
    //Unit: mole^(1)
    _species_var[53] = PFILE(68) * 1.66053872801495e-24;
    //V_LN.IL2, mwdc7b4228_a75f_4b60_8666_e570ea862287, index: 54
    //Unit: mole^(1)metre^(-3)
    _species_var[54] = PFILE(69) * 1e-06;
    //V_LN.nT1, mw07efdfe8_0d7c_4307_baab_6db9445ded3b, index: 55
    //Unit: mole^(1)
    _species_var[55] = PFILE(70) * 1.66053872801495e-24;
    //V_LN.aT1, mw04938e68_9c84_4d97_9d35_b561388eddb9, index: 56
    //Unit: mole^(1)
    _species_var[56] = PFILE(71) * 1.66053872801495e-24;
    //V_LN.T1, mw3549b82b_8481_4ad0_b80e_cebb2fe29243, index: 57
    //Unit: mole^(1)
    _species_var[57] = PFILE(72) * 1.66053872801495e-24;
    //V_LN.APC, mw5c226c59_4710_418e_aa87_95fc6ca6c4a1, index: 58
    //Unit: mole^(1)
    _species_var[58] = PFILE(73) * 1.66053872801495e-24;
    //V_LN.mAPC, mwc10261b7_c6a8_4a06_a3b3_14ac19fe865b, index: 59
    //Unit: mole^(1)
    _species_var[59] = PFILE(74) * 1.66053872801495e-24;
    //V_LN.aPD1, mwb4995a34_fea9_4f45_82c6_16126b3d15a0, index: 60
    //Unit: mole^(1)metre^(-3)
    _species_var[60] = PFILE(75) * 1e-06;
    //V_LN.aPDL1, mw743ece23_2f72_493f_805a_8898121158b5, index: 61
    //Unit: mole^(1)metre^(-3)
    _species_var[61] = PFILE(76) * 1e-06;
    //V_LN.aCTLA4, mw78659cc4_fda8_4df3_8f85_657dde9f0fc4, index: 62
    //Unit: mole^(1)metre^(-3)
    _species_var[62] = PFILE(77) * 1e-06;
    //V_LN.aTh, mwe3461f35_2487_41fa_b067_e4de2e668a1b, index: 63
    //Unit: mole^(1)
    _species_var[63] = PFILE(78) * 1.66053872801495e-24;
    //V_LN.Th, mw7cbb65f2_9515_4f73_889b_935ee85e541a, index: 64
    //Unit: mole^(1)
    _species_var[64] = PFILE(79) * 1.66053872801495e-24;
    //V_LN.cabozantinib, mw0ad48b87_c0ea_4531_b20b_24d0c88396bb, index: 65
    //Unit: mole^(1)metre^(-3)
    _species_var[65] = PFILE(80) * 1e-06;
    //V_e.P0, mwd2e89639_7aeb_4f23_beda_4a3ddad9cd6f, index: 66
    //Unit: mole^(1)metre^(-3)
    _species_var[66] = PFILE(81) * 999.9999999999999;
    //V_e.p0, mw142ce0fe_d846_446e_8aeb_059c1720c0d1, index: 67
    //Unit: mole^(1)metre^(-3)
    _species_var[67] = PFILE(82) * 999.9999999999999;
    //V_e.P1, mwb39258b8_8b3d_42b9_8f96_7e482f00fc1a, index: 68
    //Unit: mole^(1)metre^(-3)
    _species_var[68] = PFILE(83) * 999.9999999999999;
    //V_e.p1, mwa0567cba_3831_42d4_8f54_53e239500044, index: 69
    //Unit: mole^(1)metre^(-3)
    _species_var[69] = PFILE(84) * 999.9999999999999;
    //A_e.M1, mwec98f8d0_1d60_4890_8c3e_a887d0c1664b, index: 70
    //Unit: mole^(1)metre^(-2)
    _species_var[70] = PFILE(85) * 1.66053872801495e-12;
    //A_e.M1p0, mw470858e9_08e9_416a_b3b4_9cb5cc02db6f, index: 71
    //Unit: mole^(1)metre^(-2)
    _species_var[71] = PFILE(86) * 1.66053872801495e-12;
    //A_e.M1p1, mw775901d6_6e45_410e_9b1b_142a8dc7fb4c, index: 72
    //Unit: mole^(1)metre^(-2)
    _species_var[72] = PFILE(87) * 1.66053872801495e-12;
    //A_s.M1, mw0decdeae_b572_43a7_b5ad_58cfed8bfe44, index: 73
    //Unit: mole^(1)metre^(-2)
    _species_var[73] = PFILE(88) * 1.66053872801495e-12;
    //A_s.M1p0, mw1c21f02d_8d4c_4e9b_b406_bee21183ceb7, index: 74
    //Unit: mole^(1)metre^(-2)
    _species_var[74] = PFILE(89) * 1.66053872801495e-12;
    //A_s.M1p1, mw904b8884_aa26_49be_aa30_ffc7400e8297, index: 75
    //Unit: mole^(1)metre^(-2)
    _species_var[75] = PFILE(90) * 1.66053872801495e-12;
    //syn_T_C1.PD1_PDL1, mwae73c6bf_8b5d_4dc6_ad5c_56f6659885e4, index: 76
    //Unit: mole^(1)metre^(-2)
    _species_var[76] = PFILE(91) * 1.66053872801495e-12;
    //syn_T_C1.PD1_PDL2, mwbfc2bbee_6128_4870_9971_173760b55a77, index: 77
    //Unit: mole^(1)metre^(-2)
    _species_var[77] = PFILE(92) * 1.66053872801495e-12;
    //syn_T_C1.PD1, mw9ca340ff_e790_447b_bbfc_fd3b49fca939, index: 78
    //Unit: mole^(1)metre^(-2)
    _species_var[78] = PFILE(93) * 1.66053872801495e-12;
    //syn_T_C1.PDL1, mwe8a14fa4_4c4b_44d6_a45d_6cd653378e12, index: 79
    //Unit: mole^(1)metre^(-2)
    _species_var[79] = PFILE(94) * 1.66053872801495e-12;
    //syn_T_C1.PDL2, mw93c88f0a_5d5f_4ca8_92d8_a5691fd78397, index: 80
    //Unit: mole^(1)metre^(-2)
    _species_var[80] = PFILE(95) * 1.66053872801495e-12;
    //syn_T_C1.PD1_aPD1, mw16880c33_c3d5_468e_b24e_b69c9e856dd2, index: 81
    //Unit: mole^(1)metre^(-2)
    _species_var[81] = PFILE(96) * 1.66053872801495e-12;
    //syn_T_C1.PD1_aPD1_PD1, mw0d013ba1_5ba8_4c2e_bb98_66141556094b, index: 82
    //Unit: mole^(1)metre^(-2)
    _species_var[82] = PFILE(97) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_aPDL1, mwb5af13a5_63e0_415e_846f_b35ad6ede042, index: 83
    //Unit: mole^(1)metre^(-2)
    _species_var[83] = PFILE(98) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_aPDL1_PDL1, mw7ceb7ab5_5996_42b3_92df_9cf4d3d48fcf, index: 84
    //Unit: mole^(1)metre^(-2)
    _species_var[84] = PFILE(99) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1, mwebe59e39_9534_4ea1_a8c2_3f89fd0c5e1b, index: 85
    //Unit: mole^(1)metre^(-2)
    _species_var[85] = PFILE(100) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_aPDL1, mwbeb18476_5ed3_42ef_8e96_26e389039c52, index: 86
    //Unit: mole^(1)metre^(-2)
    _species_var[86] = PFILE(101) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_aPDL1_TPDL1, mw31453900_5b91_43a2_a9b5_35c97aead58a, index: 87
    //Unit: mole^(1)metre^(-2)
    _species_var[87] = PFILE(102) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD80, mwf01692b1_326b_4a22_bd4c_b755593402dc, index: 88
    //Unit: mole^(1)metre^(-2)
    _species_var[88] = PFILE(103) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD80_CD28, mwcba1327f_672e_469b_b5d7_cd36e9ed2ad6, index: 89
    //Unit: mole^(1)metre^(-2)
    _species_var[89] = PFILE(104) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD86, mw6e878748_9dfb_440b_a108_f3fdb43b5497, index: 90
    //Unit: mole^(1)metre^(-2)
    _species_var[90] = PFILE(105) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4, mwed5bb7eb_2a28_4b2d_a417_78b5935eda65, index: 91
    //Unit: mole^(1)metre^(-2)
    _species_var[91] = PFILE(106) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4_CD80, mw12e9b975_4664_401f_8a11_dc77a36a5e7f, index: 92
    //Unit: mole^(1)metre^(-2)
    _species_var[92] = PFILE(107) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_CD80_CTLA4, mwec0de2f7_cf44_491f_8794_22e59e071449, index: 93
    //Unit: mole^(1)metre^(-2)
    _species_var[93] = PFILE(108) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4_CD80_CTLA4, mwe613d3cd_8f48_4429_93d8_2373e8d21001, index: 94
    //Unit: mole^(1)metre^(-2)
    _species_var[94] = PFILE(109) * 1.66053872801495e-12;
    //syn_T_C1.CD86_CTLA4, mw9420ee0f_f375_4112_b786_57a3aef48c7c, index: 95
    //Unit: mole^(1)metre^(-2)
    _species_var[95] = PFILE(110) * 1.66053872801495e-12;
    //syn_T_C1.CD86_CTLA4_CD86, mw5fde2bbf_45ff_44dc_8168_d822071f71d7, index: 96
    //Unit: mole^(1)metre^(-2)
    _species_var[96] = PFILE(111) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_CD80, mwcc79503f_6a89_4249_b51e_d51609124bd1, index: 97
    //Unit: mole^(1)metre^(-2)
    _species_var[97] = PFILE(112) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_CD80_CD28, mwf94367bc_d601_4e78_8260_35e28ae87986, index: 98
    //Unit: mole^(1)metre^(-2)
    _species_var[98] = PFILE(113) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_CD80_CTLA4, mw91081782_089d_4cb1_ad05_4c09be9ef29e, index: 99
    //Unit: mole^(1)metre^(-2)
    _species_var[99] = PFILE(114) * 1.66053872801495e-12;
    //syn_T_C1.CD28, mwf4891490_c093_4d8f_8935_3258ed0f59a9, index: 100
    //Unit: mole^(1)metre^(-2)
    _species_var[100] = PFILE(115) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4, mw131d2c1b_32db_4771_9572_4e707ffb9f54, index: 101
    //Unit: mole^(1)metre^(-2)
    _species_var[101] = PFILE(116) * 1.66053872801495e-12;
    //syn_T_C1.CD80, mwfeb58098_b313_498f_87a6_85f1ad2ede9f, index: 102
    //Unit: mole^(1)metre^(-2)
    _species_var[102] = PFILE(117) * 1.66053872801495e-12;
    //syn_T_C1.CD80m, mw3adedb40_1cbd_44b9_a96c_d84c0e22be9c, index: 103
    //Unit: mole^(1)metre^(-2)
    _species_var[103] = PFILE(118) * 1.66053872801495e-12;
    //syn_T_C1.CD86, mw76b56ce8_1dc1_4634_a9ec_ed271f8f5d41, index: 104
    //Unit: mole^(1)metre^(-2)
    _species_var[104] = PFILE(119) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_aCTLA4, mwef99a442_d956_4a16_a953_17ef4afc23c6, index: 105
    //Unit: mole^(1)metre^(-2)
    _species_var[105] = PFILE(120) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_aCTLA4_CTLA4, mw38fb6d74_2f5e_460f_8ad9_03f04ed4aa57, index: 106
    //Unit: mole^(1)metre^(-2)
    _species_var[106] = PFILE(121) * 1.66053872801495e-12;
    //syn_T_APC.PD1_PDL1, mwd1c54167_7ecd_4877_a122_c722cdcc843d, index: 107
    //Unit: mole^(1)metre^(-2)
    _species_var[107] = PFILE(122) * 1.66053872801495e-12;
    //syn_T_APC.PD1_PDL2, mw1465f9e1_86d8_422b_accb_925530dc910f, index: 108
    //Unit: mole^(1)metre^(-2)
    _species_var[108] = PFILE(123) * 1.66053872801495e-12;
    //syn_T_APC.PD1, mwaeeb32c5_065f_4ed9_8784_c564d3959936, index: 109
    //Unit: mole^(1)metre^(-2)
    _species_var[109] = PFILE(124) * 1.66053872801495e-12;
    //syn_T_APC.PDL1, mw03e869c0_7db6_4ac7_9d87_07c8f731d967, index: 110
    //Unit: mole^(1)metre^(-2)
    _species_var[110] = PFILE(125) * 1.66053872801495e-12;
    //syn_T_APC.PDL2, mw81cb45c6_9248_4f01_b4a6_0de3b3482138, index: 111
    //Unit: mole^(1)metre^(-2)
    _species_var[111] = PFILE(126) * 1.66053872801495e-12;
    //syn_T_APC.PD1_aPD1, mw69bb9f05_256e_4b31_a4af_ed2feb9cec17, index: 112
    //Unit: mole^(1)metre^(-2)
    _species_var[112] = PFILE(127) * 1.66053872801495e-12;
    //syn_T_APC.PD1_aPD1_PD1, mw1040aff1_be4e_4118_a1f9_457f9842d74b, index: 113
    //Unit: mole^(1)metre^(-2)
    _species_var[113] = PFILE(128) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_aPDL1, mwbafb63d5_f38a_4aea_bab8_f669f06c8510, index: 114
    //Unit: mole^(1)metre^(-2)
    _species_var[114] = PFILE(129) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_aPDL1_PDL1, mw713d89d3_add9_4e69_9ca1_4739e759b573, index: 115
    //Unit: mole^(1)metre^(-2)
    _species_var[115] = PFILE(130) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1, mwbe9d09b9_e7e0_4ad8_b455_da59484143c6, index: 116
    //Unit: mole^(1)metre^(-2)
    _species_var[116] = PFILE(131) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_aPDL1, mwcaea15a3_51ae_4ad9_b6f4_917db3584d91, index: 117
    //Unit: mole^(1)metre^(-2)
    _species_var[117] = PFILE(132) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_aPDL1_TPDL1, mw47027b6a_2b5f_4644_9344_ea9defebbf29, index: 118
    //Unit: mole^(1)metre^(-2)
    _species_var[118] = PFILE(133) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD80, mw70a59707_8442_415d_a202_546c78d1a920, index: 119
    //Unit: mole^(1)metre^(-2)
    _species_var[119] = PFILE(134) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD80_CD28, mw199aa2fb_153a_4f2e_adf9_7eefc2de0322, index: 120
    //Unit: mole^(1)metre^(-2)
    _species_var[120] = PFILE(135) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD86, mw1b53f0d5_69dd_4a8d_a033_88b161fef14f, index: 121
    //Unit: mole^(1)metre^(-2)
    _species_var[121] = PFILE(136) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4, mwefaf349f_b97b_45be_b91f_59e750361eb5, index: 122
    //Unit: mole^(1)metre^(-2)
    _species_var[122] = PFILE(137) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4_CD80, mw4d6d7229_8847_4515_bb79_3e9f548ae5d6, index: 123
    //Unit: mole^(1)metre^(-2)
    _species_var[123] = PFILE(138) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_CD80_CTLA4, mwf4538f4d_a8e7_49a5_b8aa_6847f8c714e7, index: 124
    //Unit: mole^(1)metre^(-2)
    _species_var[124] = PFILE(139) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4_CD80_CTLA4, mw0cddfb1c_8690_43a2_8fa6_80fdaabc6884, index: 125
    //Unit: mole^(1)metre^(-2)
    _species_var[125] = PFILE(140) * 1.66053872801495e-12;
    //syn_T_APC.CD86_CTLA4, mwcc31a173_7c83_4d7a_8fe5_3b1dfb12f494, index: 126
    //Unit: mole^(1)metre^(-2)
    _species_var[126] = PFILE(141) * 1.66053872801495e-12;
    //syn_T_APC.CD86_CTLA4_CD86, mw15405340_f645_4e49_a72a_34806d67707f, index: 127
    //Unit: mole^(1)metre^(-2)
    _species_var[127] = PFILE(142) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_CD80, mw7844bd33_4dc9_433a_a342_518ed7ed2beb, index: 128
    //Unit: mole^(1)metre^(-2)
    _species_var[128] = PFILE(143) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_CD80_CD28, mw7eab10cc_30fa_4944_8a8f_907c4cf76603, index: 129
    //Unit: mole^(1)metre^(-2)
    _species_var[129] = PFILE(144) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_CD80_CTLA4, mw58e24515_30c4_464e_a89d_ac474362aaef, index: 130
    //Unit: mole^(1)metre^(-2)
    _species_var[130] = PFILE(145) * 1.66053872801495e-12;
    //syn_T_APC.CD28, mwd27906bd_520d_4a1d_8543_072396b98fa0, index: 131
    //Unit: mole^(1)metre^(-2)
    _species_var[131] = PFILE(146) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4, mw4fb7ddf5_f2fd_4c8c_ab95_0ab92e591776, index: 132
    //Unit: mole^(1)metre^(-2)
    _species_var[132] = PFILE(147) * 1.66053872801495e-12;
    //syn_T_APC.CD80, mwbdd65cc2_1d51_41fa_b407_6e3687e6dfe0, index: 133
    //Unit: mole^(1)metre^(-2)
    _species_var[133] = PFILE(148) * 1.66053872801495e-12;
    //syn_T_APC.CD80m, mwc0c73656_e325_4afe_86a1_829d360aa39a, index: 134
    //Unit: mole^(1)metre^(-2)
    _species_var[134] = PFILE(149) * 1.66053872801495e-12;
    //syn_T_APC.CD86, mwb76d7c96_3a66_427f_97ba_7e56d8d2ffdb, index: 135
    //Unit: mole^(1)metre^(-2)
    _species_var[135] = PFILE(150) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_aCTLA4, mw0ced11f0_fa2f_43cf_9ad8_cd34cb199608, index: 136
    //Unit: mole^(1)metre^(-2)
    _species_var[136] = PFILE(151) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_aCTLA4_CTLA4, mwe0088f76_aecf_4d6d_a474_170d0f075bea, index: 137
    //Unit: mole^(1)metre^(-2)
    _species_var[137] = PFILE(152) * 1.66053872801495e-12;
    //syn_M_C.CD47, mw1d2d953b_26ff_4c7c_96cd_3a2b9a5c13e6, index: 138
    //Unit: mole^(1)metre^(-2)
    _species_var[138] = PFILE(153) * 1.66053872801495e-12;
    //syn_M_C.SIRPa, mw7b657756_85c0_41fe_93c5_60608619b892, index: 139
    //Unit: mole^(1)metre^(-2)
    _species_var[139] = PFILE(154) * 1.66053872801495e-12;
    //syn_M_C.CD47_SIRPa, mwb868c6f7_885d_4570_94be_3c2aa8744849, index: 140
    //Unit: mole^(1)metre^(-2)
    _species_var[140] = PFILE(155) * 1.66053872801495e-12;
    //syn_M_C.PD1_PDL1, mw8ad5d321_e7f7_431f_8b54_bd3ea2c3b10b, index: 141
    //Unit: mole^(1)metre^(-2)
    _species_var[141] = PFILE(158) * 1.66053872801495e-12;
    //syn_M_C.PD1_PDL2, mwe8258224_2bdc_408a_bdf0_d5895dce6516, index: 142
    //Unit: mole^(1)metre^(-2)
    _species_var[142] = PFILE(159) * 1.66053872801495e-12;
    //syn_M_C.PD1, mwd3571b1f_5d64_43c2_908d_0140c8efc50d, index: 143
    //Unit: mole^(1)metre^(-2)
    _species_var[143] = PFILE(160) * 1.66053872801495e-12;
    //syn_M_C.PDL1, mwc64be3bf_05d1_4154_aad4_7fd95135bc9d, index: 144
    //Unit: mole^(1)metre^(-2)
    _species_var[144] = PFILE(161) * 1.66053872801495e-12;
    //syn_M_C.PDL2, mwee311bbc_1b0c_4016_8adb_620ce18cddcc, index: 145
    //Unit: mole^(1)metre^(-2)
    _species_var[145] = PFILE(162) * 1.66053872801495e-12;
    //syn_M_C.PD1_aPD1, mw4e4530a6_ea91_4c0c_9c39_576efb8bb6ba, index: 146
    //Unit: mole^(1)metre^(-2)
    _species_var[146] = PFILE(163) * 1.66053872801495e-12;
    //syn_M_C.PD1_aPD1_PD1, mw234985da_cad0_4f3a_a886_bd474dfe4693, index: 147
    //Unit: mole^(1)metre^(-2)
    _species_var[147] = PFILE(164) * 1.66053872801495e-12;
    //syn_M_C.PDL1_aPDL1, mw88db6728_e33c_4d4a_b59b_8684b6cb59d4, index: 148
    //Unit: mole^(1)metre^(-2)
    _species_var[148] = PFILE(165) * 1.66053872801495e-12;
    //syn_M_C.PDL1_aPDL1_PDL1, mw948f4fbf_6aa8_4a4b_8c8d_e9c4142ef79e, index: 149
    //Unit: mole^(1)metre^(-2)
    _species_var[149] = PFILE(166) * 1.66053872801495e-12;
    //syn_M_C.PDL1_CD80, mw2da3fcab_64c9_44fb_940f_70c390049c6b, index: 150
    //Unit: mole^(1)metre^(-2)
    _species_var[150] = PFILE(167) * 1.66053872801495e-12;
    //syn_M_C.CD80, mw16edf955_3075_447a_86ce_4edede243ab3, index: 151
    //Unit: mole^(1)metre^(-2)
    _species_var[151] = PFILE(168) * 1.66053872801495e-12;
    //syn_M_C.CD80m, mw6d7c01df_7cda_4d4d_bc0f_8cd42486753e, index: 152
    //Unit: mole^(1)metre^(-2)
    _species_var[152] = PFILE(169) * 1.66053872801495e-12;
    
    return;
}
    

void ODE_system::adjust_hybrid_variables(void){

}    

void ODE_system::setup_instance_tolerance(Param& param){

    //Tolerance
    realtype reltol = PFILE(3);
    realtype abstol_base = PFILE(4);
    N_Vector abstol = N_VNew_Serial(_neq);

    for (size_t i = 0; i < 153; i++)
    {
        NV_DATA_S(abstol)[i] = abstol_base * get_unit_conversion_species(i);
    }
    int flag = CVodeSVtolerances(_cvode_mem, reltol, abstol);
    check_flag(&flag, "CVodeSVtolerances", 1);

    
    return;
}

void ODE_system::eval_init_assignment(void){
    //Assignment Rules required before IA
    //InitialAssignment

    updateVar();
    
    return;
}
void ODE_system::setupEvents(void){

    _nevent = 3;
    _nroot = 3;

    _trigger_element_type = std::vector<EVENT_TRIGGER_ELEM_TYPE>(_nroot, TRIGGER_NON_INSTANT);
    _trigger_element_satisfied = std::vector<bool>(_nroot, false);
    _event_triggered = std::vector<bool>(_nevent, false);

    //C_total < (0.5 * cell)
    _trigger_element_type[0] = TRIGGER_NON_INSTANT;

    //V_T.C1 < (0.5 * cell)
    _trigger_element_type[1] = TRIGGER_NON_INSTANT;

    //V_T.C2 < (0.5 * cell)
    _trigger_element_type[2] = TRIGGER_NON_INSTANT;

    _event_triggered[0] = true;

    _event_triggered[1] = true;

    _event_triggered[2] = true;

    return;
}
int ODE_system::f(realtype t, N_Vector y, N_Vector ydot, void *user_data){

    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);

    //Assignment rules:

    realtype AUX_VAR_C_total = 0.0 * PARAM(10) + SPVAR(23) + SPVAR(26);

    realtype AUX_VAR_T_total = 0.0 * PARAM(10) + SPVAR(27) + SPVAR(28) + SPVAR(37);

    realtype AUX_VAR_M_total = SPVAR(44) + SPVAR(45);

    realtype AUX_VAR_T_total_LN = 0.0 * PARAM(10) + SPVAR(53) + SPVAR(57) + SPVAR(64);

    realtype AUX_VAR_k_C1_therapy = 0.0 / PARAM(11);

    realtype AUX_VAR_Tregs_ = SPVAR(27);

    realtype AUX_VAR_H_PD1_C1 = std::pow((SPVAR(76) + SPVAR(77)) / PARAM(175), PARAM(176)) / (std::pow((SPVAR(76) + SPVAR(77)) / PARAM(175), PARAM(176)) + 1.0);

    realtype AUX_VAR_H_TGFb_Teff = SPVAR(38) / (SPVAR(38) + PARAM(196));

    realtype AUX_VAR_k_C2_therapy = 0.0 / PARAM(11);

    realtype AUX_VAR_H_IL10_phago = SPVAR(47) / (PARAM(259) + SPVAR(47));

    realtype AUX_VAR_C_max = SPVAR(24);

    realtype AUX_VAR_H_CD28_APC = std::pow((SPVAR(119) + SPVAR(121) + 2.0 * SPVAR(120) + SPVAR(129)) / PARAM(177), PARAM(178)) / (std::pow((SPVAR(119) + SPVAR(121) + 2.0 * SPVAR(120) + SPVAR(129)) / PARAM(177), PARAM(178)) + 1.0);

    realtype AUX_VAR_H_APC = PARAM(81) * SPVAR(59) / (PARAM(81) * SPVAR(59) + SPVAR(51) * PARAM(26) + PARAM(10));

    realtype AUX_VAR_H_mAPC = PARAM(81) * SPVAR(59) / (PARAM(81) * SPVAR(59) + SPVAR(55) * PARAM(53) + PARAM(10));

    realtype AUX_VAR_H_APCh = PARAM(81) * SPVAR(59) / (PARAM(81) * SPVAR(59) + SPVAR(51) * PARAM(53) + PARAM(10));

    realtype AUX_VAR_pTCR_p0_MHC_tot = PARAM(98) / (PARAM(98) + PARAM(100)) * std::pow(PARAM(99) / (PARAM(98) + PARAM(99)), PARAM(101)) * 0.5 * (SPVAR(74) / PARAM(26) + PARAM(102) + PARAM(98) / PARAM(97) - PARAM(102) * std::sqrt(std::pow((SPVAR(74) / PARAM(26) + PARAM(102) + PARAM(98) / PARAM(97)) / PARAM(102), 2.0) - 4.0 * SPVAR(74) / PARAM(26) / PARAM(102)));

    realtype AUX_VAR_pTCR_p1_MHC_tot = PARAM(113) / (PARAM(113) + PARAM(115)) * std::pow(PARAM(114) / (PARAM(113) + PARAM(114)), PARAM(116)) * 0.5 * (SPVAR(75) / PARAM(53) + PARAM(117) + PARAM(113) / PARAM(112) - PARAM(117) * std::sqrt(std::pow((SPVAR(75) / PARAM(53) + PARAM(117) + PARAM(113) / PARAM(112)) / PARAM(117), 2.0) - 4.0 * SPVAR(75) / PARAM(53) / PARAM(117)));

    realtype AUX_VAR_syn_T_C1_PDL1_total = SPVAR(79) + SPVAR(76) + SPVAR(83) + 2.0 * SPVAR(84) + SPVAR(97) + SPVAR(98) + SPVAR(99);

    realtype AUX_VAR_syn_T_C1_PDL2_total = SPVAR(77) + SPVAR(80);

    realtype AUX_VAR_H_CD28_C1 = std::pow((SPVAR(88) + SPVAR(90) + 2.0 * SPVAR(89) + SPVAR(98)) / PARAM(177), PARAM(178)) / (std::pow((SPVAR(88) + SPVAR(90) + 2.0 * SPVAR(89) + SPVAR(98)) / PARAM(177), PARAM(178)) + 1.0);

    realtype AUX_VAR_syn_T_APC_PDL1_total = SPVAR(110) + SPVAR(107) + SPVAR(114) + 2.0 * SPVAR(115) + SPVAR(128) + SPVAR(129) + SPVAR(130);

    realtype AUX_VAR_syn_T_APC_PDL2_total = SPVAR(108) + SPVAR(111);

    realtype AUX_VAR_H_PD1_APC = std::pow((SPVAR(107) + SPVAR(108)) / PARAM(175), PARAM(176)) / (std::pow((SPVAR(107) + SPVAR(108)) / PARAM(175), PARAM(176)) + 1.0);

    realtype AUX_VAR_H_TGFb = SPVAR(38) / (SPVAR(38) + PARAM(195));

    realtype AUX_VAR_H_NO = SPVAR(40) / (PARAM(212) + SPVAR(40));

    realtype AUX_VAR_H_ArgI_Teff = SPVAR(41) / (PARAM(211) + SPVAR(41));

    realtype AUX_VAR_H_ArgI_Treg = SPVAR(41) / (PARAM(213) + SPVAR(41));

    realtype AUX_VAR_R_cabo = SPVAR(43) / (SPVAR(43) + PARAM(231));

    realtype AUX_VAR_H_therapy_cabo = SPVAR(43) / (SPVAR(43) + PARAM(232));

    realtype AUX_VAR_PDL1_total = SPVAR(144) + SPVAR(141) + SPVAR(148) + 2.0 * SPVAR(149) + SPVAR(150);

    realtype AUX_VAR_PDL2_total = SPVAR(142) + SPVAR(145);

    realtype AUX_VAR_H_SIRPa = std::pow(SPVAR(140) / PARAM(253), PARAM(254)) / (std::pow(SPVAR(140) / PARAM(253), PARAM(254)) + 1.0);

    realtype AUX_VAR_H_PD1_M = std::pow((SPVAR(141) + SPVAR(142)) / PARAM(175), PARAM(176)) / (std::pow((SPVAR(141) + SPVAR(142)) / PARAM(175), PARAM(176)) + 1.0);

    realtype AUX_VAR_H_IL10 = SPVAR(47) / (PARAM(246) + SPVAR(47));

    realtype AUX_VAR_H_IL12 = SPVAR(46) / (PARAM(247) + SPVAR(46));

    realtype AUX_VAR_V_T = ((SPVAR(20) + AUX_VAR_C_total) * PARAM(12) + (SPVAR(21) + SPVAR(22) + AUX_VAR_T_total) * PARAM(13)) / PARAM(14) + AUX_VAR_M_total * PARAM(250) / PARAM(14) + SPVAR(48) * PARAM(275) / PARAM(14) + SPVAR(49) * PARAM(276) / PARAM(14) + SPVAR(50) * PARAM(269) / PARAM(270);

    realtype AUX_VAR_N_aT = PARAM(47) + PARAM(48) * AUX_VAR_H_CD28_APC + PARAM(49) * SPVAR(54) / (PARAM(45) + SPVAR(54));

    realtype AUX_VAR_N_aT0 = PARAM(47) + PARAM(48) * AUX_VAR_H_CD28_APC + PARAM(50) * SPVAR(54) / (PARAM(45) + SPVAR(54));

    realtype AUX_VAR_N_aTh = PARAM(47) + PARAM(48) * AUX_VAR_H_CD28_APC + PARAM(50) * SPVAR(54) / (PARAM(45) + SPVAR(54));

    realtype AUX_VAR_H_P0 = AUX_VAR_pTCR_p0_MHC_tot / (AUX_VAR_pTCR_p0_MHC_tot + PARAM(90));

    realtype AUX_VAR_H_P1 = AUX_VAR_pTCR_p1_MHC_tot / (AUX_VAR_pTCR_p1_MHC_tot + PARAM(109));

    realtype AUX_VAR_H_MDSC = 1.0 - (1.0 - AUX_VAR_H_ArgI_Teff);

    realtype AUX_VAR_H_Mac_C = 1.0 - (1.0 - AUX_VAR_H_SIRPa) * (1.0 - AUX_VAR_H_PD1_M);

    realtype AUX_VAR_ECM_level = SPVAR(50) / AUX_VAR_V_T;

    realtype AUX_VAR_H_ECM_T_mot = AUX_VAR_ECM_level / (PARAM(274) + AUX_VAR_ECM_level);

    realtype AUX_VAR_H_ECM_T_exh = AUX_VAR_ECM_level / (PARAM(273) + AUX_VAR_ECM_level);

    realtype AUX_VAR_R_Tcell = 0.0 * PARAM(10) / PARAM(11) + (PARAM(16) + AUX_VAR_k_C1_therapy) * SPVAR(23) + PARAM(70) * SPVAR(23) * SPVAR(28) / (PARAM(73) * AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * SPVAR(28) / (SPVAR(28) + PARAM(74) * AUX_VAR_Tregs_ + PARAM(10)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC) * (1.0 - AUX_VAR_H_TGFb_Teff) * (1.0 - AUX_VAR_H_ECM_T_mot) + (PARAM(24) + AUX_VAR_k_C2_therapy) * SPVAR(26) + PARAM(70) * SPVAR(26) * SPVAR(28) / (PARAM(73) * AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * SPVAR(28) / (SPVAR(28) + PARAM(74) * AUX_VAR_Tregs_ + PARAM(10)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC) * (1.0 - AUX_VAR_H_TGFb_Teff) * (1.0 - AUX_VAR_H_ECM_T_mot) + PARAM(249) * SPVAR(23) * SPVAR(44) / (SPVAR(44) + PARAM(260) * AUX_VAR_C_total + PARAM(10)) * (1.0 - AUX_VAR_H_Mac_C) * (1.0 - AUX_VAR_H_IL10_phago) + PARAM(249) * SPVAR(26) * SPVAR(44) / (SPVAR(44) + PARAM(260) * AUX_VAR_C_total + PARAM(10)) * (1.0 - AUX_VAR_H_Mac_C) * (1.0 - AUX_VAR_H_IL10_phago);

    //Reaction fluxes:

    realtype ReactionFlux1 = PARAM(9) * SPVAR(20);

    realtype ReactionFlux2 = PARAM(9) * SPVAR(21);

    realtype ReactionFlux3 = PARAM(9) * SPVAR(22);

    realtype ReactionFlux4 = PARAM(15) * SPVAR(23) * (1.0 - AUX_VAR_C_total / AUX_VAR_C_max) * (1.0 - PARAM(234) * SPVAR(43) / (SPVAR(43) + PARAM(229) + PARAM(230)));

    realtype ReactionFlux5 = PARAM(16) * SPVAR(23);

    realtype ReactionFlux6 = PARAM(20) * AUX_VAR_C_total;

    realtype ReactionFlux7 = PARAM(21) * SPVAR(25) * AUX_VAR_V_T;

    realtype ReactionFlux8 = PARAM(18) * AUX_VAR_C_total * SPVAR(25) / (SPVAR(25) + PARAM(22)) * (1.0 - PARAM(234) * AUX_VAR_H_therapy_cabo);

    realtype ReactionFlux9 = PARAM(19) * SPVAR(24) * std::pow(std::pow(AUX_VAR_C_total / PARAM(10) * 2.5699999999999995e-06, 1.0 / 3.0), 2.0);

    realtype ReactionFlux10 = PARAM(23) * SPVAR(26) * (1.0 - AUX_VAR_C_total / AUX_VAR_C_max);

    realtype ReactionFlux11 = PARAM(24) * SPVAR(26);

    realtype ReactionFlux12 = PARAM(38) / PARAM(25);

    realtype ReactionFlux13 = PARAM(39) / PARAM(25) * SPVAR(11) / (PARAM(40) / PARAM(25) + SPVAR(11));

    realtype ReactionFlux14 = PARAM(39) / PARAM(25) * SPVAR(51) / (PARAM(40) / PARAM(25) + SPVAR(51));

    realtype ReactionFlux15 = PARAM(41) * SPVAR(11);

    realtype ReactionFlux16 = PARAM(41) * SPVAR(0);

    realtype ReactionFlux17 = PARAM(41) * SPVAR(51);

    realtype ReactionFlux18 = PARAM(36) * SPVAR(0);

    realtype ReactionFlux19 = PARAM(37) * SPVAR(11);

    realtype ReactionFlux20 = PARAM(27) * SPVAR(0);

    realtype ReactionFlux21 = PARAM(29) * SPVAR(51);

    realtype ReactionFlux22 = PARAM(30) * AUX_VAR_H_APC * AUX_VAR_H_P0 * SPVAR(51);

    realtype ReactionFlux23 = PARAM(30) * AUX_VAR_H_APC * AUX_VAR_H_P0 * SPVAR(51) * PARAM(26);

    realtype ReactionFlux24 = PARAM(31) / AUX_VAR_N_aT0 * SPVAR(52);

    realtype ReactionFlux25 = PARAM(31) / AUX_VAR_N_aT0 * std::pow(2.0, AUX_VAR_N_aT0) * SPVAR(52);

    realtype ReactionFlux26 = PARAM(32) * SPVAR(1);

    realtype ReactionFlux27 = PARAM(32) * SPVAR(12);

    realtype ReactionFlux28 = PARAM(32) * SPVAR(53);

    realtype ReactionFlux29 = PARAM(32) * SPVAR(27);

    realtype ReactionFlux30 = PARAM(9) * SPVAR(27) * (PARAM(197) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux31 = PARAM(33) * SPVAR(1);

    realtype ReactionFlux32 = PARAM(34) * SPVAR(12);

    realtype ReactionFlux33 = PARAM(35) * AUX_VAR_V_T * SPVAR(1) * (std::pow(AUX_VAR_C_total, 2.0) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux34 = PARAM(28) * SPVAR(53);

    realtype ReactionFlux35 = PARAM(42) * SPVAR(54) * PARAM(2);

    realtype ReactionFlux36 = PARAM(43) * SPVAR(57) * SPVAR(54) / (PARAM(45) + SPVAR(54));

    realtype ReactionFlux37 = PARAM(43) * SPVAR(53) * SPVAR(54) / (PARAM(46) + SPVAR(54));

    realtype ReactionFlux38 = PARAM(65) / PARAM(52);

    realtype ReactionFlux39 = PARAM(66) / PARAM(52) * SPVAR(13) / (PARAM(67) / PARAM(52) + SPVAR(13));

    realtype ReactionFlux40 = PARAM(66) / PARAM(52) * SPVAR(55) / (PARAM(67) / PARAM(52) + SPVAR(55));

    realtype ReactionFlux41 = PARAM(68) * SPVAR(13);

    realtype ReactionFlux42 = PARAM(68) * SPVAR(2);

    realtype ReactionFlux43 = PARAM(68) * SPVAR(55);

    realtype ReactionFlux44 = PARAM(63) * SPVAR(2);

    realtype ReactionFlux45 = PARAM(64) * SPVAR(13);

    realtype ReactionFlux46 = PARAM(54) * SPVAR(2);

    realtype ReactionFlux47 = PARAM(56) * SPVAR(55);

    realtype ReactionFlux48 = PARAM(57) * AUX_VAR_H_mAPC * AUX_VAR_H_P1 * SPVAR(55);

    realtype ReactionFlux49 = PARAM(57) * AUX_VAR_H_mAPC * AUX_VAR_H_P1 * SPVAR(55) * PARAM(53);

    realtype ReactionFlux50 = PARAM(58) / AUX_VAR_N_aT * SPVAR(56);

    realtype ReactionFlux51 = PARAM(58) / AUX_VAR_N_aT * std::pow(2.0, AUX_VAR_N_aT) * SPVAR(56);

    realtype ReactionFlux52 = PARAM(59) * SPVAR(3);

    realtype ReactionFlux53 = PARAM(59) * SPVAR(14);

    realtype ReactionFlux54 = PARAM(59) * SPVAR(57);

    realtype ReactionFlux55 = PARAM(59) * SPVAR(28);

    realtype ReactionFlux56 = PARAM(9) * SPVAR(28) * (PARAM(197) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux57 = PARAM(51) * SPVAR(28) * AUX_VAR_Tregs_ / (SPVAR(28) + AUX_VAR_Tregs_ + PARAM(10)) * AUX_VAR_H_IL10;

    realtype ReactionFlux58 = PARAM(69) * SPVAR(28) * AUX_VAR_C_total / (AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux59 = PARAM(71) * SPVAR(28) * AUX_VAR_H_ECM_T_exh;

    realtype ReactionFlux60 = PARAM(60) * SPVAR(3);

    realtype ReactionFlux61 = PARAM(61) * SPVAR(14);

    realtype ReactionFlux62 = PARAM(62) * AUX_VAR_V_T * SPVAR(3) * (std::pow(AUX_VAR_C_total, 2.0) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux63 = PARAM(55) * SPVAR(57);

    realtype ReactionFlux64 = PARAM(44) * SPVAR(56);

    realtype ReactionFlux65 = PARAM(72) * SPVAR(28);

    realtype ReactionFlux66 = PARAM(70) * SPVAR(23) * SPVAR(28) / (PARAM(73) * AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * SPVAR(28) / (SPVAR(28) + PARAM(74) * AUX_VAR_Tregs_ + PARAM(10)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_TGFb_Teff) * (1.0 - AUX_VAR_H_ECM_T_mot) * (1.0 - AUX_VAR_H_MDSC);

    realtype ReactionFlux67 = PARAM(70) * SPVAR(26) * SPVAR(28) / (PARAM(73) * AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * SPVAR(28) / (SPVAR(28) + PARAM(74) * AUX_VAR_Tregs_ + PARAM(10)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_TGFb_Teff) * (1.0 - AUX_VAR_H_ECM_T_mot) * (1.0 - AUX_VAR_H_MDSC);

    realtype ReactionFlux68 = PARAM(77) * (PARAM(79) * AUX_VAR_V_T - SPVAR(30));

    realtype ReactionFlux69 = PARAM(77) * (PARAM(80) * PARAM(2) - SPVAR(58));

    realtype ReactionFlux70 = PARAM(75) * SPVAR(30) * AUX_VAR_H_IL12 * (1.0 - AUX_VAR_H_IL10);

    realtype ReactionFlux71 = PARAM(76) * SPVAR(31);

    realtype ReactionFlux72 = PARAM(78) * SPVAR(31);

    realtype ReactionFlux73 = PARAM(78) * SPVAR(59);

    realtype ReactionFlux74 = PARAM(83) * SPVAR(70) * PARAM(4) - PARAM(82) * SPVAR(73) * PARAM(5);

    realtype ReactionFlux75 = PARAM(26) * (PARAM(91) * (PARAM(16) + AUX_VAR_k_C1_therapy) * SPVAR(23) + PARAM(91) * PARAM(70) * SPVAR(23) * SPVAR(28) / (PARAM(73) * AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * SPVAR(28) / (SPVAR(28) + PARAM(74) * AUX_VAR_Tregs_ + PARAM(10)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC) * (1.0 - AUX_VAR_H_TGFb_Teff) * (1.0 - AUX_VAR_H_ECM_T_mot) + PARAM(92) * (PARAM(24) + AUX_VAR_k_C2_therapy) * SPVAR(26) + PARAM(92) * PARAM(70) * SPVAR(26) * SPVAR(28) / (PARAM(73) * AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * SPVAR(28) / (SPVAR(28) + PARAM(74) * AUX_VAR_Tregs_ + PARAM(10)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC) * (1.0 - AUX_VAR_H_TGFb_Teff) * (1.0 - AUX_VAR_H_ECM_T_mot));

    realtype ReactionFlux76 = PARAM(85) * SPVAR(32) * AUX_VAR_V_T;

    realtype ReactionFlux77 = PARAM(84) * SPVAR(30) * SPVAR(32) * AUX_VAR_V_T;

    realtype ReactionFlux78 = PARAM(84) * PARAM(10) * SPVAR(32) * PARAM(3);

    realtype ReactionFlux79 = PARAM(86) * SPVAR(66) * PARAM(3);

    realtype ReactionFlux80 = PARAM(87) * SPVAR(67) * PARAM(3);

    realtype ReactionFlux81 = PARAM(88) * SPVAR(67) * SPVAR(70) * PARAM(4);

    realtype ReactionFlux82 = PARAM(89) * PARAM(88) * SPVAR(71) * PARAM(4);

    realtype ReactionFlux83 = PARAM(89) * PARAM(88) * SPVAR(74) * PARAM(5);

    realtype ReactionFlux84 = PARAM(83) * SPVAR(71) * PARAM(4);

    realtype ReactionFlux85 = PARAM(53) * (PARAM(110) * (PARAM(16) + AUX_VAR_k_C1_therapy) * SPVAR(23) + PARAM(110) * PARAM(70) * SPVAR(23) * SPVAR(28) / (PARAM(73) * AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * SPVAR(28) / (SPVAR(28) + PARAM(74) * AUX_VAR_Tregs_ + PARAM(10)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC) * (1.0 - AUX_VAR_H_TGFb_Teff) * (1.0 - AUX_VAR_H_ECM_T_mot) + PARAM(111) * (PARAM(24) + AUX_VAR_k_C2_therapy) * SPVAR(26) + PARAM(111) * PARAM(70) * SPVAR(26) * SPVAR(28) / (PARAM(73) * AUX_VAR_C_total + SPVAR(28) + PARAM(10)) * SPVAR(28) / (SPVAR(28) + PARAM(74) * AUX_VAR_Tregs_ + PARAM(10)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC) * (1.0 - AUX_VAR_H_TGFb_Teff) * (1.0 - AUX_VAR_H_ECM_T_mot));

    realtype ReactionFlux86 = PARAM(104) * SPVAR(33) * AUX_VAR_V_T;

    realtype ReactionFlux87 = PARAM(103) * SPVAR(30) * SPVAR(33) * AUX_VAR_V_T;

    realtype ReactionFlux88 = PARAM(103) * PARAM(10) * SPVAR(33) * PARAM(3);

    realtype ReactionFlux89 = PARAM(105) * SPVAR(68) * PARAM(3);

    realtype ReactionFlux90 = PARAM(106) * SPVAR(69) * PARAM(3);

    realtype ReactionFlux91 = PARAM(107) * SPVAR(69) * SPVAR(70) * PARAM(4);

    realtype ReactionFlux92 = PARAM(108) * PARAM(107) * SPVAR(72) * PARAM(4);

    realtype ReactionFlux93 = PARAM(108) * PARAM(107) * SPVAR(75) * PARAM(5);

    realtype ReactionFlux94 = PARAM(83) * SPVAR(72) * PARAM(4);

    realtype ReactionFlux95 = PARAM(118) * (SPVAR(4) / PARAM(123) - SPVAR(15) / PARAM(124));

    realtype ReactionFlux96 = PARAM(119) * (SPVAR(4) / PARAM(123) - SPVAR(34) / PARAM(125));

    realtype ReactionFlux97 = PARAM(120) * (SPVAR(4) / PARAM(123) - SPVAR(60) / PARAM(126));

    realtype ReactionFlux98 = PARAM(121) * AUX_VAR_V_T * SPVAR(34) / PARAM(125);

    realtype ReactionFlux99 = PARAM(121) * AUX_VAR_V_T * SPVAR(60) / PARAM(126);

    realtype ReactionFlux100 = PARAM(122) * SPVAR(4);

    realtype ReactionFlux101 = PARAM(127) * (SPVAR(5) / PARAM(132) - SPVAR(16) / PARAM(133));

    realtype ReactionFlux102 = PARAM(128) * (SPVAR(5) / PARAM(132) - SPVAR(35) / PARAM(134));

    realtype ReactionFlux103 = PARAM(129) * (SPVAR(5) / PARAM(132) - SPVAR(61) / PARAM(135));

    realtype ReactionFlux104 = PARAM(130) * AUX_VAR_V_T * SPVAR(35) / PARAM(134);

    realtype ReactionFlux105 = PARAM(130) * AUX_VAR_V_T * SPVAR(61) / PARAM(135);

    realtype ReactionFlux106 = PARAM(131) * SPVAR(5);

    realtype ReactionFlux107 = PARAM(136) * SPVAR(5) / (SPVAR(5) + PARAM(137));

    realtype ReactionFlux108 = PARAM(138) * (SPVAR(6) / PARAM(143) - SPVAR(17) / PARAM(144));

    realtype ReactionFlux109 = PARAM(139) * (SPVAR(6) / PARAM(143) - SPVAR(36) / PARAM(145));

    realtype ReactionFlux110 = PARAM(140) * (SPVAR(6) / PARAM(143) - SPVAR(62) / PARAM(146));

    realtype ReactionFlux111 = PARAM(141) * AUX_VAR_V_T * SPVAR(36) / PARAM(145);

    realtype ReactionFlux112 = PARAM(141) * AUX_VAR_V_T * SPVAR(62) / PARAM(146);

    realtype ReactionFlux113 = PARAM(142) * SPVAR(6);

    realtype ReactionFlux114 = PARAM(148) * SPVAR(29) / (SPVAR(29) + PARAM(201)) * (1.0 - AUX_VAR_syn_T_C1_PDL1_total / (PARAM(183) * PARAM(150) / PARAM(95)));

    realtype ReactionFlux115 = PARAM(148) * PARAM(184) * SPVAR(29) / (SPVAR(29) + PARAM(201)) * (1.0 - AUX_VAR_syn_T_C1_PDL2_total / (PARAM(183) * PARAM(150) / PARAM(95) * PARAM(184)));

    realtype ReactionFlux116 = PARAM(149) * (PARAM(183) / PARAM(95) - AUX_VAR_syn_T_C1_PDL1_total) * PARAM(6);

    realtype ReactionFlux117 = PARAM(149) * (PARAM(183) / PARAM(95) * PARAM(184) - AUX_VAR_syn_T_C1_PDL2_total) * PARAM(6);

    realtype ReactionFlux118 = (PARAM(147) * SPVAR(78) * SPVAR(79) - PARAM(161) * SPVAR(76)) * PARAM(6);

    realtype ReactionFlux119 = (PARAM(151) * SPVAR(78) * SPVAR(80) - PARAM(162) * SPVAR(77)) * PARAM(6);

    realtype ReactionFlux120 = (2.0 * PARAM(152) * (SPVAR(78) * SPVAR(34) / PARAM(125)) - PARAM(163) * SPVAR(81)) * PARAM(6);

    realtype ReactionFlux121 = (PARAM(172) * PARAM(152) * SPVAR(78) * SPVAR(81) - 2.0 * PARAM(163) * SPVAR(82)) * PARAM(6);

    realtype ReactionFlux122 = (2.0 * PARAM(153) * (SPVAR(79) * SPVAR(35) / PARAM(134)) - PARAM(164) * SPVAR(83)) * PARAM(6);

    realtype ReactionFlux123 = (PARAM(173) * PARAM(153) * SPVAR(79) * SPVAR(83) - 2.0 * PARAM(164) * SPVAR(84)) * PARAM(6);

    realtype ReactionFlux124 = (2.0 * PARAM(154) * SPVAR(100) * SPVAR(102) - PARAM(165) * SPVAR(88)) * PARAM(6);

    realtype ReactionFlux125 = (PARAM(154) * SPVAR(100) * SPVAR(88) - 2.0 * PARAM(165) * SPVAR(89)) * PARAM(6);

    realtype ReactionFlux126 = (PARAM(155) * SPVAR(100) * SPVAR(104) - PARAM(166) * SPVAR(90)) * PARAM(6);

    realtype ReactionFlux127 = (4.0 * PARAM(156) * SPVAR(101) * SPVAR(102) - PARAM(167) * SPVAR(91)) * PARAM(6);

    realtype ReactionFlux128 = (2.0 * PARAM(156) * SPVAR(101) * SPVAR(91) - 2.0 * PARAM(167) * SPVAR(93)) * PARAM(6);

    realtype ReactionFlux129 = (4.0 * PARAM(156) * SPVAR(102) * SPVAR(93) - PARAM(167) * SPVAR(94)) * PARAM(6);

    realtype ReactionFlux130 = (2.0 * PARAM(156) * SPVAR(91) * SPVAR(102) - 2.0 * PARAM(167) * SPVAR(92)) * PARAM(6);

    realtype ReactionFlux131 = (4.0 * PARAM(156) * SPVAR(101) * SPVAR(92) - PARAM(167) * SPVAR(94)) * PARAM(6);

    realtype ReactionFlux132 = (2.0 * PARAM(157) * SPVAR(101) * SPVAR(104) - PARAM(168) * SPVAR(95)) * PARAM(6);

    realtype ReactionFlux133 = (PARAM(157) * SPVAR(95) * SPVAR(104) - 2.0 * PARAM(168) * SPVAR(96)) * PARAM(6);

    realtype ReactionFlux134 = (4.0 * PARAM(159) * (SPVAR(101) * SPVAR(36) / PARAM(145)) - PARAM(170) * SPVAR(105)) * PARAM(6);

    realtype ReactionFlux135 = (2.0 * PARAM(174) * PARAM(159) * SPVAR(101) * SPVAR(105) - 2.0 * PARAM(170) * SPVAR(106)) * PARAM(6);

    realtype ReactionFlux136 = (PARAM(160) * SPVAR(103) * SPVAR(103) - PARAM(171) * SPVAR(102)) * PARAM(6);

    realtype ReactionFlux137 = (PARAM(158) * SPVAR(103) * SPVAR(79) - PARAM(169) * SPVAR(97)) * PARAM(6);

    realtype ReactionFlux138 = (PARAM(154) * SPVAR(97) * SPVAR(100) - PARAM(165) * SPVAR(98)) * PARAM(6);

    realtype ReactionFlux139 = (2.0 * PARAM(156) * SPVAR(97) * SPVAR(101) - PARAM(167) * SPVAR(99)) * PARAM(6);

    realtype ReactionFlux140 = (2.0 * PARAM(153) * (SPVAR(85) * SPVAR(35) / PARAM(134)) - PARAM(164) * SPVAR(86)) * PARAM(6);

    realtype ReactionFlux141 = (PARAM(173) * PARAM(153) * SPVAR(85) * SPVAR(86) - 2.0 * PARAM(164) * SPVAR(87)) * PARAM(6);

    realtype ReactionFlux142 = PARAM(148) * SPVAR(29) / (SPVAR(29) + PARAM(201)) * (1.0 - AUX_VAR_syn_T_APC_PDL1_total / (PARAM(187) * PARAM(150) / PARAM(95)));

    realtype ReactionFlux143 = PARAM(148) * PARAM(188) * SPVAR(29) / (SPVAR(29) + PARAM(201)) * (1.0 - AUX_VAR_syn_T_APC_PDL2_total / (PARAM(187) * PARAM(150) / PARAM(95) * PARAM(188)));

    realtype ReactionFlux144 = PARAM(149) * (PARAM(187) / PARAM(95) - AUX_VAR_syn_T_APC_PDL1_total) * PARAM(7);

    realtype ReactionFlux145 = PARAM(149) * (PARAM(187) / PARAM(95) * PARAM(188) - AUX_VAR_syn_T_APC_PDL2_total) * PARAM(7);

    realtype ReactionFlux146 = (PARAM(147) * SPVAR(109) * SPVAR(110) - PARAM(161) * SPVAR(107)) * PARAM(7);

    realtype ReactionFlux147 = (PARAM(151) * SPVAR(109) * SPVAR(111) - PARAM(162) * SPVAR(108)) * PARAM(7);

    realtype ReactionFlux148 = (2.0 * PARAM(152) * (SPVAR(109) * SPVAR(60) / PARAM(126)) - PARAM(163) * SPVAR(112)) * PARAM(7);

    realtype ReactionFlux149 = (PARAM(172) * PARAM(152) * SPVAR(109) * SPVAR(112) - 2.0 * PARAM(163) * SPVAR(113)) * PARAM(7);

    realtype ReactionFlux150 = (2.0 * PARAM(153) * (SPVAR(110) * SPVAR(61) / PARAM(135)) - PARAM(164) * SPVAR(114)) * PARAM(7);

    realtype ReactionFlux151 = (PARAM(173) * PARAM(153) * SPVAR(110) * SPVAR(114) - 2.0 * PARAM(164) * SPVAR(115)) * PARAM(7);

    realtype ReactionFlux152 = (2.0 * PARAM(154) * SPVAR(131) * SPVAR(133) - PARAM(165) * SPVAR(119)) * PARAM(7);

    realtype ReactionFlux153 = (PARAM(154) * SPVAR(131) * SPVAR(119) - 2.0 * PARAM(165) * SPVAR(120)) * PARAM(7);

    realtype ReactionFlux154 = (PARAM(155) * SPVAR(131) * SPVAR(135) - PARAM(166) * SPVAR(121)) * PARAM(7);

    realtype ReactionFlux155 = (4.0 * PARAM(156) * SPVAR(132) * SPVAR(133) - PARAM(167) * SPVAR(122)) * PARAM(7);

    realtype ReactionFlux156 = (2.0 * PARAM(156) * SPVAR(132) * SPVAR(122) - 2.0 * PARAM(167) * SPVAR(124)) * PARAM(7);

    realtype ReactionFlux157 = (4.0 * PARAM(156) * SPVAR(133) * SPVAR(124) - PARAM(167) * SPVAR(125)) * PARAM(7);

    realtype ReactionFlux158 = (2.0 * PARAM(156) * SPVAR(122) * SPVAR(133) - 2.0 * PARAM(167) * SPVAR(123)) * PARAM(7);

    realtype ReactionFlux159 = (4.0 * PARAM(156) * SPVAR(132) * SPVAR(123) - PARAM(167) * SPVAR(125)) * PARAM(7);

    realtype ReactionFlux160 = (2.0 * PARAM(157) * SPVAR(132) * SPVAR(135) - PARAM(168) * SPVAR(126)) * PARAM(7);

    realtype ReactionFlux161 = (PARAM(157) * SPVAR(126) * SPVAR(135) - 2.0 * PARAM(168) * SPVAR(127)) * PARAM(7);

    realtype ReactionFlux162 = (4.0 * PARAM(159) * (SPVAR(132) * SPVAR(62) / PARAM(146)) - PARAM(170) * SPVAR(136)) * PARAM(7);

    realtype ReactionFlux163 = (2.0 * PARAM(174) * PARAM(159) * SPVAR(132) * SPVAR(136) - 2.0 * PARAM(170) * SPVAR(137)) * PARAM(7);

    realtype ReactionFlux164 = (PARAM(160) * SPVAR(134) * SPVAR(134) - PARAM(171) * SPVAR(133)) * PARAM(7);

    realtype ReactionFlux165 = (PARAM(158) * SPVAR(134) * SPVAR(110) - PARAM(169) * SPVAR(128)) * PARAM(7);

    realtype ReactionFlux166 = (PARAM(154) * SPVAR(128) * SPVAR(131) - PARAM(165) * SPVAR(129)) * PARAM(7);

    realtype ReactionFlux167 = (2.0 * PARAM(156) * SPVAR(128) * SPVAR(132) - PARAM(167) * SPVAR(130)) * PARAM(7);

    realtype ReactionFlux168 = (2.0 * PARAM(153) * (SPVAR(116) * SPVAR(61) / PARAM(135)) - PARAM(164) * SPVAR(117)) * PARAM(7);

    realtype ReactionFlux169 = (PARAM(173) * PARAM(153) * SPVAR(116) * SPVAR(117) - 2.0 * PARAM(164) * SPVAR(118)) * PARAM(7);

    realtype ReactionFlux170 = PARAM(191) * AUX_VAR_H_APCh * AUX_VAR_H_P1 * SPVAR(51);

    realtype ReactionFlux171 = PARAM(191) * AUX_VAR_H_APCh * AUX_VAR_H_P1 * SPVAR(51) * PARAM(53);

    realtype ReactionFlux172 = PARAM(31) / AUX_VAR_N_aTh * SPVAR(63);

    realtype ReactionFlux173 = PARAM(31) / AUX_VAR_N_aTh * std::pow(2.0, AUX_VAR_N_aTh) * SPVAR(63);

    realtype ReactionFlux174 = PARAM(192) * SPVAR(37) * AUX_VAR_H_TGFb * AUX_VAR_H_ArgI_Treg;

    realtype ReactionFlux175 = PARAM(32) * SPVAR(7);

    realtype ReactionFlux176 = PARAM(32) * SPVAR(18);

    realtype ReactionFlux177 = PARAM(32) * SPVAR(64);

    realtype ReactionFlux178 = PARAM(32) * SPVAR(37);

    realtype ReactionFlux179 = PARAM(9) * SPVAR(37) * (PARAM(197) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux180 = PARAM(33) * SPVAR(7);

    realtype ReactionFlux181 = PARAM(34) * SPVAR(18);

    realtype ReactionFlux182 = PARAM(35) * AUX_VAR_V_T * SPVAR(7) * (std::pow(AUX_VAR_C_total, 2.0) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux183 = PARAM(28) * SPVAR(64);

    realtype ReactionFlux184 = PARAM(44) * SPVAR(63);

    realtype ReactionFlux185 = PARAM(194) * (PARAM(198) - SPVAR(38)) * AUX_VAR_V_T;

    realtype ReactionFlux186 = PARAM(193) * SPVAR(27);

    realtype ReactionFlux187 = PARAM(199) * SPVAR(37);

    realtype ReactionFlux188 = PARAM(200) * SPVAR(29) * AUX_VAR_V_T;

    realtype ReactionFlux189 = PARAM(202) * AUX_VAR_C_total;

    realtype ReactionFlux190 = PARAM(203) * SPVAR(42) * AUX_VAR_V_T;

    realtype ReactionFlux191 = PARAM(205) * AUX_VAR_V_T * (SPVAR(42) / (SPVAR(42) + PARAM(204)));

    realtype ReactionFlux192 = PARAM(206) * SPVAR(39);

    realtype ReactionFlux193 = PARAM(9) * SPVAR(39) * (PARAM(197) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux194 = PARAM(207) * SPVAR(40) * AUX_VAR_V_T;

    realtype ReactionFlux195 = PARAM(208) * SPVAR(41) * AUX_VAR_V_T;

    realtype ReactionFlux196 = PARAM(209) * SPVAR(39);

    realtype ReactionFlux197 = PARAM(210) * SPVAR(39);

    realtype ReactionFlux198 = PARAM(220) * PARAM(215) * SPVAR(10) * PARAM(0);

    realtype ReactionFlux199 = PARAM(220) * PARAM(214) * SPVAR(9) * PARAM(0);

    realtype ReactionFlux200 = PARAM(221) * (SPVAR(8) / PARAM(225) - SPVAR(19) / PARAM(226)) * PARAM(0);

    realtype ReactionFlux201 = PARAM(222) * (SPVAR(8) / PARAM(225) - SPVAR(43) / PARAM(227)) * PARAM(0);

    realtype ReactionFlux202 = PARAM(223) * (SPVAR(8) / PARAM(225) - SPVAR(65) / PARAM(228)) * PARAM(0);

    realtype ReactionFlux203 = PARAM(224) * SPVAR(43) / PARAM(227) * AUX_VAR_V_T;

    realtype ReactionFlux204 = PARAM(224) * SPVAR(65) / PARAM(228) * PARAM(2);

    realtype ReactionFlux205 = PARAM(216) * SPVAR(8) / (SPVAR(8) + PARAM(217)) * PARAM(0);

    realtype ReactionFlux206 = PARAM(233) * SPVAR(23) * AUX_VAR_R_cabo;

    realtype ReactionFlux207 = PARAM(233) * SPVAR(23);

    realtype ReactionFlux208 = PARAM(238) * SPVAR(45);

    realtype ReactionFlux209 = PARAM(235) * AUX_VAR_V_T * (SPVAR(42) / (SPVAR(42) + PARAM(204)));

    realtype ReactionFlux210 = PARAM(9) * SPVAR(44) * (PARAM(197) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux211 = PARAM(9) * SPVAR(45) * (PARAM(197) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux212 = PARAM(236) * SPVAR(44);

    realtype ReactionFlux213 = PARAM(236) * SPVAR(45);

    realtype ReactionFlux214 = PARAM(239) * SPVAR(31);

    realtype ReactionFlux215 = PARAM(240) * SPVAR(44);

    realtype ReactionFlux216 = PARAM(241) * SPVAR(46) * AUX_VAR_V_T;

    realtype ReactionFlux217 = PARAM(237) * SPVAR(45);

    realtype ReactionFlux218 = PARAM(242) * SPVAR(45);

    realtype ReactionFlux219 = PARAM(243) * SPVAR(47) * AUX_VAR_V_T;

    realtype ReactionFlux220 = PARAM(244) * SPVAR(44) * (SPVAR(38) / (SPVAR(38) + PARAM(195)) + SPVAR(47) / (SPVAR(47) + PARAM(246)));

    realtype ReactionFlux221 = PARAM(245) * SPVAR(45) * (SPVAR(46) / (SPVAR(46) + PARAM(247)) + SPVAR(29) / (SPVAR(29) + PARAM(248)));

    realtype ReactionFlux222 = PARAM(148) * SPVAR(29) / (SPVAR(29) + PARAM(201)) * (1.0 - AUX_VAR_PDL1_total / (PARAM(183) * PARAM(150) / PARAM(95)));

    realtype ReactionFlux223 = PARAM(148) * PARAM(184) * SPVAR(29) / (SPVAR(29) + PARAM(201)) * (1.0 - AUX_VAR_PDL2_total / (PARAM(183) * PARAM(150) / PARAM(95) * PARAM(184)));

    realtype ReactionFlux224 = PARAM(149) * (PARAM(183) / PARAM(95) - AUX_VAR_PDL1_total) * PARAM(8);

    realtype ReactionFlux225 = PARAM(149) * (PARAM(183) / PARAM(95) * PARAM(184) - AUX_VAR_PDL2_total) * PARAM(8);

    realtype ReactionFlux226 = (PARAM(251) * SPVAR(138) * SPVAR(139) - PARAM(252) * SPVAR(140)) * PARAM(8);

    realtype ReactionFlux227 = (PARAM(147) * SPVAR(143) * SPVAR(144) - PARAM(161) * SPVAR(141)) * PARAM(8);

    realtype ReactionFlux228 = (PARAM(151) * SPVAR(143) * SPVAR(145) - PARAM(162) * SPVAR(142)) * PARAM(8);

    realtype ReactionFlux229 = (2.0 * PARAM(152) * (SPVAR(143) * SPVAR(34) / PARAM(125)) - PARAM(163) * SPVAR(146)) * PARAM(8);

    realtype ReactionFlux230 = (PARAM(172) * PARAM(152) * SPVAR(143) * SPVAR(146) - 2.0 * PARAM(163) * SPVAR(147)) * PARAM(8);

    realtype ReactionFlux231 = (2.0 * PARAM(153) * (SPVAR(144) * SPVAR(35) / PARAM(134)) - PARAM(164) * SPVAR(148)) * PARAM(8);

    realtype ReactionFlux232 = (PARAM(173) * PARAM(153) * SPVAR(144) * SPVAR(148) - 2.0 * PARAM(164) * SPVAR(149)) * PARAM(8);

    realtype ReactionFlux233 = (PARAM(160) * SPVAR(152) * SPVAR(152) - PARAM(171) * SPVAR(151)) * PARAM(8);

    realtype ReactionFlux234 = (PARAM(158) * SPVAR(152) * SPVAR(144) - PARAM(169) * SPVAR(150)) * PARAM(8);

    realtype ReactionFlux235 = PARAM(249) * SPVAR(23) * SPVAR(44) / (SPVAR(44) + PARAM(260) * AUX_VAR_C_total + PARAM(10)) * (1.0 - AUX_VAR_H_Mac_C) * (1.0 - AUX_VAR_H_IL10_phago);

    realtype ReactionFlux236 = PARAM(249) * SPVAR(26) * SPVAR(44) / (SPVAR(44) + PARAM(260) * AUX_VAR_C_total + PARAM(10)) * (1.0 - AUX_VAR_H_Mac_C) * (1.0 - AUX_VAR_H_IL10_phago);

    realtype ReactionFlux237 = PARAM(238) * SPVAR(45);

    realtype ReactionFlux238 = PARAM(261) * AUX_VAR_V_T * ((PARAM(262) - SPVAR(48) / AUX_VAR_V_T) / PARAM(262));

    realtype ReactionFlux239 = PARAM(263) * SPVAR(48) * (SPVAR(38) / (SPVAR(38) + PARAM(195)));

    realtype ReactionFlux240 = PARAM(9) * SPVAR(48) * (PARAM(197) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux241 = PARAM(9) * SPVAR(49) * (PARAM(197) / (std::pow(AUX_VAR_C_total, 2.0) + PARAM(197)));

    realtype ReactionFlux242 = PARAM(271) * SPVAR(48);

    realtype ReactionFlux243 = PARAM(272) * SPVAR(49);

    realtype ReactionFlux244 = PARAM(264) * ((PARAM(268) - AUX_VAR_ECM_level) / PARAM(268)) * SPVAR(48) * (SPVAR(38) / (SPVAR(38) + PARAM(195)));

    realtype ReactionFlux245 = PARAM(265) * ((PARAM(268) - AUX_VAR_ECM_level) / PARAM(268)) * SPVAR(49) * (SPVAR(38) / (SPVAR(38) + PARAM(195)));

    realtype ReactionFlux246 = PARAM(266) * (PARAM(267) - AUX_VAR_ECM_level) * AUX_VAR_V_T;

    //dydt:

    //d(V_C.nT0)/dt
    NV_DATA_S(ydot)[0] = ReactionFlux12 - ReactionFlux16 - ReactionFlux18 + ReactionFlux19 - ReactionFlux20 + ReactionFlux21;

    //d(V_C.T0)/dt
    NV_DATA_S(ydot)[1] =  - ReactionFlux26 - ReactionFlux31 + ReactionFlux32 - ReactionFlux33 + ReactionFlux34;

    //d(V_C.nT1)/dt
    NV_DATA_S(ydot)[2] = ReactionFlux38 - ReactionFlux42 - ReactionFlux44 + ReactionFlux45 - ReactionFlux46 + ReactionFlux47;

    //d(V_C.T1)/dt
    NV_DATA_S(ydot)[3] =  - ReactionFlux52 - ReactionFlux60 + ReactionFlux61 - ReactionFlux62 + ReactionFlux63;

    //d(V_C.aPD1)/dt
    NV_DATA_S(ydot)[4] = 1/PARAM(0)*( - ReactionFlux95 - ReactionFlux96 - ReactionFlux97 + ReactionFlux99 - ReactionFlux100);

    //d(V_C.aPDL1)/dt
    NV_DATA_S(ydot)[5] = 1/PARAM(0)*( - ReactionFlux101 - ReactionFlux102 - ReactionFlux103 + ReactionFlux105 - ReactionFlux106 - ReactionFlux107);

    //d(V_C.aCTLA4)/dt
    NV_DATA_S(ydot)[6] = 1/PARAM(0)*( - ReactionFlux108 - ReactionFlux109 - ReactionFlux110 + ReactionFlux112 - ReactionFlux113);

    //d(V_C.Th)/dt
    NV_DATA_S(ydot)[7] =  - ReactionFlux175 - ReactionFlux180 + ReactionFlux181 - ReactionFlux182 + ReactionFlux183;

    //d(V_C.cabozantinib)/dt
    NV_DATA_S(ydot)[8] = 1/PARAM(0)*(ReactionFlux198 + ReactionFlux199 - ReactionFlux200 - ReactionFlux201 - ReactionFlux202 + ReactionFlux204 - ReactionFlux205);

    //d(V_C.A_site1)/dt
    NV_DATA_S(ydot)[9] = 1/PARAM(0)*( - ReactionFlux199);

    //d(V_C.A_site2)/dt
    NV_DATA_S(ydot)[10] = 1/PARAM(0)*( - ReactionFlux198);

    //d(V_P.nT0)/dt
    NV_DATA_S(ydot)[11] = ReactionFlux13 - ReactionFlux15 + ReactionFlux18 - ReactionFlux19;

    //d(V_P.T0)/dt
    NV_DATA_S(ydot)[12] =  - ReactionFlux27 + ReactionFlux31 - ReactionFlux32;

    //d(V_P.nT1)/dt
    NV_DATA_S(ydot)[13] = ReactionFlux39 - ReactionFlux41 + ReactionFlux44 - ReactionFlux45;

    //d(V_P.T1)/dt
    NV_DATA_S(ydot)[14] =  - ReactionFlux53 + ReactionFlux60 - ReactionFlux61;

    //d(V_P.aPD1)/dt
    NV_DATA_S(ydot)[15] = 1/PARAM(1)*(ReactionFlux95);

    //d(V_P.aPDL1)/dt
    NV_DATA_S(ydot)[16] = 1/PARAM(1)*(ReactionFlux101);

    //d(V_P.aCTLA4)/dt
    NV_DATA_S(ydot)[17] = 1/PARAM(1)*(ReactionFlux108);

    //d(V_P.Th)/dt
    NV_DATA_S(ydot)[18] =  - ReactionFlux176 + ReactionFlux180 - ReactionFlux181;

    //d(V_P.cabozantinib)/dt
    NV_DATA_S(ydot)[19] = 1/PARAM(1)*(ReactionFlux200);

    //d(V_T.C_x)/dt
    NV_DATA_S(ydot)[20] =  - ReactionFlux1 + ReactionFlux5 + ReactionFlux11 + ReactionFlux66 + ReactionFlux67 + ReactionFlux235 + ReactionFlux236;

    //d(V_T.T1_exh)/dt
    NV_DATA_S(ydot)[21] =  - ReactionFlux2 + ReactionFlux55 + ReactionFlux56 + ReactionFlux57 + ReactionFlux58 + ReactionFlux59;

    //d(V_T.Th_exh)/dt
    NV_DATA_S(ydot)[22] =  - ReactionFlux3 + ReactionFlux178 + ReactionFlux179;

    //d(V_T.C1)/dt
    NV_DATA_S(ydot)[23] = ReactionFlux4 - ReactionFlux5 - ReactionFlux66 - ReactionFlux206 - ReactionFlux207 - ReactionFlux235;

    //d(V_T.K)/dt
    NV_DATA_S(ydot)[24] = ReactionFlux8 - ReactionFlux9;

    //d(V_T.c_vas)/dt
    NV_DATA_S(ydot)[25] = 1/AUX_VAR_V_T*(ReactionFlux6 - ReactionFlux7 + ReactionFlux208 + ReactionFlux237);

    //d(V_T.C2)/dt
    NV_DATA_S(ydot)[26] = ReactionFlux10 - ReactionFlux11 - ReactionFlux67 + ReactionFlux206 + ReactionFlux207 - ReactionFlux236;

    //d(V_T.T0)/dt
    NV_DATA_S(ydot)[27] =  - ReactionFlux29 - ReactionFlux30 + ReactionFlux33 + ReactionFlux174;

    //d(V_T.T1)/dt
    NV_DATA_S(ydot)[28] =  - ReactionFlux55 - ReactionFlux56 - ReactionFlux57 - ReactionFlux58 - ReactionFlux59 + ReactionFlux62;

    //d(V_T.IFNg)/dt
    NV_DATA_S(ydot)[29] = 1/AUX_VAR_V_T*(ReactionFlux65 + ReactionFlux187 - ReactionFlux188);

    //d(V_T.APC)/dt
    NV_DATA_S(ydot)[30] = ReactionFlux68 - ReactionFlux70;

    //d(V_T.mAPC)/dt
    NV_DATA_S(ydot)[31] = ReactionFlux70 - ReactionFlux71 - ReactionFlux72;

    //d(V_T.P0)/dt
    NV_DATA_S(ydot)[32] = 1/AUX_VAR_V_T*(ReactionFlux75 - ReactionFlux76 - ReactionFlux77);

    //d(V_T.P1)/dt
    NV_DATA_S(ydot)[33] = 1/AUX_VAR_V_T*(ReactionFlux85 - ReactionFlux86 - ReactionFlux87);

    //d(V_T.aPD1)/dt
    NV_DATA_S(ydot)[34] = 1/AUX_VAR_V_T*(ReactionFlux96 - ReactionFlux98);

    //d(V_T.aPDL1)/dt
    NV_DATA_S(ydot)[35] = 1/AUX_VAR_V_T*(ReactionFlux102 - ReactionFlux104);

    //d(V_T.aCTLA4)/dt
    NV_DATA_S(ydot)[36] = 1/AUX_VAR_V_T*(ReactionFlux109 - ReactionFlux111);

    //d(V_T.Th)/dt
    NV_DATA_S(ydot)[37] =  - ReactionFlux174 - ReactionFlux178 - ReactionFlux179 + ReactionFlux182;

    //d(V_T.TGFb)/dt
    NV_DATA_S(ydot)[38] = 1/AUX_VAR_V_T*(ReactionFlux185 + ReactionFlux186 + ReactionFlux217);

    //d(V_T.MDSC)/dt
    NV_DATA_S(ydot)[39] = ReactionFlux191 - ReactionFlux192 - ReactionFlux193;

    //d(V_T.NO)/dt
    NV_DATA_S(ydot)[40] = 1/AUX_VAR_V_T*( - ReactionFlux194 + ReactionFlux196);

    //d(V_T.ArgI)/dt
    NV_DATA_S(ydot)[41] = 1/AUX_VAR_V_T*( - ReactionFlux195 + ReactionFlux197);

    //d(V_T.CCL2)/dt
    NV_DATA_S(ydot)[42] = 1/AUX_VAR_V_T*(ReactionFlux189 - ReactionFlux190);

    //d(V_T.cabozantinib)/dt
    NV_DATA_S(ydot)[43] = 1/AUX_VAR_V_T*(ReactionFlux201 - ReactionFlux203);

    //d(V_T.Mac_M1)/dt
    NV_DATA_S(ydot)[44] = ReactionFlux209 - ReactionFlux210 - ReactionFlux212 - ReactionFlux220 + ReactionFlux221;

    //d(V_T.Mac_M2)/dt
    NV_DATA_S(ydot)[45] =  - ReactionFlux211 - ReactionFlux213 + ReactionFlux220 - ReactionFlux221;

    //d(V_T.IL12)/dt
    NV_DATA_S(ydot)[46] = 1/AUX_VAR_V_T*(ReactionFlux214 + ReactionFlux215 - ReactionFlux216);

    //d(V_T.IL10)/dt
    NV_DATA_S(ydot)[47] = 1/AUX_VAR_V_T*(ReactionFlux218 - ReactionFlux219);

    //d(V_T.Fib)/dt
    NV_DATA_S(ydot)[48] = ReactionFlux238 - ReactionFlux239 - ReactionFlux240 - ReactionFlux242;

    //d(V_T.CAF)/dt
    NV_DATA_S(ydot)[49] = ReactionFlux239 - ReactionFlux241 - ReactionFlux243;

    //d(V_T.ECM)/dt
    NV_DATA_S(ydot)[50] = ReactionFlux244 + ReactionFlux245 + ReactionFlux246;

    //d(V_LN.nT0)/dt
    NV_DATA_S(ydot)[51] = ReactionFlux14 - ReactionFlux17 + ReactionFlux20 - ReactionFlux21 - ReactionFlux22 - ReactionFlux170;

    //d(V_LN.aT0)/dt
    NV_DATA_S(ydot)[52] = ReactionFlux23 - ReactionFlux24;

    //d(V_LN.T0)/dt
    NV_DATA_S(ydot)[53] = ReactionFlux25 - ReactionFlux28 - ReactionFlux34;

    //d(V_LN.IL2)/dt
    NV_DATA_S(ydot)[54] = 1/PARAM(2)*( - ReactionFlux35 - ReactionFlux36 - ReactionFlux37 + ReactionFlux64 + ReactionFlux184);

    //d(V_LN.nT1)/dt
    NV_DATA_S(ydot)[55] = ReactionFlux40 - ReactionFlux43 + ReactionFlux46 - ReactionFlux47 - ReactionFlux48;

    //d(V_LN.aT1)/dt
    NV_DATA_S(ydot)[56] = ReactionFlux49 - ReactionFlux50;

    //d(V_LN.T1)/dt
    NV_DATA_S(ydot)[57] = ReactionFlux51 - ReactionFlux54 - ReactionFlux63;

    //d(V_LN.APC)/dt
    NV_DATA_S(ydot)[58] = ReactionFlux69;

    //d(V_LN.mAPC)/dt
    NV_DATA_S(ydot)[59] = ReactionFlux71 - ReactionFlux73;

    //d(V_LN.aPD1)/dt
    NV_DATA_S(ydot)[60] = 1/PARAM(2)*(ReactionFlux97 + ReactionFlux98 - ReactionFlux99);

    //d(V_LN.aPDL1)/dt
    NV_DATA_S(ydot)[61] = 1/PARAM(2)*(ReactionFlux103 + ReactionFlux104 - ReactionFlux105);

    //d(V_LN.aCTLA4)/dt
    NV_DATA_S(ydot)[62] = 1/PARAM(2)*(ReactionFlux110 + ReactionFlux111 - ReactionFlux112);

    //d(V_LN.aTh)/dt
    NV_DATA_S(ydot)[63] = ReactionFlux171 - ReactionFlux172;

    //d(V_LN.Th)/dt
    NV_DATA_S(ydot)[64] = ReactionFlux173 - ReactionFlux177 - ReactionFlux183;

    //d(V_LN.cabozantinib)/dt
    NV_DATA_S(ydot)[65] = 1/PARAM(2)*(ReactionFlux202 + ReactionFlux203 - ReactionFlux204);

    //d(V_e.P0)/dt
    NV_DATA_S(ydot)[66] = 1/PARAM(3)*(ReactionFlux78 - ReactionFlux79);

    //d(V_e.p0)/dt
    NV_DATA_S(ydot)[67] = 1/PARAM(3)*(ReactionFlux79 - ReactionFlux80 - ReactionFlux81 + ReactionFlux82);

    //d(V_e.P1)/dt
    NV_DATA_S(ydot)[68] = 1/PARAM(3)*(ReactionFlux88 - ReactionFlux89);

    //d(V_e.p1)/dt
    NV_DATA_S(ydot)[69] = 1/PARAM(3)*(ReactionFlux89 - ReactionFlux90 - ReactionFlux91 + ReactionFlux92);

    //d(A_e.M1)/dt
    NV_DATA_S(ydot)[70] = 1/PARAM(4)*( - ReactionFlux74 - ReactionFlux81 + ReactionFlux82 - ReactionFlux91 + ReactionFlux92);

    //d(A_e.M1p0)/dt
    NV_DATA_S(ydot)[71] = 1/PARAM(4)*(ReactionFlux81 - ReactionFlux82 - ReactionFlux84);

    //d(A_e.M1p1)/dt
    NV_DATA_S(ydot)[72] = 1/PARAM(4)*(ReactionFlux91 - ReactionFlux92 - ReactionFlux94);

    //d(A_s.M1)/dt
    NV_DATA_S(ydot)[73] = 1/PARAM(5)*(ReactionFlux74 + ReactionFlux83 + ReactionFlux93);

    //d(A_s.M1p0)/dt
    NV_DATA_S(ydot)[74] = 1/PARAM(5)*( - ReactionFlux83 + ReactionFlux84);

    //d(A_s.M1p1)/dt
    NV_DATA_S(ydot)[75] = 1/PARAM(5)*( - ReactionFlux93 + ReactionFlux94);

    //d(syn_T_C1.PD1_PDL1)/dt
    NV_DATA_S(ydot)[76] = 1/PARAM(6)*(ReactionFlux118);

    //d(syn_T_C1.PD1_PDL2)/dt
    NV_DATA_S(ydot)[77] = 1/PARAM(6)*(ReactionFlux119);

    //d(syn_T_C1.PD1)/dt
    NV_DATA_S(ydot)[78] = 1/PARAM(6)*( - ReactionFlux118 - ReactionFlux119 - ReactionFlux120 - ReactionFlux121);

    //d(syn_T_C1.PDL1)/dt
    NV_DATA_S(ydot)[79] = 1/PARAM(6)*(ReactionFlux114 + ReactionFlux116 - ReactionFlux118 - ReactionFlux122 - ReactionFlux123 - ReactionFlux137);

    //d(syn_T_C1.PDL2)/dt
    NV_DATA_S(ydot)[80] = 1/PARAM(6)*(ReactionFlux115 + ReactionFlux117 - ReactionFlux119);

    //d(syn_T_C1.PD1_aPD1)/dt
    NV_DATA_S(ydot)[81] = 1/PARAM(6)*(ReactionFlux120 - ReactionFlux121);

    //d(syn_T_C1.PD1_aPD1_PD1)/dt
    NV_DATA_S(ydot)[82] = 1/PARAM(6)*(ReactionFlux121);

    //d(syn_T_C1.PDL1_aPDL1)/dt
    NV_DATA_S(ydot)[83] = 1/PARAM(6)*(ReactionFlux122 - ReactionFlux123);

    //d(syn_T_C1.PDL1_aPDL1_PDL1)/dt
    NV_DATA_S(ydot)[84] = 1/PARAM(6)*(ReactionFlux123);

    //d(syn_T_C1.TPDL1)/dt
    NV_DATA_S(ydot)[85] = 1/PARAM(6)*( - ReactionFlux140 - ReactionFlux141);

    //d(syn_T_C1.TPDL1_aPDL1)/dt
    NV_DATA_S(ydot)[86] = 1/PARAM(6)*(ReactionFlux140 - ReactionFlux141);

    //d(syn_T_C1.TPDL1_aPDL1_TPDL1)/dt
    NV_DATA_S(ydot)[87] = 1/PARAM(6)*(ReactionFlux141);

    //d(syn_T_C1.CD28_CD80)/dt
    NV_DATA_S(ydot)[88] = 1/PARAM(6)*(ReactionFlux124 - ReactionFlux125);

    //d(syn_T_C1.CD28_CD80_CD28)/dt
    NV_DATA_S(ydot)[89] = 1/PARAM(6)*(ReactionFlux125);

    //d(syn_T_C1.CD28_CD86)/dt
    NV_DATA_S(ydot)[90] = 1/PARAM(6)*(ReactionFlux126);

    //d(syn_T_C1.CD80_CTLA4)/dt
    NV_DATA_S(ydot)[91] = 1/PARAM(6)*(ReactionFlux127 - ReactionFlux128 - ReactionFlux130);

    //d(syn_T_C1.CD80_CTLA4_CD80)/dt
    NV_DATA_S(ydot)[92] = 1/PARAM(6)*(ReactionFlux130 - ReactionFlux131);

    //d(syn_T_C1.CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[93] = 1/PARAM(6)*(ReactionFlux128 - ReactionFlux129);

    //d(syn_T_C1.CD80_CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[94] = 1/PARAM(6)*(ReactionFlux129 + ReactionFlux131);

    //d(syn_T_C1.CD86_CTLA4)/dt
    NV_DATA_S(ydot)[95] = 1/PARAM(6)*(ReactionFlux132 - ReactionFlux133);

    //d(syn_T_C1.CD86_CTLA4_CD86)/dt
    NV_DATA_S(ydot)[96] = 1/PARAM(6)*(ReactionFlux133);

    //d(syn_T_C1.PDL1_CD80)/dt
    NV_DATA_S(ydot)[97] = 1/PARAM(6)*(ReactionFlux137 - ReactionFlux138 - ReactionFlux139);

    //d(syn_T_C1.PDL1_CD80_CD28)/dt
    NV_DATA_S(ydot)[98] = 1/PARAM(6)*(ReactionFlux138);

    //d(syn_T_C1.PDL1_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[99] = 1/PARAM(6)*(ReactionFlux139);

    //d(syn_T_C1.CD28)/dt
    NV_DATA_S(ydot)[100] = 1/PARAM(6)*( - ReactionFlux124 - ReactionFlux125 - ReactionFlux126 - ReactionFlux138);

    //d(syn_T_C1.CTLA4)/dt
    NV_DATA_S(ydot)[101] = 1/PARAM(6)*( - ReactionFlux127 - ReactionFlux128 - ReactionFlux131 - ReactionFlux132 - ReactionFlux134 - ReactionFlux135 - ReactionFlux139);

    //d(syn_T_C1.CD80)/dt
    NV_DATA_S(ydot)[102] = 1/PARAM(6)*( - ReactionFlux124 - ReactionFlux127 - ReactionFlux129 - ReactionFlux130 + ReactionFlux136);

    //d(syn_T_C1.CD80m)/dt
    NV_DATA_S(ydot)[103] = 1/PARAM(6)*( - ReactionFlux136 - ReactionFlux136 - ReactionFlux137);

    //d(syn_T_C1.CD86)/dt
    NV_DATA_S(ydot)[104] = 1/PARAM(6)*( - ReactionFlux126 - ReactionFlux132 - ReactionFlux133);

    //d(syn_T_C1.CTLA4_aCTLA4)/dt
    NV_DATA_S(ydot)[105] = 1/PARAM(6)*(ReactionFlux134 - ReactionFlux135);

    //d(syn_T_C1.CTLA4_aCTLA4_CTLA4)/dt
    NV_DATA_S(ydot)[106] = 1/PARAM(6)*(ReactionFlux135);

    //d(syn_T_APC.PD1_PDL1)/dt
    NV_DATA_S(ydot)[107] = 1/PARAM(7)*(ReactionFlux146);

    //d(syn_T_APC.PD1_PDL2)/dt
    NV_DATA_S(ydot)[108] = 1/PARAM(7)*(ReactionFlux147);

    //d(syn_T_APC.PD1)/dt
    NV_DATA_S(ydot)[109] = 1/PARAM(7)*( - ReactionFlux146 - ReactionFlux147 - ReactionFlux148 - ReactionFlux149);

    //d(syn_T_APC.PDL1)/dt
    NV_DATA_S(ydot)[110] = 1/PARAM(7)*(ReactionFlux142 + ReactionFlux144 - ReactionFlux146 - ReactionFlux150 - ReactionFlux151 - ReactionFlux165);

    //d(syn_T_APC.PDL2)/dt
    NV_DATA_S(ydot)[111] = 1/PARAM(7)*(ReactionFlux143 + ReactionFlux145 - ReactionFlux147);

    //d(syn_T_APC.PD1_aPD1)/dt
    NV_DATA_S(ydot)[112] = 1/PARAM(7)*(ReactionFlux148 - ReactionFlux149);

    //d(syn_T_APC.PD1_aPD1_PD1)/dt
    NV_DATA_S(ydot)[113] = 1/PARAM(7)*(ReactionFlux149);

    //d(syn_T_APC.PDL1_aPDL1)/dt
    NV_DATA_S(ydot)[114] = 1/PARAM(7)*(ReactionFlux150 - ReactionFlux151);

    //d(syn_T_APC.PDL1_aPDL1_PDL1)/dt
    NV_DATA_S(ydot)[115] = 1/PARAM(7)*(ReactionFlux151);

    //d(syn_T_APC.TPDL1)/dt
    NV_DATA_S(ydot)[116] = 1/PARAM(7)*( - ReactionFlux168 - ReactionFlux169);

    //d(syn_T_APC.TPDL1_aPDL1)/dt
    NV_DATA_S(ydot)[117] = 1/PARAM(7)*(ReactionFlux168 - ReactionFlux169);

    //d(syn_T_APC.TPDL1_aPDL1_TPDL1)/dt
    NV_DATA_S(ydot)[118] = 1/PARAM(7)*(ReactionFlux169);

    //d(syn_T_APC.CD28_CD80)/dt
    NV_DATA_S(ydot)[119] = 1/PARAM(7)*(ReactionFlux152 - ReactionFlux153);

    //d(syn_T_APC.CD28_CD80_CD28)/dt
    NV_DATA_S(ydot)[120] = 1/PARAM(7)*(ReactionFlux153);

    //d(syn_T_APC.CD28_CD86)/dt
    NV_DATA_S(ydot)[121] = 1/PARAM(7)*(ReactionFlux154);

    //d(syn_T_APC.CD80_CTLA4)/dt
    NV_DATA_S(ydot)[122] = 1/PARAM(7)*(ReactionFlux155 - ReactionFlux156 - ReactionFlux158);

    //d(syn_T_APC.CD80_CTLA4_CD80)/dt
    NV_DATA_S(ydot)[123] = 1/PARAM(7)*(ReactionFlux158 - ReactionFlux159);

    //d(syn_T_APC.CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[124] = 1/PARAM(7)*(ReactionFlux156 - ReactionFlux157);

    //d(syn_T_APC.CD80_CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[125] = 1/PARAM(7)*(ReactionFlux157 + ReactionFlux159);

    //d(syn_T_APC.CD86_CTLA4)/dt
    NV_DATA_S(ydot)[126] = 1/PARAM(7)*(ReactionFlux160 - ReactionFlux161);

    //d(syn_T_APC.CD86_CTLA4_CD86)/dt
    NV_DATA_S(ydot)[127] = 1/PARAM(7)*(ReactionFlux161);

    //d(syn_T_APC.PDL1_CD80)/dt
    NV_DATA_S(ydot)[128] = 1/PARAM(7)*(ReactionFlux165 - ReactionFlux166 - ReactionFlux167);

    //d(syn_T_APC.PDL1_CD80_CD28)/dt
    NV_DATA_S(ydot)[129] = 1/PARAM(7)*(ReactionFlux166);

    //d(syn_T_APC.PDL1_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[130] = 1/PARAM(7)*(ReactionFlux167);

    //d(syn_T_APC.CD28)/dt
    NV_DATA_S(ydot)[131] = 1/PARAM(7)*( - ReactionFlux152 - ReactionFlux153 - ReactionFlux154 - ReactionFlux166);

    //d(syn_T_APC.CTLA4)/dt
    NV_DATA_S(ydot)[132] = 1/PARAM(7)*( - ReactionFlux155 - ReactionFlux156 - ReactionFlux159 - ReactionFlux160 - ReactionFlux162 - ReactionFlux163 - ReactionFlux167);

    //d(syn_T_APC.CD80)/dt
    NV_DATA_S(ydot)[133] = 1/PARAM(7)*( - ReactionFlux152 - ReactionFlux155 - ReactionFlux157 - ReactionFlux158 + ReactionFlux164);

    //d(syn_T_APC.CD80m)/dt
    NV_DATA_S(ydot)[134] = 1/PARAM(7)*( - ReactionFlux164 - ReactionFlux164 - ReactionFlux165);

    //d(syn_T_APC.CD86)/dt
    NV_DATA_S(ydot)[135] = 1/PARAM(7)*( - ReactionFlux154 - ReactionFlux160 - ReactionFlux161);

    //d(syn_T_APC.CTLA4_aCTLA4)/dt
    NV_DATA_S(ydot)[136] = 1/PARAM(7)*(ReactionFlux162 - ReactionFlux163);

    //d(syn_T_APC.CTLA4_aCTLA4_CTLA4)/dt
    NV_DATA_S(ydot)[137] = 1/PARAM(7)*(ReactionFlux163);

    //d(syn_M_C.CD47)/dt
    NV_DATA_S(ydot)[138] = 1/PARAM(8)*( - ReactionFlux226);

    //d(syn_M_C.SIRPa)/dt
    NV_DATA_S(ydot)[139] = 1/PARAM(8)*( - ReactionFlux226);

    //d(syn_M_C.CD47_SIRPa)/dt
    NV_DATA_S(ydot)[140] = 1/PARAM(8)*(ReactionFlux226);

    //d(syn_M_C.PD1_PDL1)/dt
    NV_DATA_S(ydot)[141] = 1/PARAM(8)*(ReactionFlux227);

    //d(syn_M_C.PD1_PDL2)/dt
    NV_DATA_S(ydot)[142] = 1/PARAM(8)*(ReactionFlux228);

    //d(syn_M_C.PD1)/dt
    NV_DATA_S(ydot)[143] = 1/PARAM(8)*( - ReactionFlux227 - ReactionFlux228 - ReactionFlux229 - ReactionFlux230);

    //d(syn_M_C.PDL1)/dt
    NV_DATA_S(ydot)[144] = 1/PARAM(8)*(ReactionFlux222 + ReactionFlux224 - ReactionFlux227 - ReactionFlux231 - ReactionFlux232 - ReactionFlux234);

    //d(syn_M_C.PDL2)/dt
    NV_DATA_S(ydot)[145] = 1/PARAM(8)*(ReactionFlux223 + ReactionFlux225 - ReactionFlux228);

    //d(syn_M_C.PD1_aPD1)/dt
    NV_DATA_S(ydot)[146] = 1/PARAM(8)*(ReactionFlux229 - ReactionFlux230);

    //d(syn_M_C.PD1_aPD1_PD1)/dt
    NV_DATA_S(ydot)[147] = 1/PARAM(8)*(ReactionFlux230);

    //d(syn_M_C.PDL1_aPDL1)/dt
    NV_DATA_S(ydot)[148] = 1/PARAM(8)*(ReactionFlux231 - ReactionFlux232);

    //d(syn_M_C.PDL1_aPDL1_PDL1)/dt
    NV_DATA_S(ydot)[149] = 1/PARAM(8)*(ReactionFlux232);

    //d(syn_M_C.PDL1_CD80)/dt
    NV_DATA_S(ydot)[150] = 1/PARAM(8)*(ReactionFlux234);

    //d(syn_M_C.CD80)/dt
    NV_DATA_S(ydot)[151] = 1/PARAM(8)*(ReactionFlux233);

    //d(syn_M_C.CD80m)/dt
    NV_DATA_S(ydot)[152] = 1/PARAM(8)*( - ReactionFlux233 - ReactionFlux233 - ReactionFlux234);

    return(0);
}
int ODE_system::g(realtype t, N_Vector y, realtype *gout, void *user_data){

    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);

    //Assignment rules:

    //C_total < (0.5 * cell)
    gout[0] = 0.5 * PARAM(10) - - (SPVAR(23)) - (SPVAR(26));

    //V_T.C1 < (0.5 * cell)
    gout[1] = 0.5 * PARAM(10) - (SPVAR(23));

    //V_T.C2 < (0.5 * cell)
    gout[2] = 0.5 * PARAM(10) - (SPVAR(26));

    return(0);
}

bool ODE_system::triggerComponentEvaluate(int i, realtype t, bool curr) {

    bool discrete = false;
    realtype diff = 0;
    bool eval = false;
    //Assignment rules:

    switch(i)
    {
    case 0:
        //C_total < (0.5 * cell)
        diff = 0.5 * _class_parameter[10] - (NV_DATA_S(_y)[23]) - (NV_DATA_S(_y)[26]);
        break;
    case 1:
        //V_T.C1 < (0.5 * cell)
        diff = 0.5 * _class_parameter[10] - (NV_DATA_S(_y)[23]);
        break;
    case 2:
        //V_T.C2 < (0.5 * cell)
        diff = 0.5 * _class_parameter[10] - (NV_DATA_S(_y)[26]);
        break;
    default:
        break;
    }
    if (!discrete){
        eval = diff == 0 ? curr : (diff > 0);
    }
    return eval;
}

bool ODE_system::eventEvaluate(int i) {
    bool eval = false;
    switch(i)
    {
    case 0:
        eval = getSatisfied(0);
        break;
    case 1:
        eval = getSatisfied(1);
        break;
    case 2:
        eval = getSatisfied(2);
        break;
    default:
        break;
    }
    return eval;
}

bool ODE_system::eventExecution(int i, bool delayed, realtype& dt){

    bool setDelay = false;

    //Assignment rules:

    switch(i)
    {
    case 0:
        NV_DATA_S(_y)[24] = 0.01 * _class_parameter[10];
        break;
    case 1:
        NV_DATA_S(_y)[23] = 0.01 * _class_parameter[10];
        break;
    case 2:
        NV_DATA_S(_y)[26] = 0.01 * _class_parameter[10];
        break;
    default:
        break;
    }
    return setDelay;
}
void ODE_system::update_y_other(void){

    realtype AUX_VAR_PDL1_total = NV_DATA_S(_y)[144] + NV_DATA_S(_y)[141] + NV_DATA_S(_y)[148] + 2.0 * NV_DATA_S(_y)[149] + NV_DATA_S(_y)[150];

    realtype AUX_VAR_PDL2_total = NV_DATA_S(_y)[142] + NV_DATA_S(_y)[145];

    //syn_M_C.PDL1_total
    _species_other[0] = AUX_VAR_PDL1_total;

    //syn_M_C.PDL2_total
    _species_other[1] = AUX_VAR_PDL2_total;

    return;
}
std::string ODE_system::getHeader(){

    std::string s = "";
    s += ",V_C.nT0";
    s += ",V_C.T0";
    s += ",V_C.nT1";
    s += ",V_C.T1";
    s += ",V_C.aPD1";
    s += ",V_C.aPDL1";
    s += ",V_C.aCTLA4";
    s += ",V_C.Th";
    s += ",V_C.cabozantinib";
    s += ",V_C.A_site1";
    s += ",V_C.A_site2";
    s += ",V_P.nT0";
    s += ",V_P.T0";
    s += ",V_P.nT1";
    s += ",V_P.T1";
    s += ",V_P.aPD1";
    s += ",V_P.aPDL1";
    s += ",V_P.aCTLA4";
    s += ",V_P.Th";
    s += ",V_P.cabozantinib";
    s += ",V_T.C_x";
    s += ",V_T.T1_exh";
    s += ",V_T.Th_exh";
    s += ",V_T.C1";
    s += ",V_T.K";
    s += ",V_T.c_vas";
    s += ",V_T.C2";
    s += ",V_T.T0";
    s += ",V_T.T1";
    s += ",V_T.IFNg";
    s += ",V_T.APC";
    s += ",V_T.mAPC";
    s += ",V_T.P0";
    s += ",V_T.P1";
    s += ",V_T.aPD1";
    s += ",V_T.aPDL1";
    s += ",V_T.aCTLA4";
    s += ",V_T.Th";
    s += ",V_T.TGFb";
    s += ",V_T.MDSC";
    s += ",V_T.NO";
    s += ",V_T.ArgI";
    s += ",V_T.CCL2";
    s += ",V_T.cabozantinib";
    s += ",V_T.Mac_M1";
    s += ",V_T.Mac_M2";
    s += ",V_T.IL12";
    s += ",V_T.IL10";
    s += ",V_T.Fib";
    s += ",V_T.CAF";
    s += ",V_T.ECM";
    s += ",V_LN.nT0";
    s += ",V_LN.aT0";
    s += ",V_LN.T0";
    s += ",V_LN.IL2";
    s += ",V_LN.nT1";
    s += ",V_LN.aT1";
    s += ",V_LN.T1";
    s += ",V_LN.APC";
    s += ",V_LN.mAPC";
    s += ",V_LN.aPD1";
    s += ",V_LN.aPDL1";
    s += ",V_LN.aCTLA4";
    s += ",V_LN.aTh";
    s += ",V_LN.Th";
    s += ",V_LN.cabozantinib";
    s += ",V_e.P0";
    s += ",V_e.p0";
    s += ",V_e.P1";
    s += ",V_e.p1";
    s += ",A_e.M1";
    s += ",A_e.M1p0";
    s += ",A_e.M1p1";
    s += ",A_s.M1";
    s += ",A_s.M1p0";
    s += ",A_s.M1p1";
    s += ",syn_T_C1.PD1_PDL1";
    s += ",syn_T_C1.PD1_PDL2";
    s += ",syn_T_C1.PD1";
    s += ",syn_T_C1.PDL1";
    s += ",syn_T_C1.PDL2";
    s += ",syn_T_C1.PD1_aPD1";
    s += ",syn_T_C1.PD1_aPD1_PD1";
    s += ",syn_T_C1.PDL1_aPDL1";
    s += ",syn_T_C1.PDL1_aPDL1_PDL1";
    s += ",syn_T_C1.TPDL1";
    s += ",syn_T_C1.TPDL1_aPDL1";
    s += ",syn_T_C1.TPDL1_aPDL1_TPDL1";
    s += ",syn_T_C1.CD28_CD80";
    s += ",syn_T_C1.CD28_CD80_CD28";
    s += ",syn_T_C1.CD28_CD86";
    s += ",syn_T_C1.CD80_CTLA4";
    s += ",syn_T_C1.CD80_CTLA4_CD80";
    s += ",syn_T_C1.CTLA4_CD80_CTLA4";
    s += ",syn_T_C1.CD80_CTLA4_CD80_CTLA4";
    s += ",syn_T_C1.CD86_CTLA4";
    s += ",syn_T_C1.CD86_CTLA4_CD86";
    s += ",syn_T_C1.PDL1_CD80";
    s += ",syn_T_C1.PDL1_CD80_CD28";
    s += ",syn_T_C1.PDL1_CD80_CTLA4";
    s += ",syn_T_C1.CD28";
    s += ",syn_T_C1.CTLA4";
    s += ",syn_T_C1.CD80";
    s += ",syn_T_C1.CD80m";
    s += ",syn_T_C1.CD86";
    s += ",syn_T_C1.CTLA4_aCTLA4";
    s += ",syn_T_C1.CTLA4_aCTLA4_CTLA4";
    s += ",syn_T_APC.PD1_PDL1";
    s += ",syn_T_APC.PD1_PDL2";
    s += ",syn_T_APC.PD1";
    s += ",syn_T_APC.PDL1";
    s += ",syn_T_APC.PDL2";
    s += ",syn_T_APC.PD1_aPD1";
    s += ",syn_T_APC.PD1_aPD1_PD1";
    s += ",syn_T_APC.PDL1_aPDL1";
    s += ",syn_T_APC.PDL1_aPDL1_PDL1";
    s += ",syn_T_APC.TPDL1";
    s += ",syn_T_APC.TPDL1_aPDL1";
    s += ",syn_T_APC.TPDL1_aPDL1_TPDL1";
    s += ",syn_T_APC.CD28_CD80";
    s += ",syn_T_APC.CD28_CD80_CD28";
    s += ",syn_T_APC.CD28_CD86";
    s += ",syn_T_APC.CD80_CTLA4";
    s += ",syn_T_APC.CD80_CTLA4_CD80";
    s += ",syn_T_APC.CTLA4_CD80_CTLA4";
    s += ",syn_T_APC.CD80_CTLA4_CD80_CTLA4";
    s += ",syn_T_APC.CD86_CTLA4";
    s += ",syn_T_APC.CD86_CTLA4_CD86";
    s += ",syn_T_APC.PDL1_CD80";
    s += ",syn_T_APC.PDL1_CD80_CD28";
    s += ",syn_T_APC.PDL1_CD80_CTLA4";
    s += ",syn_T_APC.CD28";
    s += ",syn_T_APC.CTLA4";
    s += ",syn_T_APC.CD80";
    s += ",syn_T_APC.CD80m";
    s += ",syn_T_APC.CD86";
    s += ",syn_T_APC.CTLA4_aCTLA4";
    s += ",syn_T_APC.CTLA4_aCTLA4_CTLA4";
    s += ",syn_M_C.CD47";
    s += ",syn_M_C.SIRPa";
    s += ",syn_M_C.CD47_SIRPa";
    s += ",syn_M_C.PD1_PDL1";
    s += ",syn_M_C.PD1_PDL2";
    s += ",syn_M_C.PD1";
    s += ",syn_M_C.PDL1";
    s += ",syn_M_C.PDL2";
    s += ",syn_M_C.PD1_aPD1";
    s += ",syn_M_C.PD1_aPD1_PD1";
    s += ",syn_M_C.PDL1_aPDL1";
    s += ",syn_M_C.PDL1_aPDL1_PDL1";
    s += ",syn_M_C.PDL1_CD80";
    s += ",syn_M_C.CD80";
    s += ",syn_M_C.CD80m";
    s += ",syn_M_C.PDL1_total";
    s += ",syn_M_C.PDL2_total";
    return s;
}
realtype ODE_system::get_unit_conversion_species(int i) const{

    static std::vector<realtype> scalor = {
        //sp_var
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.0000000000000002e-06,
        1.0000000000000002e-06,
        1.0000000000000002e-06,
        1.66053872801495e-24,
        1.0000000000000002e-06,
        1.0000000000000002e-06,
        1.0000000000000002e-06,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.0000000000000002e-06,
        1.0000000000000002e-06,
        1.0000000000000002e-06,
        1.66053872801495e-24,
        1.0000000000000002e-06,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1e-09,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1e-06,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000.0,
        1000.0,
        1e-06,
        1e-06,
        1e-06,
        1.66053872801495e-24,
        1e-06,
        1.66053872801495e-24,
        1e-06,
        1e-06,
        1e-06,
        1e-06,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1e-06,
        1e-06,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1e-09,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1e-06,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1e-06,
        1e-06,
        1e-06,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1e-06,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        //sp_other
        1.66053872801495e-12,
        1.66053872801495e-12,
    };
    return scalor[i];
}
realtype ODE_system::get_unit_conversion_nspvar(int i) const{

    static std::vector<realtype> scalor = {
    };
    return scalor[i];
}
};
