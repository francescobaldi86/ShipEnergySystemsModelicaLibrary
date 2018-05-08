within ShipEnergySystemsModelicaLibrary;

package HeatExchanger
  model Evaporator "Counter current Moving Boundary model: Fluid enters subcooled and exits in super-heated conditons. The model consider the fluid in one side and the metal wall. The secondary fluid is a constant specific heat fluid"
    // Note that this model is taken "as is" from the ThermoCycle Library //
    /**************** MEDIUM ***************************/
    replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialMedium annotation(
       choicesAllMatching = true);
    /***************** PORTS ******************/
    Modelica.Fluid.Interfaces.FluidPort_a InFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-108, -74}, {-88, -54}})));
    Modelica.Fluid.Interfaces.FluidPort_b OutFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{90, -74}, {110, -54}})));
    ShipEnergySystemsModelicaLibrary.Other.FlangeCdot_B OutFlow_sf annotation(
      Placement(transformation(extent = {{-110, 50}, {-90, 70}})));
    ShipEnergySystemsModelicaLibrary.Other.FlangeCdot_A InFlow_sf annotation(
      Placement(transformation(extent = {{90, 50}, {110, 70}})));
    /***********THERMODYNAMIC STATES *********************/
    Medium.ThermodynamicState SubCooled "Subcooled thermodynamic state";
    Medium.SaturationProperties sat;
    Medium.ThermodynamicState SuperHeated "Superheated thermodynamic state";
    Medium.ThermodynamicState outlet "Outlet thermodynamic state";
    /************* Tube pressure *****************/
    Modelica.SIunits.Pressure p(start = 100000) "Pressure in the heat-exchanger: No pressure drop considered";
    /*********** GENERAL VARIABLES AND PARAMETER ***********/
    parameter Modelica.SIunits.Area A "Cross-sectional area for both side of the heat exchanger";
    parameter Modelica.SIunits.Length L "Total length of the exchanger";
    final parameter Modelica.SIunits.Length D = 2 * sqrt(A / pi) "Diameter";
    constant Real pi = Modelica.Constants.pi;
    /*************** Wall between the fluids *************************/
    parameter Modelica.SIunits.Mass M_tot "Total mass of metal walll between the fluids";
    parameter Modelica.SIunits.SpecificHeatCapacity c_wall "Heat capacity of the wall - assumed constant in the whole heat exchanger";
    parameter Modelica.SIunits.Density rho_wall "Density of the wall - assumed constant in the whole heat exchanger";
    /********** MASS FLOWS *******************/
    //parameter Modelica.SIunits.MassFlowRate Mdotnom
    //    "Nominal fluid flow rate of the primary fluid";
    Modelica.SIunits.MassFlowRate M_dot_su "Mass flow at the inlet of the heat exchanger";
    Modelica.SIunits.MassFlowRate M_dot_A "Mass flow at the interface between sub-cooled and two-phase";
    Modelica.SIunits.MassFlowRate M_dot_B "Mass flow at the interface between sub-cooled and two-phase";
    Modelica.SIunits.MassFlowRate M_dot_ex "Mass flow at the outlet of the heat exchanger";
    Modelica.SIunits.MassFlowRate Mdot_sf "Mass flow of the secondary fluid of the heat exchanger";
    /*****************HEAT FLOW *************************/
    Modelica.SIunits.Power Q_SB "Thermal energy transfer with the wall";
    Modelica.SIunits.Power Q_TP "Thermal energy transfer with the wall in the two-phase region";
    Modelica.SIunits.Power Q_SH "Thermal energy transfer with the wall in the superheated region";
    Modelica.SIunits.Power Q_tot "Total Thermal energy transfer with the wall from the primary fluid";
    Modelica.SIunits.Power Q_tot_sf "Total Thermal energy transfer with the wall from the secondary fluid";
    Modelica.SIunits.Power Qsf_SH;
    Modelica.SIunits.Power Qsf_TP;
    Modelica.SIunits.Power Qsf_SB;
    /***************** TEMPERATURES *************************/
    /** primary fluid **/
    Modelica.SIunits.Temperature T_SU "Temperatue at the inlet of the heat exchanger";
    Modelica.SIunits.Temperature T_SB "Mean Temperatue of the subcooled region";
    Modelica.SIunits.Temperature T_TP "Mean Temperatue of the two-phase region";
    Modelica.SIunits.Temperature T_SH "Mean Temperatue of the superheated region";
    Modelica.SIunits.Temperature TwA "Wall Temperature at the interface between sub-cooled and two-phase";
    Modelica.SIunits.Temperature TwB "Wall Temperature at the interface between two-phase and super-heated";
    Modelica.SIunits.Temperature T_out "Temperature at the outlet of the moving boundary model";
    /** Wall **/
    Modelica.SIunits.Temperature TwSB(start = 300, fixed = true) "Mean Wall temperature in the sub-cooled region";
    Modelica.SIunits.Temperature TwTP(start = 300, fixed = true) "Mean Wall Temperature in the two-phase region";
    Modelica.SIunits.Temperature TwSH(start = 300, fixed = true) "Mean Wall Temperature in the super-heated region";
    /** Secondary fluid **/
    parameter Modelica.SIunits.Temperature dTsf_start = 10 "|Initialization|Temperature change";
    parameter Modelica.SIunits.Temperature Tsf_SU_start = 273.15 + 25 "|Initialization|Temperature inlet";
    Modelica.SIunits.Temperature Tsf_SU(start = Tsf_SU_start, fixed = true) "Tempearture at the inlet of the secondary fluid";
    Modelica.SIunits.Temperature Tsf_SH(start = Tsf_SU_start - 1 * dTsf_start / 6, fixed = true) "Mean Temperatue of the superheated region secondary fluid";
    Modelica.SIunits.Temperature Tsf_B(start = Tsf_SU_start - 2 * dTsf_start / 6, fixed = true) "Secondary fluid temperature at the interface between two-phase and super-heated";
    Modelica.SIunits.Temperature Tsf_TP(start = Tsf_SU_start - 3 * dTsf_start / 6, fixed = true) "Mean Temperatue of the two-phase region secondary fluid";
    Modelica.SIunits.Temperature Tsf_A(start = Tsf_SU_start - 4 * dTsf_start / 6, fixed = true) "Secondary fluid temperature at the interface between two-phase and super-heated";
    Modelica.SIunits.Temperature Tsf_SB(start = Tsf_SU_start - 5 * dTsf_start / 6, fixed = true) "Mean Temperatue of the subcooled region secondary fluid";
    Modelica.SIunits.Temperature Tsf_EX(start = Tsf_SU_start - dTsf_start, fixed = true) "Secondary fluid temperature at the exit of the heat exchanger";
    /*************** HEAT TRANSFER PARAMETER ****************/
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_SB "Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_TP "Heat transfer coefficient two-phase side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_SH "Heat transfer coefficient Superheated side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer Usf "Secondary fluid Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.QualityFactor ETA = 1 "Fin efficiency";
    /*************Tube region length ********************/
    Modelica.SIunits.Length L_SB(start = 0.8 * L, min = 0) "Length of the subcooled region";
    Modelica.SIunits.Length L_TP(start = 0.1 * L, min = 0) "Length of the two-phase region";
    Modelica.SIunits.Length L_SH(start = 0.1 * L, min = 0) "Length of the Superheated region";
    /********** Densities ***************/
    Modelica.SIunits.Density rho_SB "Mean density of the subcooled region";
    Modelica.SIunits.Density rho_LS "Bubble density";
    Modelica.SIunits.Density rho_VS "Dew density";
    Modelica.SIunits.Density rho_TP "Mean density of the two-phase region";
    Modelica.SIunits.Density rho_SH "Mean density of the superheated region";
    Modelica.SIunits.Density rho_sf "Density of the secondary fluid at the inlet";
    /********** Specific heat capacity ***************/
    Modelica.SIunits.SpecificHeatCapacity cp_sf "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_SB "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_TP "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_SH "Specific heat capacity of the secondary fluid at the inlet";
    /************** Enthalpies ************************/
    Modelica.SIunits.SpecificEnthalpy h_SB "Mean Enthalpy of the subcooled region";
    Modelica.SIunits.SpecificEnthalpy h_LS "Bubble enthalpy";
    Modelica.SIunits.SpecificEnthalpy h_VS "Dew enthalpy";
    //Modelica.SIunits.SpecificEnthalpy h_TP "Mean Enthalpy of the two-phase region";
    Modelica.SIunits.SpecificEnthalpy h_SH "Mean Enthalpy of the Superheated region";
    Modelica.SIunits.SpecificEnthalpy h_SU "Enthalpy at the inlet of the heat exchanger of the primary fluid";
    Modelica.SIunits.SpecificEnthalpy h_EX(start = 0) "Enthalpy at the outlet of the heat exchanger of the primary fluid";
    /***** Enthalpy&Density****/
    Real rhoh_TP "Density times enthalpy in the two-phase region";
    /******************* Density derivatives ********************/
    Modelica.SIunits.DerDensityByEnthalpy drdh_SB "Density derivative with respect to enthalpy at constant pressure in the subcooled region";
    Modelica.SIunits.DerDensityByPressure drdp_SB "Density derivative with respect to pressure at constant enthalpy in the subcooled region";
    Modelica.SIunits.DerDensityByPressure drdp_LS "Bubble point density derivative with respect to pressure";
    Modelica.SIunits.DerDensityByPressure drdp_VS "Dew point density derivative with respect to pressure";
    Modelica.SIunits.DerDensityByEnthalpy drdh_SH "Density derivative with respect to enthalpy at constant pressure in the superheated region";
    Modelica.SIunits.DerDensityByPressure drdp_SH "Density derivative with respect to pressure at constant enthalpy in the superheated region";
    /*************** Enthalpy derivatives *******************/
    Modelica.SIunits.DerEnthalpyByPressure dhdp_LS "Bubble point specific enthalpy derivative with respect to pressure";
    Modelica.SIunits.DerEnthalpyByPressure dhdp_VS "Dew point specific enthalpy derivative with respect to pressure";
    /************ Derivative with respect to time **************/
    Real drdt_SB "Density derivative with respect to time in the subcooled region";
    Real dhdt_SB "Enthalpy derivative with respect to time in the subcooled region";
    Real drdt_TP "Density derivative with respect to time in the two-phase region";
    Real drhdt_TP "Density and Enthalpy derivative with respect to time in the two-phase region";
    Real drdt_SH "Density derivative with respect to time in the superheated region";
    Real dhdt_SH "Enthalpy derivative with respect to time in the superheated region";
    /********************** VOID fraction *****************************/
    parameter Real Void "Void fraction assumed constant for now";
    parameter Real dVoid_dp = 0 "Void fraction derivative with respect to pressure";
    parameter Real dVoid_dh = 0 "Void fraction derivative with respect to enthalpy";
    /********************** Epsilon-NTU method *********************/
    /******* Secondary fluid - Wall ****************/
    Real epsilon_SB_sf "Epsilon Sub-cooled secondary fluid - wall Side";
    Real NTU_SB_sf "NTU Sub-cooled secondary fluid-Wall side";
    Real epsilon_TP_sf "Epsilon Two-Phase secondary fluid - wall Side";
    Real NTU_TP_sf "NTU Two-Phase secondary fluid-Wall side";
    Real epsilon_SH_sf "Epsilon Super-heated secondary fluid - wall Side";
    Real NTU_SH_sf "NTU Super-heated secondary fluid-Wall side";
    /*********** Wall - Primary fluid ************/
    Real epsilon_SB_pf "Epsilon Sub-cooled wall side - primary fluid";
    Real NTU_SB_pf "NTU Sub-cooled wall side - primary fluid";
    Real epsilon_TP_pf "Epsilon Two-Phase wall side - primary fluid";
    Real NTU_TP_pf "NTU Two-Phase wall side - primary fluid";
    //
    //
    Real epsilon_SH_pf "Epsilon Super-heated wall side - primary fluid";
    Real NTU_SH_pf "NTU Super-heated secondary wall side - primary fluid";
  equation
/****************** BOUNDARY EQUATION *************************/
    L = L_SB + L_TP + L_SH;
/********** Enthalpies ************/
    h_SU = inStream(InFlow.h_outflow);
    h_EX = OutFlow.h_outflow;
    h_SU = InFlow.h_outflow;
/*********** Temperatures *********/
    OutFlow_sf.T = Tsf_EX;
    InFlow_sf.T = Tsf_SU;
/**********Specific Heat Capacity *************/
    OutFlow_sf.cp = cp_sf;
    InFlow_sf.cp = cp_sf;
/*********** Density ***********/
    OutFlow_sf.rho = rho_sf;
    InFlow_sf.rho = rho_sf;
/*************** Pressure **************/
    InFlow.p = p;
    OutFlow.p = p;
/***********Mass Flow ******************/
    InFlow.m_flow = M_dot_su;
    OutFlow.m_flow = -M_dot_ex;
    OutFlow_sf.Mdot = Mdot_sf;
    InFlow_sf.Mdot = Mdot_sf;
/*************** CONSTITUTIVE EQUATIONS ******************/
    h_SB = (h_SU + h_LS) / 2 "Mean enthalpy of the subcooled region";
    rho_TP = rho_VS * Void + rho_LS * (1 - Void);
    rhoh_TP = rho_VS * h_VS * Void + rho_LS * h_LS * (1 - Void);
    h_SH = (h_EX + h_VS) / 2 "Mean enthalpy of the superheated region";
    TwA = (TwSB * L_TP + TwTP * L_SB) / (L_TP + L_SB) "Wall Temperature at the border between subcooled and two-phase region";
    TwB = (TwSH * L_TP + TwTP * L_SH) / (L_TP + L_SH) "Wall Temperature at the border between two-phase and superheated region";
    Tsf_SH = (Tsf_SU + Tsf_B) / 2 "Mean temperature of the secondary fluid in the superheated region";
// Tsf_B = (L_TP*Tsf_SH + L_SH*Tsf_TP)/(L_SH+L_TP)
//     "Secondary fluid temperature at the border between two-phase and superheated region";
    Tsf_TP = (Tsf_A + Tsf_B) / 2 "Mean temperature of the secondary fluid in the two phase region";
// Tsf_A = (L_TP*Tsf_SB + L_SB*Tsf_TP)/(L_SB+L_TP)
//     "Secondary fluid temperature at the border between two-phase and subcooled";
    Tsf_SB = (Tsf_A + Tsf_EX) / 2 "Mean temperature of the secondary fluid in the two phase region";
/**** SECONDARY FLUID - WALL HEAT TRANSFER ****/
/** Subcooled zone **/
    NTU_SB_sf = Usf * pi * D * L_SB / (cp_sf * Mdot_sf);
    epsilon_SB_sf = 1 - exp(-NTU_SB_sf);
    Qsf_SB = cp_sf * Mdot_sf * epsilon_SB_sf * (Tsf_A - TwSB);
/** Two-phase zone **/
    NTU_TP_sf = Usf * pi * D * L_TP / (cp_sf * Mdot_sf);
    epsilon_TP_sf = 1 - exp(-NTU_TP_sf);
    Qsf_TP = cp_sf * Mdot_sf * epsilon_TP_sf * (Tsf_B - TwTP);
/** SuperHeated zone **/
    NTU_SH_sf = Usf * pi * D * L_SH / (cp_sf * Mdot_sf);
    epsilon_SH_sf = 1 - exp(-NTU_SH_sf);
    Qsf_SH = cp_sf * Mdot_sf * epsilon_SH_sf * (Tsf_SU - TwSH);
    Qsf_SB = Mdot_sf * cp_sf * (Tsf_A - Tsf_EX);
    Qsf_TP = Mdot_sf * cp_sf * (Tsf_B - Tsf_A);
    Qsf_SH = Mdot_sf * cp_sf * (Tsf_SU - Tsf_B);
/*** WALL - PRIMARY FLUID HEAT TRANSFER ***/
/** Sub-cooled zone **/
    NTU_SB_pf = U_SB * pi * D * L_SB / (cp_SB * (M_dot_su + M_dot_A) / 2);
    epsilon_SB_pf = 1 - exp(-NTU_SB_pf);
    Q_SB = cp_SB * (M_dot_su + M_dot_A) / 2 * epsilon_SB_pf * (TwSB - T_SU) "Heat from the wall sub-cooled region";
/** Two-phase zone **/
/** I say that Epsilon NTU in the TP zone is useless **/
    NTU_TP_pf = U_TP * pi * D * L_TP / (cp_TP * (M_dot_A + M_dot_B) / 2);
    epsilon_TP_pf = 1 - exp(-NTU_TP_pf);
    Q_TP = pi * D * L_TP * U_TP * (TwTP - T_TP) "Heat from the wall two-phase region";
/** SuperHeated zone **/
    NTU_SH_pf = U_SH * pi * D * L_SH / (cp_SH * (M_dot_B + M_dot_ex) / 2);
    epsilon_SH_pf = 1 - exp(-NTU_SH_pf);
    Q_SH = cp_SH * (M_dot_B + M_dot_ex) / 2 * (TwSH - T_TP) "Heat from the wall super-heated region";
/** Total energy exchange from the primary fluid and the secondary fluid**/
    Q_tot = Q_SB + Q_TP + Q_SH;
    Q_tot_sf = Qsf_SB + Qsf_TP + Qsf_SH;
/************ FLUID PROPERTIES ************************/
    SubCooled = Medium.setState_ph(p, h_SB);
    sat = Medium.setSat_p(p);
    SuperHeated = Medium.setState_ph(p, h_SH);
    outlet = Medium.setState_ph(p, h_EX);
/**Temperatures **/
    T_SU = Medium.temperature_ph(p, h_SU);
    T_SB = Medium.temperature(SubCooled);
    T_TP = Medium.saturationTemperature_sat(sat);
    T_SH = Medium.temperature(SuperHeated);
    T_out = Medium.temperature(outlet);
/** Densities **/
    rho_SB = Medium.density(SubCooled);
    rho_LS = Medium.bubbleDensity(sat);
    rho_VS = Medium.dewDensity(sat);
    rho_SH = Medium.density(SuperHeated);
/** SpecificHeatCapacity **/
    cp_SB = Medium.specificHeatCapacityCp(SubCooled);
    cp_TP = (cp_SB + cp_SH) / 2;
    cp_SH = Medium.specificHeatCapacityCp(SuperHeated);
/** Enthalpies **/
    h_LS = Medium.bubbleEnthalpy(sat);
    h_VS = Medium.dewEnthalpy(sat);
/**Derivatives**/
    drdh_SB = Medium.density_derh_p(SubCooled);
    drdp_SB = Medium.density_derp_h(SubCooled);
    drdp_LS = Medium.dBubbleDensity_dPressure(sat);
    drdp_VS = Medium.dDewDensity_dPressure(sat);
    dhdp_LS = Medium.dBubbleEnthalpy_dPressure(sat);
    dhdp_VS = Medium.dDewEnthalpy_dPressure(sat);
    drdh_SH = Medium.density_derh_p(SuperHeated);
    drdp_SH = Medium.density_derp_h(SuperHeated);
/************ MassBalance SUB-COOLED ****************/
    A * (L_SB * drdt_SB + (rho_SB - rho_LS) * der(L_SB)) = M_dot_su - M_dot_A;
    drdt_SB = (drdp_SB + 1 / 2 * drdh_SB * dhdp_LS) * der(p) + 1 / 2 * drdh_SB * der(h_SU);
/********** EnergyBalance SUB-COOLED ******************/
    A * L_SB * (rho_SB * dhdt_SB + h_SB * drdt_SB - der(p)) + A * (rho_SB * h_SB - rho_LS * h_LS) * der(L_SB) = M_dot_su * h_SU - M_dot_A * h_LS + Q_SB;
    dhdt_SB = 1 / 2 * (dhdp_LS * der(p) + der(h_SU));
/*******  WallEnergyBalance SUB-COOLED **********/
    c_wall * M_tot / L * (L_SB * der(TwSB) + (TwSB - TwA) * der(L_SB)) = Qsf_SB - Q_SB;
/**************************************************************************************************************/
/************ MassBalance TWO-PHASE REGION ******************/
    A * (L_TP * drdt_TP + (rho_TP - rho_VS) * der(L_TP) + (rho_LS - rho_VS) * der(L_SB)) = M_dot_A - M_dot_B;
    drdt_TP = (Void * drdp_VS + (1 - Void) * drdp_LS) * der(p);
//(rho_VS - rho_LS)*(dVoid_dp*der(p) + dVoid_dh*der(h_EX));
/************ EnergyBalance TWO-PHASE REGION ******************/
    A * (L_TP * drhdt_TP + (rhoh_TP - rho_VS * h_VS) * der(L_TP) + (rho_LS * h_LS - rho_VS * h_VS) * der(L_SB) - L_TP * der(p)) = M_dot_A * h_LS - M_dot_B * h_VS + Q_TP;
    drhdt_TP = (Void * (drdp_VS * h_VS + rho_VS * drdp_VS) + (1 - Void) * (h_LS * drdp_LS + rho_LS * dhdp_LS)) * der(p);
//(rho_VS*h_VS -rho_LS*h_LS)*(dVoid_dp*der(p) + dVoid_dh*der(h_EX));
/***************** WallEnergyBalance TWO-PHASE REGION **************/
    c_wall * M_tot / L * (L_TP * der(TwTP) + (TwA - TwB) * der(L_SB) + (TwTP - TwB) * der(L_TP)) = Qsf_TP - Q_TP;
/**************************************************************************************************************/
/************ MassBalance SUPER-HEATED ****************/
    A * (L_SH * drdt_SH + (rho_VS - rho_SH) * der(L_SH)) = M_dot_B - M_dot_ex;
    drdt_SH = (drdp_SH + 1 / 2 * drdh_SH * dhdp_VS) * der(p) + 1 / 2 * drdh_SH * der(h_EX);
/********** EnergyBalance SUPER-HEATED ******************/
    A * L_SH * (rho_SH * dhdt_SH + h_SH * drdt_SH - der(p)) + A * (rho_VS * h_VS - rho_SH * h_SH) * der(L_SH) = M_dot_B * h_VS - M_dot_ex * h_EX + Q_SH;
    dhdt_SH = 1 / 2 * (dhdp_VS * der(p) + der(h_EX));
/*******  WallEnergyBalance SUPER-HEATED **********/
    c_wall * M_tot / L * (L_SH * der(TwSH) + (TwB - TwSH) * der(L_SH)) = Qsf_SH - Q_SH;
/******************************* SUMMARY ***********************************/
  protected
    parameter Integer N = 3;
    Modelica.SIunits.Temperature[N] T_fluid_;
    Modelica.SIunits.Temperature[N] T_wall_;
    Modelica.SIunits.Temperature[N] T_sf_;
  public
    record SummaryBase
      replaceable Arrays T_profile;

      record Arrays
        parameter Integer n;
        Modelica.SIunits.Temperature[n] Twf;
        Modelica.SIunits.Temperature[n] Twall;
        Modelica.SIunits.Temperature[n] Tsf;
      end Arrays;

      //Modelica.SIunits.Pressure p_sf;
      Modelica.SIunits.Pressure p_wf;
      //Modelica.SIunits.Power Q_sf;
      //Modelica.SIunits.Power Q_wf;
    end SummaryBase;

    replaceable record SummaryClass = SummaryBase;
    SummaryClass Summary(T_profile(n = N, Twf = T_fluid_[:], Twall = T_wall_[:], Tsf = T_sf_[:]), p_wf = p);
  equation
/* Fluid temperature */
    T_fluid_[1] = T_SU;
    T_fluid_[2] = T_TP;
    T_fluid_[3] = T_out;
/* Wall temperature */
    T_wall_[1] = TwSB;
    T_wall_[2] = TwTP;
    T_wall_[3] = TwSH;
/* Secondary fluid temperature */
    T_sf_[1] = Tsf_SU;
    T_sf_[2] = Tsf_TP;
    T_sf_[3] = Tsf_EX;
    annotation(
      Diagram(graphics),
      Icon(graphics = {Rectangle(extent = {{-100, -20}, {-32, -100}}, lineColor = {60, 121, 182}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-100, 20}, {100, -20}}, lineColor = {135, 135, 135}, fillPattern = FillPattern.CrossDiag), Rectangle(extent = {{-100, 100}, {100, 20}}, lineColor = {135, 135, 135}, fillColor = {255, 85, 85}, fillPattern = FillPattern.Solid), Line(points = {{-100, 56}, {-80, 36}, {-60, 56}, {-40, 36}, {-20, 56}, {0, 36}, {20, 56}, {40, 36}, {60, 56}, {80, 36}, {100, 56}}, color = {255, 0, 0}, smooth = Smooth.None, thickness = 0.5), Polygon(points = {{-18, 84}, {-18, 64}, {-32, 74}, {-18, 84}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-18, 74}, {34, 74}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{-18, -72}, {-10, -80}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{-100, -100}, {40, -100}, {-16, -20}, {-100, -20}, {-100, -24}, {-100, -100}}, lineColor = {0, 128, 255}, smooth = Smooth.None, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{32, -100}, {44, -100}, {102, -100}, {100, -20}, {-16, -20}, {32, -100}}, lineColor = {170, 255, 255}, smooth = Smooth.None, fillColor = {170, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{-97, -63}, {-75, -43}, {-57, -63}, {-37, -43}, {-17, -63}, {3, -43}, {23, -63}, {43, -43}, {63, -63}, {83, -43}, {103, -63}}, color = {0, 0, 255}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{4, -30}, {12, -38}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-2, -54}, {6, -62}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{42, -70}, {42, -90}, {56, -80}, {42, -70}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-8, -80}, {44, -80}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{28, -28}, {34, -34}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-14, -26}, {-2, -40}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{12, -64}, {22, -76}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{26, -78}, {38, -90}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{16, -42}, {28, -54}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{36, -62}, {42, -68}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid)}),
      Documentation(info = "<HTML>
          
         <p><big>Model <b>MBeva</b> is a moving boundaries heat exchanger model for evaporators.
         The model considers a fictitius heat transfer channel split up into three different section based on the phase state of the working fluid:
         <ul><li> Sub-cooled zone (SB)
         <li> Two-phase zone (TP)
         <li> SuperHeated zone (SH)
         </ul>
          <p><big><b>Pressure</b> and <b>enthalpy</b> are selected as state variables for each of the three zones. 
          
        <p><big>  The name moving boundary is derived from the fact that the interfaces between these sections do not
have a fixed spatial position but merely a fixed thermodynamic location depending
on the presence of liquid and gaseous fluid, respectively. The actual existence of
a certain section and its length are determined based on the fluid state resulting
in variable sectioning. A fixed total length superimposes the required boundary
condition to calculate the length of each section.
 
 <p><big>The assumptions for this model are:
         <ul><li> The tube is cylindrical with a constant cross sectional area
         <li> The velocity of the fluid is uniform on the cross sectional area
         <li> The enthalpy of the fluid is linear in each region of the tube (sub-cooled, two-phase, super-heated)
         <li> Pressure is considered constant
         <li> Constant void fraction
         <li> Thermal energy accumulation in the metal wall is taken into account
         <li> The secondary fluid is treated as a constant heat capacity fluid
         <li> The heat transfer between the secondary fluid and the metal wall and between the metal wall and the primary fluid is computed using the epsilon-NTU method
         </ul>
        </HTML>"));
  end Evaporator;

  model HEX_1f "Counter current simple HEX model. No fluid undergo any phase change. The primary fluid is a proper Modelica fluid, and should be in liquid state (hence the ''secondary fluid - to - liquid'' denomination)"
    /**************** MEDIUM ***************************/
    replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialMedium annotation(
       choicesAllMatching = true);
    /***************** PORTS ******************/
    Modelica.Fluid.Interfaces.FluidPort_a InFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-108, -74}, {-88, -54}})));
    Modelica.Fluid.Interfaces.FluidPort_b OutFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{90, -74}, {110, -54}})));
    ShipEnergySystems.Other.FlangeCdot_B OutFlow_sf annotation(
      Placement(transformation(extent = {{-110, 50}, {-90, 70}})));
    ShipEnergySystems.Other.FlangeCdot_A InFlow_sf annotation(
      Placement(transformation(extent = {{90, 50}, {110, 70}})));
    /***********THERMODYNAMIC STATES *********************/
    Medium.ThermodynamicState PrimaryFluidState "Primary fluid thermodynamic state";
    Medium.ThermodynamicState inlet "Primary fluid thermodynamic state";
    Medium.ThermodynamicState outlet "Primary fluid thermodynamic state";
    /************* Tube pressure *****************/
    Modelica.SIunits.Pressure p "Pressure in the heat-exchanger: No pressure drop considered";
    /*********** GENERAL VARIABLES AND PARAMETER ***********/
    parameter Modelica.SIunits.Area A_ht "Heat transfer area";
    parameter Modelica.SIunits.Area A_cs_pf "Cross sectional area, working fluid side";
    parameter Modelica.SIunits.Area A_cs_sf "Cross sectional area, secondary fluid side";
    parameter Modelica.SIunits.Length L "Total length of the exchanger";
    final parameter Modelica.SIunits.Volume V = L * A_cs_pf "Volume";
    constant Real pi = Modelica.Constants.pi;
    /*************** Wall between the fluids *************************/
    parameter Modelica.SIunits.Mass M_wall "Total mass of metal walll between the fluids";
    parameter Modelica.SIunits.SpecificHeatCapacity c_wall "Heat capacity of the wall - assumed constant in the whole heat exchanger";
    parameter Modelica.SIunits.Density rho_wall "Density of the wall - assumed constant in the whole heat exchanger";
    parameter Modelica.SIunits.ThermalConductivity k_wall "Thermal conductivity of the metal wall";
    /********** MASS FLOWS *******************/
    //parameter Modelica.SIunits.MassFlowRate Mdotnom
    //    "Nominal fluid flow rate of the primary fluid";
    Modelica.SIunits.MassFlowRate mdot_pf_in "Mass flow at the inlet of the heat exchanger";
    Modelica.SIunits.MassFlowRate mdot_pf_out "Mass flow at the outlet of the heat exchanger";
    Modelica.SIunits.MassFlowRate mdot_sf "Mass flow of the secondary fluid of the heat exchanger";
    /*****************HEAT FLOW *************************/
    Modelica.SIunits.Power Q_pf "Thermal energy transfer with the wall";
    Modelica.SIunits.Power Q_sf "Total Thermal energy transfer with the wall from the secondary fluid";
    /***************** TEMPERATURES *************************/
    /** primary fluid **/
    Modelica.SIunits.Temperature T_pf_in "Temperatue at the inlet of the heat exchanger";
    Modelica.SIunits.Temperature T_pf "Mean Temperatue of the primary fluid";
    Modelica.SIunits.Temperature T_pf_out "Temperature at the outlet of the moving boundary model";
    /** Wall **/
    Modelica.SIunits.Temperature T_w "Mean Wall temperature";
    Modelica.SIunits.Temperature T_w_pf_in "Wall temperature at the inlet of the primary fluid";
    Modelica.SIunits.Temperature T_w_sf_in "Wall temperature at the inlet of the secondary fluid";
    Real dT_w(start = 10) "Wall temperature difference between inlet and outlet";
    /** Secondary fluid **/
    parameter Modelica.SIunits.Temperature dTsf_start = 10 "|Initialization|Temperature change";
    parameter Modelica.SIunits.Temperature Tsf_SU_start = 273.15 + 25 "|Initialization|Temperature inlet";
    Modelica.SIunits.Temperature T_sf_in "Tempearture at the inlet of the secondary fluid";
    Modelica.SIunits.Temperature T_sf "Mean Temperatue of the subcooled region secondary fluid";
    Modelica.SIunits.Temperature T_sf_out "Secondary fluid temperature at the exit of the heat exchanger";
    /*************** HEAT TRANSFER PARAMETER ****************/
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_pf "Heat transfer coefficient of the primary fluid";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_sf "Secondary fluid Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.QualityFactor eta_fin = 1 "Fin efficiency";
    /********** Densities ***************/
    Modelica.SIunits.Density rho_pf "Mean density of the subcooled region";
    Modelica.SIunits.Density rho_sf "Density of the secondary fluid at the inlet";
    /********** Specific heat capacity ***************/
    Modelica.SIunits.SpecificHeatCapacity cp_sf "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_pf "Specific heat capacity of the secondary fluid at the inlet";
    /************** Enthalpies ************************/
    Modelica.SIunits.SpecificEnthalpy h_pf "Mean Enthalpy of the primary fluid within the heat exchanger";
    Modelica.SIunits.SpecificEnthalpy h_pf_in "Enthalpy at the inlet of the heat exchanger of the primary fluid";
    Modelica.SIunits.SpecificEnthalpy h_pf_out "Enthalpy at the outlet of the heat exchanger of the primary fluid";
    /***** Enthalpy&Density****/
    /******************* Density derivatives ********************/
    Modelica.SIunits.DerDensityByEnthalpy drdh_pf "Density derivative with respect to enthalpy at constant pressure for the primary fluid";
    Modelica.SIunits.DerDensityByPressure drdp_pf "Density derivative with respect to pressure at constant enthalpy for the primary fluid";
    /************ Derivative with respect to time **************/
    Real drdt_pf "Density derivative with respect to time of the primary fluid";
    Real dhdt_pf "Enthalpy derivative with respect to time of the primary fluid";
    /********************** Epsilon-NTU method *********************/
    /******* Secondary fluid - Wall ****************/
    Real epsilon_sf "Epsilon Sub-cooled secondary fluid - wall Side";
    Real NTU_sf "NTU Sub-cooled secondary fluid-Wall side";
    /*********** Wall - Primary fluid ************/
    Real epsilon_pf "Epsilon Sub-cooled wall side - primary fluid";
    Real NTU_pf "NTU Sub-cooled wall side - primary fluid";
  initial equation
    h_pf = h_pf_in;
    h_pf = h_pf_out;
  equation
/********** Enthalpies ************/
    h_pf_in = inStream(InFlow.h_outflow);
    h_pf_out = OutFlow.h_outflow;
    h_pf_in = InFlow.h_outflow;
/*********** Temperatures *********/
    OutFlow_sf.T = T_sf_out;
    InFlow_sf.T = T_sf_in;
/**********Specific Heat Capacity *************/
    OutFlow_sf.cp = cp_sf;
    InFlow_sf.cp = cp_sf;
/*********** Density ***********/
    OutFlow_sf.rho = rho_sf;
    InFlow_sf.rho = rho_sf;
/*************** Pressure **************/
    InFlow.p = p;
    OutFlow.p = p;
/***********Mass Flow ******************/
    InFlow.m_flow = mdot_pf_in;
    OutFlow.m_flow = -mdot_pf_out;
    OutFlow_sf.Mdot = mdot_sf;
    InFlow_sf.Mdot = mdot_sf;
/*************** CONSTITUTIVE EQUATIONS ******************/
    h_pf = (h_pf_in + h_pf_out) / 2 "Mean enthalpy of the primary fluid";
    T_sf = (T_sf_in + T_sf_out) / 2 "Mean temperature of the secondary fluid";
    T_w = (T_w_pf_in + T_w_sf_in) / 2 "Mean wall temperature";
    dT_w = T_w_sf_in - T_w_pf_in;
    M_wall * c_wall * der(dT_w) = A_ht * U_pf * (T_pf_out - T_pf_in - dT_w) + A_ht * eta_fin * U_sf * (T_sf_in - T_sf_out - dT_w);
/**** SECONDARY FLUID - WALL HEAT TRANSFER ****/
    NTU_sf = U_sf * A_ht * eta_fin / (cp_sf * mdot_sf);
    epsilon_sf = 1 - exp(-NTU_sf);
    Q_sf = cp_sf * mdot_sf * epsilon_sf * (T_sf - T_w_pf_in);
    Q_sf = mdot_sf * cp_sf * (T_sf_in - T_sf_out);
/*** WALL - PRIMARY FLUID HEAT TRANSFER ***/
    NTU_pf = U_pf * A_ht / (cp_pf * (mdot_pf_in + mdot_pf_out) / 2);
    epsilon_pf = 1 - exp(-NTU_pf);
    Q_pf = cp_pf * (mdot_pf_in + mdot_pf_out) / 2 * epsilon_pf * (T_w_sf_in - T_pf_in) "Heat from the wall to the primary fluid";
/************ FLUID PROPERTIES ************************/
    PrimaryFluidState = Medium.setState_ph(p, h_pf);
    inlet = Medium.setState_ph(p, h_pf_in);
    outlet = Medium.setState_ph(p, h_pf_out);
/**Temperatures **/
    T_pf_in = Medium.temperature(inlet);
    T_pf = Medium.temperature(PrimaryFluidState);
    T_pf_out = Medium.temperature(outlet);
/** Densities **/
    rho_pf = Medium.density(PrimaryFluidState);
/** SpecificHeatCapacity **/
    cp_pf = Medium.specificHeatCapacityCp(PrimaryFluidState);
/**Derivatives**/
    drdh_pf = Medium.density_derh_p(PrimaryFluidState);
    drdp_pf = Medium.density_derp_h(PrimaryFluidState);
/************ MassBalance Primary fluid side ****************/
    V * drdt_pf = mdot_pf_in - mdot_pf_out;
    drdt_pf = drdp_pf * der(p) + drdh_pf * der(h_pf);
/********** EnergyBalance Primary fluid side ******************/
    V * (rho_pf * dhdt_pf + h_pf * drdt_pf - der(p)) = mdot_pf_in * h_pf_in - mdot_pf_out * h_pf_out + Q_pf;
    dhdt_pf = der(h_pf);
//Q_pf = (mdot_pf_in/2 + mdot_pf_out/2) * cp_pf * dhdt_pf ;//
/*******  WallEnergyBalance  **********/
    c_wall * M_wall * der(T_w) = Q_sf - Q_pf;
//c_wall * M_wall * der(T_w) = Q_sf - Q_pf + (T_w_sf_in - 2*T_w + T_w_pf_in) / (rho_wall * L^2) * M_wall * k_wall;//
/******************************* SUMMARY ***********************************/
  protected
    parameter Integer N = 3;
    Modelica.SIunits.Temperature[N] T_fluid_;
    Modelica.SIunits.Temperature[N] T_wall_;
    Modelica.SIunits.Temperature[N] T_sf_;
  public
    record SummaryBase
      replaceable Arrays T_profile;

      record Arrays
        parameter Integer n;
        Modelica.SIunits.Temperature[n] Twf;
        Modelica.SIunits.Temperature[n] Twall;
        Modelica.SIunits.Temperature[n] Tsf;
      end Arrays;

      //Modelica.SIunits.Pressure p_sf;
      Modelica.SIunits.Pressure p_wf;
      //Modelica.SIunits.Power Q_sf;
      //Modelica.SIunits.Power Q_wf;
    end SummaryBase;

    replaceable record SummaryClass = SummaryBase;
    SummaryClass Summary(T_profile(n = N, Twf = T_fluid_[:], Twall = T_wall_[:], Tsf = T_sf_[:]), p_wf = p);
  equation
/* Fluid temperature */
    T_fluid_[1] = T_pf_in;
    T_fluid_[2] = T_pf;
    T_fluid_[3] = T_pf_out;
/* Wall temperature */
    T_wall_[1] = T_w_pf_in;
    T_wall_[2] = T_w;
    T_wall_[3] = T_w_sf_in;
/* Secondary fluid temperature */
    T_sf_[1] = T_sf_out;
    T_sf_[2] = T_sf;
    T_sf_[3] = T_sf_in;
    annotation(
      Diagram(graphics),
      Icon(graphics = {Rectangle(extent = {{-100, -20}, {-32, -100}}, lineColor = {60, 121, 182}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-100, 20}, {100, -20}}, lineColor = {135, 135, 135}, fillPattern = FillPattern.CrossDiag), Rectangle(extent = {{-100, 100}, {100, 20}}, lineColor = {135, 135, 135}, fillColor = {255, 85, 85}, fillPattern = FillPattern.Solid), Line(points = {{-100, 56}, {-80, 36}, {-60, 56}, {-40, 36}, {-20, 56}, {0, 36}, {20, 56}, {40, 36}, {60, 56}, {80, 36}, {100, 56}}, color = {255, 0, 0}, smooth = Smooth.None, thickness = 0.5), Polygon(points = {{-18, 84}, {-18, 64}, {-32, 74}, {-18, 84}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-18, 74}, {34, 74}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{-18, -72}, {-10, -80}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{-100, -100}, {40, -100}, {-16, -20}, {-100, -20}, {-100, -24}, {-100, -100}}, lineColor = {0, 128, 255}, smooth = Smooth.None, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{32, -100}, {44, -100}, {102, -100}, {100, -20}, {-16, -20}, {32, -100}}, lineColor = {170, 255, 255}, smooth = Smooth.None, fillColor = {170, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{-97, -63}, {-75, -43}, {-57, -63}, {-37, -43}, {-17, -63}, {3, -43}, {23, -63}, {43, -43}, {63, -63}, {83, -43}, {103, -63}}, color = {0, 0, 255}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{4, -30}, {12, -38}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-2, -54}, {6, -62}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{42, -70}, {42, -90}, {56, -80}, {42, -70}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-8, -80}, {44, -80}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{28, -28}, {34, -34}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-14, -26}, {-2, -40}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{12, -64}, {22, -76}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{26, -78}, {38, -90}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{16, -42}, {28, -54}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{36, -62}, {42, -68}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid)}),
      Documentation(info = "<HTML>
            
           <p><big>Model <b>MBeva</b> is a moving boundaries heat exchanger model for evaporators.
           The model considers a fictitius heat transfer channel split up into three different section based on the phase state of the working fluid:
           <ul><li> Sub-cooled zone (SB)
           <li> Two-phase zone (TP)
           <li> SuperHeated zone (SH)
           </ul>
            <p><big><b>Pressure</b> and <b>enthalpy</b> are selected as state variables for each of the three zones. 
            
          <p><big>  The name moving boundary is derived from the fact that the interfaces between these sections do not
  have a fixed spatial position but merely a fixed thermodynamic location depending
  on the presence of liquid and gaseous fluid, respectively. The actual existence of
  a certain section and its length are determined based on the fluid state resulting
  in variable sectioning. A fixed total length superimposes the required boundary
  condition to calculate the length of each section.
   
   <p><big>The assumptions for this model are:
           <ul><li> The tube is cylindrical with a constant cross sectional area
           <li> The velocity of the fluid is uniform on the cross sectional area
           <li> The enthalpy of the fluid is linear in each region of the tube (sub-cooled, two-phase, super-heated)
           <li> Pressure is considered constant
           <li> Constant void fraction
           <li> Thermal energy accumulation in the metal wall is taken into account
           <li> The secondary fluid is treated as a constant heat capacity fluid
           <li> The heat transfer between the secondary fluid and the metal wall and between the metal wall and the primary fluid is computed using the epsilon-NTU method
           </ul>
          </HTML>"));
  end HEX_1f;

  model HEX_2f "Counter current simple HEX model. No fluid undergo any phase change. The primary fluid is a proper Modelica fluid, and should be in liquid state (hence the ''secondary fluid - to - liquid'' denomination)"
    /**************** MEDIUM ***************************/
    replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialMedium annotation(
       choicesAllMatching = true);
    /***************** PORTS ******************/
    Modelica.Fluid.Interfaces.FluidPort_a InFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-108, -74}, {-88, -54}})));
    Modelica.Fluid.Interfaces.FluidPort_b OutFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{90, -74}, {110, -54}})));
    ShipEnergySystems.Other.FlangeCdot_B OutFlow_sf annotation(
      Placement(transformation(extent = {{-110, 50}, {-90, 70}})));
    ShipEnergySystems.Other.FlangeCdot_A InFlow_sf annotation(
      Placement(transformation(extent = {{90, 50}, {110, 70}})));
    /***********THERMODYNAMIC STATES *********************/
    Medium.ThermodynamicState state_sl "Primary fluid thermodynamic state, subcooled liquid";
    Medium.SaturationProperties sat "Primary fluid thermodynamic state, evaporation";
    Medium.ThermodynamicState state_ev "Primary fluid thermodynamic state, evaporation";
    Medium.ThermodynamicState inlet "Primary fluid thermodynamic state";
    Medium.ThermodynamicState outlet "Primary fluid thermodynamic state";
    /************* Tube pressure *****************/
    Modelica.SIunits.Pressure p "Pressure in the heat-exchanger: No pressure drop considered";
    /*********** GENERAL VARIABLES AND PARAMETER ***********/
    parameter Modelica.SIunits.Area A "Heat transfer area";
    Modelica.SIunits.Area A_sl "Heat transfer area";
    Modelica.SIunits.Area A_ev "Heat transfer area";
    parameter Modelica.SIunits.Length L "Total length of the exchanger";
    Modelica.SIunits.Length L_sl "Length of the HEX, subcooled liquid part";
    Modelica.SIunits.Length L_ev "Length of the HEX, evaporation part";
    final parameter Modelica.SIunits.Length D = A / pi / L "Diameter";
    constant Real pi = Modelica.Constants.pi;
    /*************** Wall between the fluids *************************/
    parameter Modelica.SIunits.Mass M_wall "Total mass of metal walll between the fluids";
    parameter Modelica.SIunits.SpecificHeatCapacity c_wall "Heat capacity of the wall - assumed constant in the whole heat exchanger";
    parameter Modelica.SIunits.Density rho_wall "Density of the wall - assumed constant in the whole heat exchanger";
    parameter Modelica.SIunits.ThermalConductivity k_wall "Thermal conductivity of the metal wall";
    /********** MASS FLOWS *******************/
    Modelica.SIunits.MassFlowRate mdot_pf_in "Mass flow at the inlet of the heat exchanger";
    Modelica.SIunits.MassFlowRate mdot_pf_out "Mass flow at the outlet of the heat exchanger";
    Modelica.SIunits.MassFlowRate mdot_pf_sl2ev "Mass flow of the primary fluid at the interface between the two regions";
    Modelica.SIunits.MassFlowRate mdot_sf "Mass flow of the secondary fluid of the heat exchanger";
    /*****************HEAT FLOW *************************/
    Modelica.SIunits.Power Q_pf_sl "Thermal energy transfer with the wall, subcooled liquid part";
    Modelica.SIunits.Power Q_pf_ev "Thermal energy transfer with the wall, two-phase part";
    Modelica.SIunits.Power Q_sf_sl "Energy transfer with the wall from the secondary fluid, subcooled liquid part";
    Modelica.SIunits.Power Q_sf_ev "Energy transfer with the wall from the secondary fluid, two-phase part";
    /***************** TEMPERATURES *************************/
    /** primary fluid **/
    Modelica.SIunits.Temperature T_pf_in "Temperatue at the inlet of the heat exchanger";
    Modelica.SIunits.Temperature T_pf_sl "Mean Temperatue of the primary fluid in the subcooled liquid region";
    Modelica.SIunits.Temperature T_pf_sat "Mean Temperatue of the primary fluid in the evaporation region";
    Modelica.SIunits.Temperature T_pf_out "Temperature at the outlet of the moving boundary model";
    /** Wall **/
    Modelica.SIunits.Temperature T_w_sl "Mean Wall temperature, subcooled liquid region";
    Modelica.SIunits.Temperature T_w_ev "Mean Wall temperature, evaporation region";
    Modelica.SIunits.Temperature T_w_sl2ev "Wall temperature at the interface";
    Modelica.SIunits.Temperature T_w_pf_in "Wall temperature at the inlet of the primary fluid";
    Modelica.SIunits.Temperature T_w_sf_in "Wall temperature at the inlet of the secondary fluid";
    Real dT_w_sl(start = 10) "Wall temperature difference along the subcooled liquid region";
    Real dT_w_ev(start = 10) "Wall temperature difference along the evaporation region";
    /** Secondary fluid **/
    Modelica.SIunits.Temperature T_sf_in "Tempearture at the inlet of the secondary fluid";
    Modelica.SIunits.Temperature T_sf_out "Secondary fluid temperature at the exit of the heat exchanger";
    Modelica.SIunits.Temperature T_sf_sl "Mean Temperatue of the secondary fluid, subcooled liquid region";
    Modelica.SIunits.Temperature T_sf_ev "Mean Temperatue of the secondary fluid, evaporation region";
    Modelica.SIunits.Temperature T_sf_sl2ev "Temperatue of the secondary fluid, interface";
    /*************** HEAT TRANSFER PARAMETER ****************/
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_pf_sl "Heat transfer coefficient of the primary fluid, subcooled region";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_pf_ev "Heat transfer coefficient of the primary fluid, evaporation";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_sf "Secondary fluid Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.QualityFactor eta_fin = 1 "Fin efficiency";
    /********** Densities ***************/
    Modelica.SIunits.Density rho_pf_sl "Mean density of the primary fluid, subcooled liquid region";
    Modelica.SIunits.Density rho_pf_ls "Density of the primary fluid, saturated liquid";
    Modelica.SIunits.Density rho_pf_vs "Density of the primary fluid, saturated vapor";
    Modelica.SIunits.Density rho_sf "Density of the secondary fluid at the inlet";
    /********** Specific heat capacity ***************/
    Modelica.SIunits.SpecificHeatCapacity cp_sf "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_pf_sl "Specific heat capacity of the primary fluid in the subcooled liquid region";
    /************** Enthalpies ************************/
    Modelica.SIunits.SpecificEnthalpy h_pf_in "Enthalpy at the inlet of the heat exchanger of the primary fluid";
    Modelica.SIunits.SpecificEnthalpy h_pf_out "Enthalpy at the outlet of the heat exchanger of the primary fluid";
    Modelica.SIunits.SpecificEnthalpy h_pf_sl "Mean Enthalpy of the primary fluid in the subcooled liquid region";
    Modelica.SIunits.SpecificEnthalpy h_pf_ev "Mean Enthalpy of the primary fluid in the evaporation region";
    Modelica.SIunits.SpecificEnthalpy h_pf_ls "Enthalpy of the primary fluid in saturated liquid conditions";
    Modelica.SIunits.SpecificEnthalpy h_pf_vs "Enthalpy of the primary fluid in saturated liquid conditions";
    /******************* Density derivatives ********************/
    Modelica.SIunits.DerDensityByEnthalpy drdh_pf_sl "Density derivative with respect to enthalpy at constant pressure for the primary fluid in the subcooled liquid region";
    Modelica.SIunits.DerDensityByPressure drdp_pf_sl "Density derivative with respect to pressure at constant enthalpy for the primary fluid in the subcooled liquid region";
    Modelica.SIunits.DerDensityByPressure drdp_pf_ls "Density derivative with respect to pressure at constant enthalpy for the primary fluid in along the saturated liquid curve";
    Modelica.SIunits.DerDensityByPressure drdp_pf_vs "Density derivative with respect to pressure at constant enthalpy for the primary fluid along the saturated vapor curve";
    Modelica.SIunits.DerEnthalpyByPressure dhdp_pf_ls "Enthalpy derivative with respect to pressure for the primary fluid in along the saturated liquid curve";
    Modelica.SIunits.DerEnthalpyByPressure dhdp_pf_vs "Enthalpy derivative with respect to pressure for the primary fluid along the saturated vapor curve";
    /***********  Void fraction calculation **************/
    Real void "Void fraction in the two-phase exchange region";
    Real mu "Density ratio of saturated liquid over saturated vapor";
    /********************** Epsilon-NTU method *********************/
    /******* Secondary fluid - Wall ****************/
    Real epsilon_sf_sl "Epsilon secondary fluid - wall Side, subcooled liquid region";
    Real NTU_sf_sl "NTU secondary fluid-Wall side, subcooled liquid region";
    Real epsilon_sf_ev "Epsilon secondary fluid - wall Side, evaporation region";
    Real NTU_sf_ev "NTU secondary fluid-Wall side, evaporation region";
    /*********** Wall - Primary fluid ************/
    Real epsilon_pf_sl "Epsilon wall side - primary fluid, , subcooled liquid region";
    Real NTU_pf_sl "NTU wall side - primary fluid, subcooled liquid region";
  equation
/********** Surface balance ************/
    L_sl + L_ev - L = 0;
    A_sl + A_ev - A = 0;
    A_sl = A * L_sl / L;
/********** Enthalpies ************/
    h_pf_in = inStream(InFlow.h_outflow);
    h_pf_out = OutFlow.h_outflow;
    h_pf_in = InFlow.h_outflow;
/*********** Temperatures *********/
    OutFlow_sf.T = T_sf_out;
    InFlow_sf.T = T_sf_in;
/**********Specific Heat Capacity *************/
    OutFlow_sf.cp = cp_sf;
    InFlow_sf.cp = cp_sf;
/*********** Density ***********/
    OutFlow_sf.rho = rho_sf;
    InFlow_sf.rho = rho_sf;
/*************** Pressure **************/
    InFlow.p = p;
    OutFlow.p = p;
/***********Mass Flow ******************/
    InFlow.m_flow = mdot_pf_in;
    OutFlow.m_flow = -mdot_pf_out;
    OutFlow_sf.Mdot = mdot_sf;
    InFlow_sf.Mdot = mdot_sf;
/*************** CONSTITUTIVE EQUATIONS ******************/
    h_pf_sl = (h_pf_in + h_pf_ls) / 2 "Mean enthalpy of the primary fluid, subcooled liquid region";
    h_pf_ev = (h_pf_ls + h_pf_out) / 2 "Mean enthalpy of the primary fluid, evaporation region";
    T_sf_sl = (T_sf_in + T_sf_sl2ev) / 2 "Mean temperature of the secondary fluid";
    T_sf_ev = (T_sf_sl2ev + T_sf_out) / 2 "Mean temperature of the secondary fluid";
    T_w_sl = (T_w_pf_in + T_w_sl2ev) / 2 "Mean wall temperature";
    T_w_ev = (T_w_sl2ev + T_w_sf_in) / 2 "Mean wall temperature";
/************* WALL NON-UNIFORM TEMPERATURE FORMULATION **********/
    dT_w_sl = T_w_sl2ev - T_w_pf_in;
    dT_w_ev = T_w_sf_in - T_w_sl2ev;
    M_wall * c_wall * L_sl / L * der(dT_w_sl) = A_sl * U_pf_sl * (T_pf_sat - T_pf_in - dT_w_sl) + A_sl * eta_fin * U_sf * (T_sf_sl2ev - T_sf_out - dT_w_sl);
    M_wall * c_wall * L_ev / L * der(dT_w_ev) = A_ev * U_pf_ev * (T_pf_out - T_pf_sat - dT_w_ev) + A_ev * eta_fin * U_sf * (T_sf_in - T_sf_sl2ev - dT_w_ev);
/**** SECONDARY FLUID - WALL HEAT TRANSFER ****/
// Subcooled liquid region //
    NTU_sf_sl = U_sf * A_sl * eta_fin / (cp_sf * mdot_sf);
    epsilon_sf_sl = 1 - exp(-NTU_sf_sl);
    Q_sf_sl = cp_sf * mdot_sf * epsilon_sf_sl * (T_sf_sl2ev - T_w_pf_in);
    Q_sf_sl = mdot_sf * cp_sf * (T_sf_sl2ev - T_sf_out);
// Evaporation region //
    NTU_sf_ev = U_sf * A_ev * eta_fin / (cp_sf * mdot_sf);
    epsilon_sf_ev = 1 - exp(-NTU_sf_ev);
    Q_sf_ev = cp_sf * mdot_sf * epsilon_sf_ev * (T_sf_in - T_w_sl2ev);
    Q_sf_ev = mdot_sf * cp_sf * (T_sf_in - T_sf_sl2ev);
/*** WALL - PRIMARY FLUID HEAT TRANSFER ***/
// Subcooled liquid region //
    NTU_pf_sl = U_pf_sl * A_sl / (cp_pf_sl * mdot_pf_in / 2);
    epsilon_pf_sl = 1 - exp(-NTU_pf_sl);
    Q_pf_sl = cp_pf_sl * (mdot_pf_in + mdot_pf_sl2ev) / 2 * epsilon_pf_sl * (T_w_sl2ev - T_pf_in);
// Evaporation region //
    Q_pf_ev = U_pf_ev * A_ev * (T_w_ev - T_pf_sat);
/************ FLUID PROPERTIES ************************/
    state_sl = Medium.setState_ph(p, h_pf_sl);
    state_ev = Medium.setState_ph(p, h_pf_ev);
    sat = Medium.setSat_p(p);
    inlet = Medium.setState_ph(p, h_pf_in);
    outlet = Medium.setState_ph(p, h_pf_out);
/**Temperatures **/
    T_pf_in = Medium.temperature(inlet);
    T_pf_sl = Medium.temperature(state_sl);
    T_pf_sat = Medium.saturationTemperature_sat(sat);
    T_pf_out = Medium.temperature(outlet);
/** Densities **/
    rho_pf_sl = Medium.density(state_sl);
    rho_pf_ls = Medium.bubbleDensity(sat);
    rho_pf_vs = Medium.dewDensity(sat);
    mu = rho_pf_ls / rho_pf_vs;
/** Enthalpies **/
    h_pf_ls = Medium.bubbleEnthalpy(sat);
    h_pf_vs = Medium.dewEnthalpy(sat);
/** SpecificHeatCapacity **/
    cp_pf_sl = Medium.specificHeatCapacityCp(state_sl);
/**Derivatives**/
    drdh_pf_sl = Medium.density_derh_p(state_sl);
    drdp_pf_sl = Medium.density_derp_h(state_sl);
    drdp_pf_ls = Medium.dBubbleDensity_dPressure(sat);
    dhdp_pf_ls = Medium.dBubbleEnthalpy_dPressure(sat);
    dhdp_pf_vs = Medium.dDewEnthalpy_dPressure(sat);
    drdp_pf_vs = Medium.dDewDensity_dPressure(sat);
/************ MassBalance Primary fluid side ****************/
// Subcooled liquid region //
    pi * D ^ 2 / 4 * (L_sl * (drdp_pf_sl + drdh_pf_sl * dhdp_pf_ls) * der(p) + (rho_pf_sl - rho_pf_ls) * der(L_sl) + 0.5 * L_sl * drdh_pf_sl * der(h_pf_in)) = mdot_pf_in - mdot_pf_sl2ev;
// Updated
// Evaporation region //
    pi * D ^ 2 / 4 * ((rho_pf_ls - rho_pf_vs) * der(L_sl) + (1 - void) * (rho_pf_ls - rho_pf_vs) * der(L_ev) + L_ev * (void * drdp_pf_vs + (1 - void) * drdp_pf_ls) * der(p)) = mdot_pf_sl2ev - mdot_pf_out;
// Updated
/********** EnergyBalance Primary fluid side ******************/
// Subcooled liquid region //
    0.5 * pi * D ^ 2 / 4 * ((rho_pf_sl * (h_pf_in + h_pf_ls) - 2 * rho_pf_sl * h_pf_ls) * der(L_sl) + (rho_pf_sl * L_sl + drdh_pf_sl) * der(h_pf_in) + L_sl * (rho_pf_sl * dhdp_pf_ls + (h_pf_in + h_pf_ls) * (drdp_pf_sl + 0.5 * drdh_pf_sl * dhdp_pf_ls - 2))) * der(p) = mdot_pf_in * h_pf_in - mdot_pf_sl2ev * h_pf_ls + Q_pf_sl;
// Updated
// Evaporation region //
    0.5 * pi * D ^ 2 / 4 * (L_ev * (void * (h_pf_vs * drdp_pf_vs + rho_pf_vs * dhdp_pf_vs) + (1 - void) * (h_pf_ls * drdp_pf_ls + rho_pf_ls * dhdp_pf_ls) - 1) * der(p) + (void * rho_pf_vs * h_pf_vs + (1 - void) * rho_pf_ls * h_pf_ls) * der(L_sl) + (1 - void) * (rho_pf_ls * h_pf_ls - rho_pf_vs * h_pf_vs) * der(L_ev)) = mdot_pf_sl2ev * h_pf_ls - mdot_pf_out * h_pf_out + Q_pf_sl;
/*******  WallEnergyBalance  **********/
    c_wall * M_wall * (L_sl * der(T_w_sl) + T_w_sl * der(L_sl)) / L = Q_sf_sl - Q_pf_sl;
    c_wall * M_wall * (L_ev * der(T_w_ev) + T_w_ev * der(L_ev)) / L = Q_sf_ev - Q_pf_ev;
/****** Calculation of the void fraction ******/
    1 - void = (1 + mu ^ (-2 / 3) * (2 / 3 * log(1 / mu) - 1)) / (mu ^ (-2 / 3) - 1) ^ 2;
/******************************* SUMMARY ***********************************/
  protected
    parameter Integer N = 5;
    Modelica.SIunits.Temperature[N] T_fluid_;
    Modelica.SIunits.Temperature[N] T_wall_;
    Modelica.SIunits.Temperature[N] T_sf_;
  public
    record SummaryBase
      replaceable Arrays T_profile;

      record Arrays
        parameter Integer n;
        Modelica.SIunits.Temperature[n] Twf;
        Modelica.SIunits.Temperature[n] Twall;
        Modelica.SIunits.Temperature[n] Tsf;
      end Arrays;

      //Modelica.SIunits.Pressure p_sf;
      Modelica.SIunits.Pressure p_wf;
      //Modelica.SIunits.Power Q_sf;
      //Modelica.SIunits.Power Q_wf;
    end SummaryBase;

    replaceable record SummaryClass = SummaryBase;
    SummaryClass Summary(T_profile(n = N, Twf = T_fluid_[:], Twall = T_wall_[:], Tsf = T_sf_[:]), p_wf = p);
  equation
/* Fluid temperature */
    T_fluid_[1] = T_pf_in;
    T_fluid_[2] = T_pf_sl;
    T_fluid_[3] = T_pf_sat;
    T_fluid_[4] = T_pf_sat;
    T_fluid_[5] = T_pf_out;
/* Wall temperature */
    T_wall_[1] = T_w_pf_in;
    T_wall_[2] = T_w_sl;
    T_wall_[3] = T_w_sl2ev;
    T_wall_[4] = T_w_ev;
    T_wall_[5] = T_w_sf_in;
/* Secondary fluid temperature */
    T_sf_[1] = T_sf_out;
    T_sf_[2] = T_sf_sl;
    T_sf_[3] = T_sf_sl2ev;
    T_sf_[4] = T_sf_ev;
    T_sf_[5] = T_sf_in;
    annotation(
      Diagram(graphics),
      Icon(graphics = {Rectangle(extent = {{-100, -20}, {-32, -100}}, lineColor = {60, 121, 182}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-100, 20}, {100, -20}}, lineColor = {135, 135, 135}, fillPattern = FillPattern.CrossDiag), Rectangle(extent = {{-100, 100}, {100, 20}}, lineColor = {135, 135, 135}, fillColor = {255, 85, 85}, fillPattern = FillPattern.Solid), Line(points = {{-100, 56}, {-80, 36}, {-60, 56}, {-40, 36}, {-20, 56}, {0, 36}, {20, 56}, {40, 36}, {60, 56}, {80, 36}, {100, 56}}, color = {255, 0, 0}, smooth = Smooth.None, thickness = 0.5), Polygon(points = {{-18, 84}, {-18, 64}, {-32, 74}, {-18, 84}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-18, 74}, {34, 74}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{-18, -72}, {-10, -80}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{-100, -100}, {40, -100}, {-16, -20}, {-100, -20}, {-100, -24}, {-100, -100}}, lineColor = {0, 128, 255}, smooth = Smooth.None, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{32, -100}, {44, -100}, {102, -100}, {100, -20}, {-16, -20}, {32, -100}}, lineColor = {170, 255, 255}, smooth = Smooth.None, fillColor = {170, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{-97, -63}, {-75, -43}, {-57, -63}, {-37, -43}, {-17, -63}, {3, -43}, {23, -63}, {43, -43}, {63, -63}, {83, -43}, {103, -63}}, color = {0, 0, 255}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{4, -30}, {12, -38}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-2, -54}, {6, -62}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{42, -70}, {42, -90}, {56, -80}, {42, -70}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-8, -80}, {44, -80}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{28, -28}, {34, -34}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-14, -26}, {-2, -40}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{12, -64}, {22, -76}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{26, -78}, {38, -90}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{16, -42}, {28, -54}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{36, -62}, {42, -68}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid)}),
      Documentation(info = "<HTML>
          
         <p><big>Model <b>MBeva</b> is a moving boundaries heat exchanger model for evaporators.
         The model considers a fictitius heat transfer channel split up into three different section based on the phase state of the working fluid:
         <ul><li> Sub-cooled zone (SB)
         <li> Two-phase zone (TP)
         <li> SuperHeated zone (SH)
         </ul>
          <p><big><b>Pressure</b> and <b>enthalpy</b> are selected as state variables for each of the three zones. 
          
        <p><big>  The name moving boundary is derived from the fact that the interfaces between these sections do not
have a fixed spatial position but merely a fixed thermodynamic location depending
on the presence of liquid and gaseous fluid, respectively. The actual existence of
a certain section and its length are determined based on the fluid state resulting
in variable sectioning. A fixed total length superimposes the required boundary
condition to calculate the length of each section.
 
 <p><big>The assumptions for this model are:
         <ul><li> The tube is cylindrical with a constant cross sectional area
         <li> The velocity of the fluid is uniform on the cross sectional area
         <li> The enthalpy of the fluid is linear in each region of the tube (sub-cooled, two-phase, super-heated)
         <li> Pressure is considered constant
         <li> Constant void fraction
         <li> Thermal energy accumulation in the metal wall is taken into account
         <li> The secondary fluid is treated as a constant heat capacity fluid
         <li> The heat transfer between the secondary fluid and the metal wall and between the metal wall and the primary fluid is computed using the epsilon-NTU method
         </ul>
        </HTML>"));
  end HEX_2f;

  model HEX_3f "Counter current Moving Boundary model: Fluid enters subcooled and exits in super-heated conditons. The model consider the fluid in one side and the metal wall. The secondary fluid is a constant specific heat fluid"
    /**************** MEDIUM ***************************/
    replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialMedium annotation(
       choicesAllMatching = true);
    /***************** PORTS ******************/
    Modelica.Fluid.Interfaces.FluidPort_a InFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-108, -74}, {-88, -54}})));
    Modelica.Fluid.Interfaces.FluidPort_b OutFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{90, -74}, {110, -54}})));
    ShipEnergySystemsModelicaLibrary.Other.FlangeCdot_B OutFlow_sf annotation(
      Placement(transformation(extent = {{-110, 50}, {-90, 70}})));
    ShipEnergySystemsModelicaLibrary.Other.FlangeCdot_A InFlow_sf annotation(
      Placement(transformation(extent = {{90, 50}, {110, 70}})));
    /***********THERMODYNAMIC STATES *********************/
    Medium.ThermodynamicState SubCooled "Subcooled thermodynamic state";
    Medium.SaturationProperties sat;
    Medium.ThermodynamicState SuperHeated "Superheated thermodynamic state";
    Medium.ThermodynamicState outlet "Outlet thermodynamic state";
    /************* Tube pressure *****************/
    Modelica.SIunits.Pressure p(start = 0) "Pressure in the heat-exchanger: No pressure drop considered";
    /*********** GENERAL VARIABLES AND PARAMETER ***********/
    parameter Modelica.SIunits.Area A "Cross-sectional area for both side of the heat exchanger";
    parameter Modelica.SIunits.Length L "Total length of the exchanger";
    final parameter Modelica.SIunits.Length D = 2 * sqrt(A / pi) "Diameter";
    constant Real pi = Modelica.Constants.pi;
    /*************** Wall between the fluids *************************/
    parameter Modelica.SIunits.Mass M_tot "Total mass of metal walll between the fluids";
    parameter Modelica.SIunits.SpecificHeatCapacity c_wall "Heat capacity of the wall - assumed constant in the whole heat exchanger";
    parameter Modelica.SIunits.Density rho_wall "Density of the wall - assumed constant in the whole heat exchanger";
    /********** MASS FLOWS *******************/
    //parameter Modelica.SIunits.MassFlowRate Mdotnom
    //    "Nominal fluid flow rate of the primary fluid";
    Modelica.SIunits.MassFlowRate M_dot_su "Mass flow at the inlet of the heat exchanger";
    Modelica.SIunits.MassFlowRate M_dot_A "Mass flow at the interface between sub-cooled and two-phase";
    Modelica.SIunits.MassFlowRate M_dot_B "Mass flow at the interface between sub-cooled and two-phase";
    Modelica.SIunits.MassFlowRate M_dot_ex "Mass flow at the outlet of the heat exchanger";
    Modelica.SIunits.MassFlowRate Mdot_sf "Mass flow of the secondary fluid of the heat exchanger";
    /*****************HEAT FLOW *************************/
    Modelica.SIunits.Power Q_SB "Thermal energy transfer with the wall";
    Modelica.SIunits.Power Q_TP "Thermal energy transfer with the wall in the two-phase region";
    Modelica.SIunits.Power Q_SH "Thermal energy transfer with the wall in the superheated region";
    Modelica.SIunits.Power Q_tot "Total Thermal energy transfer with the wall from the primary fluid";
    Modelica.SIunits.Power Q_tot_sf "Total Thermal energy transfer with the wall from the secondary fluid";
    Modelica.SIunits.Power Qsf_SH;
    Modelica.SIunits.Power Qsf_TP;
    Modelica.SIunits.Power Qsf_SB;
    /***************** TEMPERATURES *************************/
    /** primary fluid **/
    Modelica.SIunits.Temperature T_SU "Temperatue at the inlet of the heat exchanger";
    Modelica.SIunits.Temperature T_SB "Mean Temperatue of the subcooled region";
    Modelica.SIunits.Temperature T_TP "Mean Temperatue of the two-phase region";
    Modelica.SIunits.Temperature T_SH "Mean Temperatue of the superheated region";
    Modelica.SIunits.Temperature TwA "Wall Temperature at the interface between sub-cooled and two-phase";
    Modelica.SIunits.Temperature TwB "Wall Temperature at the interface between two-phase and super-heated";
    Modelica.SIunits.Temperature T_out "Temperature at the outlet of the moving boundary model";
    /** Wall **/
    Modelica.SIunits.Temperature TwSB(start = 0) "Mean Wall temperature in the sub-cooled region";
    Modelica.SIunits.Temperature TwTP(start = 0) "Mean Wall Temperature in the two-phase region";
    Modelica.SIunits.Temperature TwSH(start = 0) "Mean Wall Temperature in the super-heated region";
    /** Secondary fluid **/
    parameter Modelica.SIunits.Temperature dTsf_start = 10 "|Initialization|Temperature change";
    parameter Modelica.SIunits.Temperature Tsf_SU_start = 273.15 + 25 "|Initialization|Temperature inlet";
    Modelica.SIunits.Temperature Tsf_SU "Tempearture at the inlet of the secondary fluid";
    Modelica.SIunits.Temperature Tsf_SH "Mean Temperatue of the superheated region secondary fluid";
    Modelica.SIunits.Temperature Tsf_B "Secondary fluid temperature at the interface between two-phase and super-heated";
    Modelica.SIunits.Temperature Tsf_TP "Mean Temperatue of the two-phase region secondary fluid";
    Modelica.SIunits.Temperature Tsf_A "Secondary fluid temperature at the interface between two-phase and super-heated";
    Modelica.SIunits.Temperature Tsf_SB "Mean Temperatue of the subcooled region secondary fluid";
    Modelica.SIunits.Temperature Tsf_EX(start = Tsf_SU_start - dTsf_start) "Secondary fluid temperature at the exit of the heat exchanger";
    /*************** HEAT TRANSFER PARAMETER ****************/
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_SB "Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_TP "Heat transfer coefficient two-phase side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_SH "Heat transfer coefficient Superheated side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer Usf "Secondary fluid Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.QualityFactor ETA = 1 "Fin efficiency";
    /*************Tube region length ********************/
    Modelica.SIunits.Length L_SB(start = 0) "Length of the subcooled region";
    Modelica.SIunits.Length L_TP(start = 0) "Length of the two-phase region";
    Modelica.SIunits.Length L_SH "Length of the Superheated region";
    /********** Densities ***************/
    Modelica.SIunits.Density rho_SB "Mean density of the subcooled region";
    Modelica.SIunits.Density rho_LS "Bubble density";
    Modelica.SIunits.Density rho_VS "Dew density";
    Modelica.SIunits.Density rho_TP "Mean density of the two-phase region";
    Modelica.SIunits.Density rho_SH "Mean density of the superheated region";
    Modelica.SIunits.Density rho_sf "Density of the secondary fluid at the inlet";
    /********** Specific heat capacity ***************/
    Modelica.SIunits.SpecificHeatCapacity cp_sf "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_SB "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_TP "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_SH "Specific heat capacity of the secondary fluid at the inlet";
    /************** Enthalpies ************************/
    Modelica.SIunits.SpecificEnthalpy h_SB "Mean Enthalpy of the subcooled region";
    Modelica.SIunits.SpecificEnthalpy h_LS "Bubble enthalpy";
    Modelica.SIunits.SpecificEnthalpy h_VS "Dew enthalpy";
    //Modelica.SIunits.SpecificEnthalpy h_TP "Mean Enthalpy of the two-phase region";
    Modelica.SIunits.SpecificEnthalpy h_SH "Mean Enthalpy of the Superheated region";
    Modelica.SIunits.SpecificEnthalpy h_SU "Enthalpy at the inlet of the heat exchanger of the primary fluid";
    Modelica.SIunits.SpecificEnthalpy h_EX(start = 0) "Enthalpy at the outlet of the heat exchanger of the primary fluid";
    /***** Enthalpy&Density****/
    Real rhoh_TP "Density times enthalpy in the two-phase region";
    /******************* Density derivatives ********************/
    Modelica.SIunits.DerDensityByEnthalpy drdh_SB "Density derivative with respect to enthalpy at constant pressure in the subcooled region";
    Modelica.SIunits.DerDensityByPressure drdp_SB "Density derivative with respect to pressure at constant enthalpy in the subcooled region";
    Modelica.SIunits.DerDensityByPressure drdp_LS "Bubble point density derivative with respect to pressure";
    Modelica.SIunits.DerDensityByPressure drdp_VS "Dew point density derivative with respect to pressure";
    Modelica.SIunits.DerDensityByEnthalpy drdh_SH "Density derivative with respect to enthalpy at constant pressure in the superheated region";
    Modelica.SIunits.DerDensityByPressure drdp_SH "Density derivative with respect to pressure at constant enthalpy in the superheated region";
    /*************** Enthalpy derivatives *******************/
    Modelica.SIunits.DerEnthalpyByPressure dhdp_LS "Bubble point specific enthalpy derivative with respect to pressure";
    Modelica.SIunits.DerEnthalpyByPressure dhdp_VS "Dew point specific enthalpy derivative with respect to pressure";
    /************ Derivative with respect to time **************/
    Real drdt_SB "Density derivative with respect to time in the subcooled region";
    Real dhdt_SB "Enthalpy derivative with respect to time in the subcooled region";
    Real drdt_TP "Density derivative with respect to time in the two-phase region";
    Real drhdt_TP "Density and Enthalpy derivative with respect to time in the two-phase region";
    Real drdt_SH "Density derivative with respect to time in the superheated region";
    Real dhdt_SH "Enthalpy derivative with respect to time in the superheated region";
    /********************** VOID fraction *****************************/
    parameter Real Void "Void fraction assumed constant for now";
    parameter Real dVoid_dp = 0 "Void fraction derivative with respect to pressure";
    parameter Real dVoid_dh = 0 "Void fraction derivative with respect to enthalpy";
    /********************** Epsilon-NTU method *********************/
    /******* Secondary fluid - Wall ****************/
    Real epsilon_SB_sf "Epsilon Sub-cooled secondary fluid - wall Side";
    Real NTU_SB_sf "NTU Sub-cooled secondary fluid-Wall side";
    Real epsilon_TP_sf "Epsilon Two-Phase secondary fluid - wall Side";
    Real NTU_TP_sf "NTU Two-Phase secondary fluid-Wall side";
    Real epsilon_SH_sf "Epsilon Super-heated secondary fluid - wall Side";
    Real NTU_SH_sf "NTU Super-heated secondary fluid-Wall side";
    /*********** Wall - Primary fluid ************/
    Real epsilon_SB_pf "Epsilon Sub-cooled wall side - primary fluid";
    Real NTU_SB_pf "NTU Sub-cooled wall side - primary fluid";
    Real epsilon_TP_pf "Epsilon Two-Phase wall side - primary fluid";
    Real NTU_TP_pf "NTU Two-Phase wall side - primary fluid";
    //
    //
    Real epsilon_SH_pf "Epsilon Super-heated wall side - primary fluid";
    Real NTU_SH_pf "NTU Super-heated secondary wall side - primary fluid";
    //initial equation
    // M_dot_A = M_dot_su;
    // M_dot_B = M_dot_su;
    // M_dot_ex = M_dot_su;
  equation
/****************** BOUNDARY EQUATION *************************/
    L = L_SB + L_TP + L_SH;
/********** Enthalpies ************/
    h_SU = inStream(InFlow.h_outflow);
    h_EX = OutFlow.h_outflow;
    h_SU = InFlow.h_outflow;
/*********** Temperatures *********/
    OutFlow_sf.T = Tsf_EX;
    InFlow_sf.T = Tsf_SU;
/**********Specific Heat Capacity *************/
    OutFlow_sf.cp = cp_sf;
    InFlow_sf.cp = cp_sf;
/*********** Density ***********/
    OutFlow_sf.rho = rho_sf;
    InFlow_sf.rho = rho_sf;
/*************** Pressure **************/
    InFlow.p = p;
    OutFlow.p = p;
/***********Mass Flow ******************/
    InFlow.m_flow = M_dot_su;
    OutFlow.m_flow = -M_dot_ex;
    OutFlow_sf.Mdot = Mdot_sf;
    InFlow_sf.Mdot = Mdot_sf;
/*************** CONSTITUTIVE EQUATIONS ******************/
    h_SB = (h_SU + h_LS) / 2 "Mean enthalpy of the subcooled region";
    rho_TP = rho_VS * Void + rho_LS * (1 - Void);
    rhoh_TP = rho_VS * h_VS * Void + rho_LS * h_LS * (1 - Void);
    h_SH = (h_EX + h_VS) / 2 "Mean enthalpy of the superheated region";
    TwA = (TwSB * L_TP + TwTP * L_SB) / (L_TP + L_SB) "Wall Temperature at the border between subcooled and two-phase region";
    TwB = (TwSH * L_TP + TwTP * L_SH) / (L_TP + L_SH) "Wall Temperature at the border between two-phase and superheated region";
    Tsf_SH = (Tsf_SU + Tsf_B) / 2 "Mean temperature of the secondary fluid in the superheated region";
// Tsf_B = (L_TP*Tsf_SH + L_SH*Tsf_TP)/(L_SH+L_TP)
//     "Secondary fluid temperature at the border between two-phase and superheated region";
    Tsf_TP = (Tsf_A + Tsf_B) / 2 "Mean temperature of the secondary fluid in the two phase region";
// Tsf_A = (L_TP*Tsf_SB + L_SB*Tsf_TP)/(L_SB+L_TP)
//     "Secondary fluid temperature at the border between two-phase and subcooled";
    Tsf_SB = (Tsf_A + Tsf_EX) / 2 "Mean temperature of the secondary fluid in the two phase region";
/**** SECONDARY FLUID - WALL HEAT TRANSFER ****/
/** Subcooled zone **/
    NTU_SB_sf = Usf * pi * D * L_SB / (cp_sf * Mdot_sf);
    epsilon_SB_sf = 1 - exp(-NTU_SB_sf);
    Qsf_SB = cp_sf * Mdot_sf * epsilon_SB_sf * (Tsf_A - TwSB);
/** Two-phase zone **/
    NTU_TP_sf = Usf * pi * D * L_TP / (cp_sf * Mdot_sf);
    epsilon_TP_sf = 1 - exp(-NTU_TP_sf);
    Qsf_TP = cp_sf * Mdot_sf * epsilon_TP_sf * (Tsf_B - TwTP);
/** SuperHeated zone **/
    NTU_SH_sf = Usf * pi * D * L_SH / (cp_sf * Mdot_sf);
    epsilon_SH_sf = 1 - exp(-NTU_SH_sf);
    Qsf_SH = cp_sf * Mdot_sf * epsilon_SH_sf * (Tsf_SU - TwSH);
    Qsf_SB = Mdot_sf * cp_sf * (Tsf_A - Tsf_EX);
    Qsf_TP = Mdot_sf * cp_sf * (Tsf_B - Tsf_A);
    Qsf_SH = Mdot_sf * cp_sf * (Tsf_SU - Tsf_B);
/*** WALL - PRIMARY FLUID HEAT TRANSFER ***/
/** Sub-cooled zone **/
    NTU_SB_pf = U_SB * pi * D * L_SB / (cp_SB * (M_dot_su + M_dot_A) / 2);
    epsilon_SB_pf = 1 - exp(-NTU_SB_pf);
    Q_SB = cp_SB * (M_dot_su + M_dot_A) / 2 * epsilon_SB_pf * (TwSB - T_SU) "Heat from the wall sub-cooled region";
/** Two-phase zone **/
/** I say that Epsilon NTU in the TP zone is useless **/
    NTU_TP_pf = U_TP * pi * D * L_TP / (cp_TP * (M_dot_A + M_dot_B) / 2);
    epsilon_TP_pf = 1 - exp(-NTU_TP_pf);
    Q_TP = pi * D * L_TP * U_TP * (TwTP - T_TP) "Heat from the wall two-phase region";
/** SuperHeated zone **/
    NTU_SH_pf = U_SH * pi * D * L_SH / (cp_SH * (M_dot_B + M_dot_ex) / 2);
    epsilon_SH_pf = 1 - exp(-NTU_SH_pf);
    Q_SH = cp_SH * (M_dot_B + M_dot_ex) / 2 * (TwSH - T_TP) "Heat from the wall super-heated region";
//* epsilon_SH_pf
/** Total energy exchange from the primary fluid and the secondary fluid**/
    Q_tot = Q_SB + Q_TP + Q_SH;
    Q_tot_sf = Qsf_SB + Qsf_TP + Qsf_SH;
/************ FLUID PROPERTIES ************************/
    SubCooled = Medium.setState_ph(p, h_SB);
    sat = Medium.setSat_p(p);
    SuperHeated = Medium.setState_ph(p, h_SH);
    outlet = Medium.setState_ph(p, h_EX);
/**Temperatures **/
    T_SU = Medium.temperature_ph(p, h_SU);
    T_SB = Medium.temperature(SubCooled);
    T_TP = Medium.saturationTemperature_sat(sat);
    T_SH = Medium.temperature(SuperHeated);
    T_out = Medium.temperature(outlet);
/** Densities **/
    rho_SB = Medium.density(SubCooled);
    rho_LS = Medium.bubbleDensity(sat);
    rho_VS = Medium.dewDensity(sat);
    rho_SH = Medium.density(SuperHeated);
/** SpecificHeatCapacity **/
    cp_SB = Medium.specificHeatCapacityCp(SubCooled);
    cp_TP = (cp_SB + cp_SH) / 2;
    cp_SH = Medium.specificHeatCapacityCp(SuperHeated);
/** Enthalpies **/
    h_LS = Medium.bubbleEnthalpy(sat);
    h_VS = Medium.dewEnthalpy(sat);
/**Derivatives**/
    drdh_SB = Medium.density_derh_p(SubCooled);
    drdp_SB = Medium.density_derp_h(SubCooled);
    drdp_LS = Medium.dBubbleDensity_dPressure(sat);
    drdp_VS = Medium.dDewDensity_dPressure(sat);
    dhdp_LS = Medium.dBubbleEnthalpy_dPressure(sat);
    dhdp_VS = Medium.dDewEnthalpy_dPressure(sat);
    drdh_SH = Medium.density_derh_p(SuperHeated);
    drdp_SH = Medium.density_derp_h(SuperHeated);
/************ MassBalance SUB-COOLED ****************/
    A * (L_SB * drdt_SB + (rho_SB - rho_LS) * der(L_SB)) = M_dot_su - M_dot_A;
    drdt_SB = (drdp_SB + 1 / 2 * drdh_SB * dhdp_LS) * der(p) + 1 / 2 * drdh_SB * der(h_SU);
/********** EnergyBalance SUB-COOLED ******************/
    A * L_SB * (rho_SB * dhdt_SB + h_SB * drdt_SB - der(p)) + A * (rho_SB * h_SB - rho_LS * h_LS) * der(L_SB) = M_dot_su * h_SU - M_dot_A * h_LS + Q_SB;
    dhdt_SB = 1 / 2 * (dhdp_LS * der(p) + der(h_SU));
/*******  WallEnergyBalance SUB-COOLED **********/
    c_wall * M_tot / L * (L_SB * der(TwSB) + (TwSB - TwA) * der(L_SB)) = Qsf_SB - Q_SB;
/**************************************************************************************************************/
/************ MassBalance TWO-PHASE REGION ******************/
    A * (L_TP * drdt_TP + (rho_TP - rho_VS) * der(L_TP) + (rho_LS - rho_VS) * der(L_SB)) = M_dot_A - M_dot_B;
    drdt_TP = (Void * drdp_VS + (1 - Void) * drdp_LS) * der(p);
//(rho_VS - rho_LS)*(dVoid_dp*der(p) + dVoid_dh*der(h_EX));
/************ EnergyBalance TWO-PHASE REGION ******************/
    A * (L_TP * drhdt_TP + (rhoh_TP - rho_VS * h_VS) * der(L_TP) + (rho_LS * h_LS - rho_VS * h_VS) * der(L_SB) - L_TP * der(p)) = M_dot_A * h_LS - M_dot_B * h_VS + Q_TP;
    drhdt_TP = (Void * (drdp_VS * h_VS + rho_VS * drdp_VS) + (1 - Void) * (h_LS * drdp_LS + rho_LS * dhdp_LS)) * der(p);
//(rho_VS*h_VS -rho_LS*h_LS)*(dVoid_dp*der(p) + dVoid_dh*der(h_EX));
/***************** WallEnergyBalance TWO-PHASE REGION **************/
    c_wall * M_tot / L * (L_TP * der(TwTP) + (TwA - TwB) * der(L_SB) + (TwTP - TwB) * der(L_TP)) = Qsf_TP - Q_TP;
/**************************************************************************************************************/
/************ MassBalance SUPER-HEATED ****************/
    A * (L_SH * drdt_SH + (rho_VS - rho_SH) * der(L_SH)) = M_dot_B - M_dot_ex;
    drdt_SH = (drdp_SH + 1 / 2 * drdh_SH * dhdp_VS) * der(p) + 1 / 2 * drdh_SH * der(h_EX);
/********** EnergyBalance SUPER-HEATED ******************/
    A * L_SH * (rho_SH * dhdt_SH + h_SH * drdt_SH - der(p)) + A * (rho_VS * h_VS - rho_SH * h_SH) * der(L_SH) = M_dot_B * h_VS - M_dot_ex * h_EX + Q_SH;
    dhdt_SH = 1 / 2 * (dhdp_VS * der(p) + der(h_EX));
/*******  WallEnergyBalance SUPER-HEATED **********/
    c_wall * M_tot / L * (L_SH * der(TwSH) + (TwB - TwSH) * der(L_SH)) = Qsf_SH - Q_SH;
/******************************* SUMMARY ***********************************/
  protected
    parameter Integer N = 3;
    Modelica.SIunits.Temperature[N] T_fluid_;
    Modelica.SIunits.Temperature[N] T_wall_;
    Modelica.SIunits.Temperature[N] T_sf_;
  public
    record SummaryBase
      replaceable Arrays T_profile;

      record Arrays
        parameter Integer n;
        Modelica.SIunits.Temperature[n] Twf;
        Modelica.SIunits.Temperature[n] Twall;
        Modelica.SIunits.Temperature[n] Tsf;
      end Arrays;

      //Modelica.SIunits.Pressure p_sf;
      Modelica.SIunits.Pressure p_wf;
      //Modelica.SIunits.Power Q_sf;
      //Modelica.SIunits.Power Q_wf;
    end SummaryBase;

    replaceable record SummaryClass = SummaryBase;
    SummaryClass Summary(T_profile(n = N, Twf = T_fluid_[:], Twall = T_wall_[:], Tsf = T_sf_[:]), p_wf = p);
  equation
/* Fluid temperature */
    T_fluid_[1] = T_SU;
    T_fluid_[2] = T_TP;
    T_fluid_[3] = T_out;
/* Wall temperature */
    T_wall_[1] = TwSB;
    T_wall_[2] = TwTP;
    T_wall_[3] = TwSH;
/* Secondary fluid temperature */
    T_sf_[1] = Tsf_SU;
    T_sf_[2] = Tsf_TP;
    T_sf_[3] = Tsf_EX;
    annotation(
      Diagram(graphics),
      Icon(graphics = {Rectangle(extent = {{-100, -20}, {-32, -100}}, lineColor = {60, 121, 182}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-100, 20}, {100, -20}}, lineColor = {135, 135, 135}, fillPattern = FillPattern.CrossDiag), Rectangle(extent = {{-100, 100}, {100, 20}}, lineColor = {135, 135, 135}, fillColor = {255, 85, 85}, fillPattern = FillPattern.Solid), Line(points = {{-100, 56}, {-80, 36}, {-60, 56}, {-40, 36}, {-20, 56}, {0, 36}, {20, 56}, {40, 36}, {60, 56}, {80, 36}, {100, 56}}, color = {255, 0, 0}, smooth = Smooth.None, thickness = 0.5), Polygon(points = {{-18, 84}, {-18, 64}, {-32, 74}, {-18, 84}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-18, 74}, {34, 74}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{-18, -72}, {-10, -80}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{-100, -100}, {40, -100}, {-16, -20}, {-100, -20}, {-100, -24}, {-100, -100}}, lineColor = {0, 128, 255}, smooth = Smooth.None, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{32, -100}, {44, -100}, {102, -100}, {100, -20}, {-16, -20}, {32, -100}}, lineColor = {170, 255, 255}, smooth = Smooth.None, fillColor = {170, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{-97, -63}, {-75, -43}, {-57, -63}, {-37, -43}, {-17, -63}, {3, -43}, {23, -63}, {43, -43}, {63, -63}, {83, -43}, {103, -63}}, color = {0, 0, 255}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{4, -30}, {12, -38}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-2, -54}, {6, -62}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{42, -70}, {42, -90}, {56, -80}, {42, -70}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-8, -80}, {44, -80}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Ellipse(extent = {{28, -28}, {34, -34}}, lineColor = {85, 170, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-14, -26}, {-2, -40}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{12, -64}, {22, -76}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{26, -78}, {38, -90}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{16, -42}, {28, -54}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{36, -62}, {42, -68}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid)}),
      Documentation(info = "<HTML>
          
         <p><big>Model <b>MBeva</b> is a moving boundaries heat exchanger model for evaporators.
         The model considers a fictitius heat transfer channel split up into three different section based on the phase state of the working fluid:
         <ul><li> Sub-cooled zone (SB)
         <li> Two-phase zone (TP)
         <li> SuperHeated zone (SH)
         </ul>
          <p><big><b>Pressure</b> and <b>enthalpy</b> are selected as state variables for each of the three zones. 
          
        <p><big>  The name moving boundary is derived from the fact that the interfaces between these sections do not
have a fixed spatial position but merely a fixed thermodynamic location depending
on the presence of liquid and gaseous fluid, respectively. The actual existence of
a certain section and its length are determined based on the fluid state resulting
in variable sectioning. A fixed total length superimposes the required boundary
condition to calculate the length of each section.
 
 <p><big>The assumptions for this model are:
         <ul><li> The tube is cylindrical with a constant cross sectional area
         <li> The velocity of the fluid is uniform on the cross sectional area
         <li> The enthalpy of the fluid is linear in each region of the tube (sub-cooled, two-phase, super-heated)
         <li> Pressure is considered constant
         <li> Constant void fraction
         <li> Thermal energy accumulation in the metal wall is taken into account
         <li> The secondary fluid is treated as a constant heat capacity fluid
         <li> The heat transfer between the secondary fluid and the metal wall and between the metal wall and the primary fluid is computed using the epsilon-NTU method
         </ul>
        </HTML>"));
  end HEX_3f;


  model HEX_3f_ss "Counter current STEADY-STATE Moving Boundary model: Fluid enters subcooled and exits in super-heated conditons. The model consider the fluid in one side and the metal wall. The secondary fluid is a constant specific heat fluid"
    /**************** MEDIUM ***************************/
    replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialMedium annotation(
       choicesAllMatching = true);
    /***************** PORTS ******************/
    Modelica.Fluid.Interfaces.FluidPort_a InFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-108, -74}, {-88, -54}})));
    Modelica.Fluid.Interfaces.FluidPort_b OutFlow(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{90, -74}, {110, -54}})));
    ShipEnergySystemsModelicaLibrary.Other.FlangeCdot_B OutFlow_sf annotation(
      Placement(transformation(extent = {{-110, 50}, {-90, 70}})));
    ShipEnergySystemsModelicaLibrary.Other.FlangeCdot_A InFlow_sf annotation(
      Placement(transformation(extent = {{90, 50}, {110, 70}})));
    /***********THERMODYNAMIC STATES *********************/
    Medium.ThermodynamicState SubCooled "Subcooled thermodynamic state";
    Medium.SaturationProperties sat;
    Medium.ThermodynamicState SuperHeated "Superheated thermodynamic state";
    /************* Tube pressure *****************/
    Modelica.SIunits.Pressure p(start = 0) "Pressure in the heat-exchanger: No pressure drop considered";
    /*********** GENERAL VARIABLES AND PARAMETER ***********/
    parameter Modelica.SIunits.Area A "Cross-sectional area for both side of the heat exchanger";
    parameter Modelica.SIunits.Length L "Total length of the exchanger";
    final parameter Modelica.SIunits.Length D = 2 * sqrt(A / pi) "Diameter";
    constant Real pi = Modelica.Constants.pi;
  
    /********** MASS FLOWS *******************/
    //parameter Modelica.SIunits.MassFlowRate Mdotnom
    //    "Nominal fluid flow rate of the primary fluid";
    Modelica.SIunits.MassFlowRate Mdot_wf "Mass flow of the working fluid of the heat exchanger";
    Modelica.SIunits.MassFlowRate Mdot_sf "Mass flow of the secondary fluid of the heat exchanger";
    /*****************HEAT FLOW *************************/
    Modelica.SIunits.Power Q_SB "Thermal energy transfer with the wall";
    Modelica.SIunits.Power Q_TP "Thermal energy transfer with the wall in the two-phase region";
    Modelica.SIunits.Power Q_SH "Thermal energy transfer with the wall in the superheated region";
    Modelica.SIunits.Power Q_tot "Total Thermal energy transfer with the wall from the primary fluid";
  
    /***************** TEMPERATURES *************************/
    /** primary fluid **/
    Modelica.SIunits.Temperature T_SU "Temperatue at the inlet of the heat exchanger";
    Modelica.SIunits.Temperature T_TP "Mean Temperatue of the two-phase region";
    Modelica.SIunits.Temperature T_out "Outlet Temperature of the heat exchanger, working fluid side";
    
    /** Secondary fluid **/
    Modelica.SIunits.Temperature Tsf_SU "Tempearture at the inlet of the secondary fluid";
    Modelica.SIunits.Temperature Tsf_B "Secondary fluid temperature at the interface between two-phase and super-heated";
    Modelica.SIunits.Temperature Tsf_A "Secondary fluid temperature at the interface between subcooled and two-phase";
    Modelica.SIunits.Temperature Tsf_EX "Secondary fluid temperature at the exit of the heat exchanger";
    /*************** HEAT TRANSFER PARAMETER ****************/
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_wf_SB "Working fluid Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_wf_TP "Working fluid Heat transfer coefficient two-phase side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_wf_SH "Working fluid Heat transfer coefficient Superheated side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_sf "Secondary fluid Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_SB "Total Heat transfer coefficient sub-cooled side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_TP "Total Heat transfer coefficient two-phase side";
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U_SH "Total Heat transfer coefficient Superheated side"; 
    /*************Tube region length ********************/
    Modelica.SIunits.Length L_SB(start = 0) "Length of the subcooled region";
    Modelica.SIunits.Length L_TP(start = 0) "Length of the two-phase region";
    Modelica.SIunits.Length L_SH "Length of the Superheated region";
  
    /********** Specific heat capacity ***************/
    Modelica.SIunits.SpecificHeatCapacity cp_sf "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_SB "Specific heat capacity of the secondary fluid at the inlet";
    Modelica.SIunits.SpecificHeatCapacity cp_SH "Specific heat capacity of the secondary fluid at the inlet";
    /************** Enthalpies ************************/
    Modelica.SIunits.SpecificEnthalpy h_SU "Enthalpy at the inlet of the heat exchanger of the primary fluid";
    Modelica.SIunits.Temperature h_sf_A "Primary fluid enthalpy at the interface between subcooled and two-phase";
    Modelica.SIunits.SpecificEnthalpy h_sf_SB "Subcooled region average enthalpy";
    Modelica.SIunits.SpecificEnthalpy h_sf_LS "Bubble enthalpy";
    Modelica.SIunits.Temperature h_sf_B "Primary fluid enthalpy at the interface between two-phase and super-heated";
    Modelica.SIunits.SpecificEnthalpy h_sf_VS "Dew enthalpy";
    Modelica.SIunits.SpecificEnthalpy h_sf_SH "Superheatd region average enthalpy";
    
    Modelica.SIunits.SpecificEnthalpy h_sf_EX "Enthalpy at the outlet of the heat exchanger of the primary fluid";
    
    /***** Enthalpy&Density****/
    Modelica.SIunits.Density rho_sf "Density of the secondary fluid";
    Real MC_min_SB "Minimum Mdot*cp in the subcooled region";
    Real MC_min_SH "Minimum Mdot*cp in the superheated region";
  
    /********************** VOID fraction *****************************/
    
    /********************** Epsilon-NTU method *********************/
    Real epsilon_SB "Epsilon Sub-cooled secondary fluid - wall Side";
    Real NTU_SB "NTU Sub-cooled secondary fluid-Wall side";
    Real epsilon_TP "Epsilon Two-Phase secondary fluid - wall Side";
    Real NTU_TP "NTU Two-Phase secondary fluid-Wall side";
    Real epsilon_SH "Epsilon Super-heated secondary fluid - wall Side";
    Real NTU_SH "NTU Super-heated secondary fluid-Wall side";
  
  equation
  /****************** BOUNDARY EQUATION *************************/
    L = L_SB + L_TP + L_SH;
  /********** Enthalpies ************/
    h_SU = inStream(InFlow.h_outflow);
    h_EX = OutFlow.h_outflow;
    h_SU = InFlow.h_outflow;
  /*********** Temperatures *********/
    OutFlow_sf.T = Tsf_EX;
    InFlow_sf.T = Tsf_SU;
  /**********Specific Heat Capacity *************/
    OutFlow_sf.cp = cp_sf;
    InFlow_sf.cp = cp_sf;
  /*********** Density ***********/
    OutFlow_sf.rho = rho_sf;
    InFlow_sf.rho = rho_sf;
  /*************** Pressure **************/
    InFlow.p = p;
    OutFlow.p = p;
  /***********Mass Flow ******************/
    InFlow.m_flow = Mdot_wf;
    OutFlow.m_flow = -Mdot_wf;
    OutFlow_sf.Mdot = Mdot_sf;
    InFlow_sf.Mdot = Mdot_sf;
    
  /**** EPSILON - NTU method ****/
  /** Subcooled zone **/
    U_SB = 1/(1/U_wf_SB + 1/U_sf);
    if cp_sf * Mdot_sf > Mdot_wf * cp_SB then
      MC_min_SB = Mdot_wf * cp_SB;
    else
      MC_min_SB = cp_sf * Mdot_sf;
    end if;
    NTU_SB = U_SB * pi * D * L_SB / MC_min_SB;
    epsilon_SB = 1 - exp(-NTU_SB);
    Q_SB = MC_min_SB * epsilon_SB * (Tsf_A - T_SU);
    
  /** Two-phase zone **/
    U_TP = 1/(1/U_wf_TP + 1/U_sf);
    NTU_TP = U_TP * pi * D * L_TP / (cp_sf * Mdot_sf);
    epsilon_TP = 1 - exp(-NTU_TP);
    Q_TP = cp_sf * Mdot_sf * epsilon_TP * (Tsf_B - T_TP);
  
  /** SuperHeated zone **/
    U_SH = 1/(1/U_wf_SH + 1/U_sf);
    if cp_sf * Mdot_sf > Mdot_wf * cp_SH then
      MC_min_SH = Mdot_wf * cp_SH;
    else
      MC_min_SH = cp_sf * Mdot_sf;
    end if;
    NTU_SH = U_SH * pi * D * L_SH / MC_min_SH;
    epsilon_SH = 1 - exp(-NTU_SH);
    Q_SH = MC_min_SH * epsilon_SH * (Tsf_SU - T_TP);
    
    // Average enthalpies in the subcooled and superheatd regions //
    h_SB = (h_wf_SU + h_wf_A)/2 ;
    h_SH = (h_wf_EX + h_wf_VS)/2 ;
    
  
  /************ FLUID PROPERTIES ************************/
    sat = Medium.setSat_p(p);
    SubCooled = Medium.setState_ph(p, h_wf_SB);
    SuperHeated = Medium.setState_ph(p, h_wf_SH);
    //outlet = Medium.setState_ph(p, h_EX);
  /**Temperatures **/
    T_SU = Medium.temperature_ph(p, h_wf_SU);
    T_TP = Medium.saturationTemperature_sat(sat);
    if L_SH > 0 then
      T_out = Medium.temperature_ph(p, h_wf_VS + Q_SH/Mdot_wf);
    elseif L_TP > 0 then
      T_out = T_TP;
    else
      T_out = Medium.temperature_ph(p, h_wf_SU + Q_SB/Mdot_wf);
    end if;
  /** SpecificHeatCapacity **/
    cp_SB = Medium.specificHeatCapacityCp(SubCooled);
    cp_SH = Medium.specificHeatCapacityCp(SuperHeated);
  /** Enthalpies **/
    h_LS = Medium.bubbleEnthalpy(sat);
    h_VS = Medium.dewEnthalpy(sat);
    
  /************ Energy balances ****************/
  0 = Mdot_wf * (h_wf_SU - h_wf_A) + Q_SB; 
  0 = Mdot_wf * (h_wf_LS - h_wf_VS) + Q_TP;
  0 = Mdot_wf * (h_wf_VS - h_wf_EX) + Q_SH;
  Q_SB = Mdot_sf * cp_sf * (T_sf_A - T_sf_EX);
  Q_TP = Mdot_sf * cp_sf * (T_sf_B - T_sf_A);
  Q_SH = Mdot_sf * cp_sf * (T_sf_SU - T_sf_B);
  
  end HEX_3f_ss;
  
  
  
  
  
  model HEX_tf_Uk "Counter current transfer-function (tf) type heat exchanger model. The idea is to simply model the heat exchanger as a transfer function. A constant U (Uk) is used"
  /**************** MEDIUM ***************************/
  replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialMedium annotation(
     choicesAllMatching = true);
  /***************** PORTS ******************/
  Modelica.Fluid.Interfaces.FluidPort_a InFlow(redeclare package Medium = Medium) annotation(
    Placement(transformation(extent = {{-108, -74}, {-88, -54}})));
  Modelica.Fluid.Interfaces.FluidPort_b OutFlow(redeclare package Medium = Medium) annotation(
    Placement(transformation(extent = {{90, -74}, {110, -54}})));
  ShipEnergySystemsModelicaLibrary.Other.FlangeCdot_B OutFlow_sf annotation(
    Placement(transformation(extent = {{-110, 50}, {-90, 70}})));
  ShipEnergySystemsModelicaLibrary.Other.FlangeCdot_A InFlow_sf annotation(
    Placement(transformation(extent = {{90, 50}, {110, 70}})));
  /***********THERMODYNAMIC STATES *********************/
  /************* Tube pressure *****************/
  Modelica.SIunits.Pressure p(start = 0) "Pressure in the heat-exchanger: No pressure drop considered";
  /*********** GENERAL VARIABLES AND PARAMETER ***********/
  parameter Modelica.SIunits.Area A "Cross-sectional area for both side of the heat exchanger";
  /********** MASS FLOWS *******************/
  Modelica.SIunits.MassFlowRate Mdot_wf "Mass flow of the working fluid of the heat exchanger";
  Modelica.SIunits.MassFlowRate Mdot_sf "Mass flow of the secondary fluid of the heat exchanger";
  /*****************HEAT FLOW *************************/
  Modelica.SIunits.Power Q_wf "Thermal energy transfer from the primary fluid";
  Modelica.SIunits.Power Q_sf "Thermal energy transfer from the secondary fluid";

  /***************** TEMPERATURES *************************/
  Modelica.SIunits.Temperature dT_ml "Mean logarithmic temperature difference";
  /** primary fluid **/
  Modelica.SIunits.Temperature T_wf_in "Temperatue at the inlet of the heat exchanger";
  Modelica.SIunits.Temperature T_wf "Mean Temperatue of the working fluid";
  Modelica.SIunits.Temperature T_wf_out "Outlet Temperature of the heat exchanger, working fluid side";
  /** Secondary fluid **/
  Modelica.SIunits.Temperature Tsf_in "Tempearture at the inlet of the secondary fluid";
  Modelica.SIunits.Temperature Tsf_B "Mean Secondary fluid temperature";
  Modelica.SIunits.Temperature Tsf_out "Secondary fluid temperature at the exit of the heat exchanger";
  /*************** HEAT TRANSFER PARAMETER ****************/
  parameter Modelica.SIunits.CoefficientOfHeatTransfer U_wf "Working fluid Heat transfer coefficient";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer U_sf "Secondary fluid Heat transfer coefficient sub-cooled side";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer U_tot "Total Heat transfer coefficient";
  
  /********** Specific heat capacity ***************/
  Modelica.SIunits.SpecificHeatCapacity cp_sf "Specific heat capacity of the secondary";
  Modelica.SIunits.SpecificHeatCapacity cp_wf "Specific heat capacity of the working fluid";

  /************** Enthalpies ************************/
  Modelica.SIunits.SpecificEnthalpy h_wf_in "Enthalpy at the inlet of the heat exchanger of the primary fluid";
  Modelica.SIunits.SpecificEnthalpy h_wf "Mean Enthalpy of the heat exchanger of the primary fluid";
  Modelica.SIunits.SpecificEnthalpy h_wf_out "Enthalpy at the outlet of the heat exchanger of the primary fluid";
  Modelica.SIunits.SpecificEnthalpy h_sf_in "Enthalpy at the inlet of the heat exchanger of the secondary fluid";
  Modelica.SIunits.SpecificEnthalpy h_sf "Mean Enthalpy of the heat exchanger of the secondary fluid";
  Modelica.SIunits.SpecificEnthalpy h_sf_out "Enthalpy at the outlet of the heat exchanger of the secondary fluid";
  
  
  /***** Enthalpy&Density****/
  Modelica.SIunits.Density rho_sf "Density of the secondary fluid";
  Real MC_min "Minimum Mdot*cp ";

equation

/********** Enthalpies ************/
  h_wf_in = inStream(InFlow.h_outflow);
  h_wf = OutFlow.h_outflow;
  h_wf_out = InFlow.h_outflow;
/*********** Temperatures *********/
  OutFlow_sf.T = T_sf_out;
  InFlow_sf.T = T_sf_in;
/**********Specific Heat Capacity *************/
  OutFlow_sf.cp = cp_sf;
  InFlow_sf.cp = cp_sf;
/*********** Density ***********/
  OutFlow_sf.rho = rho_sf;
  InFlow_sf.rho = rho_sf;
/*************** Pressure **************/
  InFlow.p = p;
  OutFlow.p = p;
/***********Mass Flow ******************/
  InFlow.m_flow = Mdot_wf;
  OutFlow.m_flow = -Mdot_wf;
  OutFlow_sf.Mdot = Mdot_sf;
  InFlow_sf.Mdot = Mdot_sf;
  
/*********** Q = UA dT_ml ********/
  Qdot_sf = U_sf * A * dT_ml;
  Qdot_wf = U_wf * A * dT_ml;
  TC * der(dT_ml) = Qdot_sf * Qdot_wf;

/************ FLUID PROPERTIES ************************/
  sat = Medium.setSat_p(p);
  SubCooled = Medium.setState_ph(p, h_wf_SB);
  SuperHeated = Medium.setState_ph(p, h_wf_SH);
  //outlet = Medium.setState_ph(p, h_EX);
/**Temperatures **/
  T_SU = Medium.temperature_ph(p, h_wf_SU);
  T_TP = Medium.saturationTemperature_sat(sat);
  if L_SH > 0 then
    T_out = Medium.temperature_ph(p, h_wf_VS + Q_SH/Mdot_wf);
  elseif L_TP > 0 then
    T_out = T_TP;
  else
    T_out = Medium.temperature_ph(p, h_wf_SU + Q_SB/Mdot_wf);
  end if;
/** SpecificHeatCapacity **/
  cp_SB = Medium.specificHeatCapacityCp(SubCooled);
  cp_SH = Medium.specificHeatCapacityCp(SuperHeated);
/** Enthalpies **/
  h_LS = Medium.bubbleEnthalpy(sat);
  h_VS = Medium.dewEnthalpy(sat);
  
/************ Energy balances ****************/
0 = Mdot_wf * (h_wf_SU - h_wf_A) + Q_SB; 
0 = Mdot_wf * (h_wf_LS - h_wf_VS) + Q_TP;
0 = Mdot_wf * (h_wf_VS - h_wf_EX) + Q_SH;
Q_SB = Mdot_sf * cp_sf * (T_sf_A - T_sf_EX);
Q_TP = Mdot_sf * cp_sf * (T_sf_B - T_sf_A);
Q_SH = Mdot_sf * cp_sf * (T_sf_SU - T_sf_B);

end HEX_tf_Uk;














































































end HeatExchanger;
