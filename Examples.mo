within ShipEnergySystemsModelicaLibrary;

package Examples
  model HullAndpropellerTest
    DCmotor electricMotor annotation(
      Placement(visible = true, transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp InputSignal(duration = 200, height = 1000) annotation(
      Placement(visible = true, transformation(origin = {-64, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Components.Fixed fixed1 annotation(
      Placement(visible = true, transformation(origin = {86, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShipEnergySystemsModelicaLibrary.Propeller.FixedEfficiency Propeller(K_Q = 0.035, K_T = 0.25, diameter = 2, useSupportR = false, useSupportT = false) annotation(
      Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShipHull.FixedC Hull(C1 = 100) annotation(
      Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Hull.flange_water, fixed1.flange) annotation(
      Line(points = {{76, 0}, {86, 0}, {86, -54}, {86, -54}}, color = {0, 127, 0}));
    connect(Propeller.flangeT, Hull.flange_prop) annotation(
      Line(points = {{36, 0}, {56, 0}, {56, 0}, {56, 0}}, color = {0, 127, 0}));
    connect(Propeller.flangeR, electricMotor.flange_a) annotation(
      Line(points = {{16, 0}, {-4, 0}, {-4, 0}, {-4, 0}}));
    connect(InputSignal.y, electricMotor.u) annotation(
      Line(points = {{-52, 0}, {-22, 0}, {-22, 0}, {-22, 0}}, color = {0, 0, 127}));
    annotation(
      uses(Modelica(version = "3.2.2")));
  end HullAndpropellerTest;

  model baseTest
    DCmotor electricMotor annotation(
      Placement(visible = true, transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp InputSignal annotation(
      Placement(visible = true, transformation(origin = {-64, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Components.IdealRollingWheel idealRollingWheel1 annotation(
      Placement(visible = true, transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(electricMotor.flange_a, idealRollingWheel1.flangeR) annotation(
      Line(points = {{-4, 0}, {20, 0}, {20, 0}, {20, 0}}));
    connect(InputSignal.y, electricMotor.u) annotation(
      Line(points = {{-52, 0}, {-22, 0}, {-22, 0}, {-22, 0}}, color = {0, 0, 127}));
    annotation(
      uses(Modelica(version = "3.2.2")));
  end baseTest;

  model propellerTest
    DCmotor electricMotor annotation(
      Placement(visible = true, transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp InputSignal(duration = 200, height = 1000) annotation(
      Placement(visible = true, transformation(origin = {-64, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Components.Mass Inertia(m = 1000) annotation(
      Placement(visible = true, transformation(origin = {58, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Components.Fixed fixed1 annotation(
      Placement(visible = true, transformation(origin = {78, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShipEnergySystemsModelicaLibrary.Propeller.FixedEfficiency Propeller(K_Q = 0.035, K_T = 0.25, diameter = 2, useSupportR = false, useSupportT = false) annotation(
      Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Components.Damper damper1(d = 1000) annotation(
      Placement(visible = true, transformation(origin = {78, -22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  equation
    connect(damper1.flange_b, fixed1.flange) annotation(
      Line(points = {{78, -32}, {78, -32}, {78, -54}, {78, -54}}, color = {0, 127, 0}));
    connect(Inertia.flange_b, damper1.flange_a) annotation(
      Line(points = {{68, 0}, {78, 0}, {78, -12}, {78, -12}}, color = {0, 127, 0}));
    connect(Inertia.flange_a, Propeller.flangeT) annotation(
      Line(points = {{48, 0}, {36, 0}}, color = {0, 127, 0}));
    connect(Propeller.flangeR, electricMotor.flange_a) annotation(
      Line(points = {{16, 0}, {-4, 0}, {-4, 0}, {-4, 0}}));
    connect(InputSignal.y, electricMotor.u) annotation(
      Line(points = {{-52, 0}, {-22, 0}, {-22, 0}, {-22, 0}}, color = {0, 0, 127}));
    annotation(
      uses(Modelica(version = "3.2.2")));
  end propellerTest;

  model SpeedControlTest
    DCmotor electricMotor annotation(
      Placement(visible = true, transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp InputSignal(duration = 200, height = 6) annotation(
      Placement(visible = true, transformation(origin = {-86, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Components.Fixed fixed1 annotation(
      Placement(visible = true, transformation(origin = {86, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShipEnergySystemsModelicaLibrary.Propeller.FixedEfficiency Propeller(K_Q = 0.035, K_T = 0.25, diameter = 2, useSupportR = false, useSupportT = false) annotation(
      Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShipHull.FixedC Hull(C1 = 100) annotation(
      Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor1 annotation(
      Placement(visible = true, transformation(origin = {32, -42}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Continuous.PID PID(Ti = 10, k = 1000) annotation(
      Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Feedback feedback1 annotation(
      Placement(visible = true, transformation(origin = {-40, -46}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
    Modelica.Blocks.Math.Gain gain1(k = -1) annotation(
      Placement(visible = true, transformation(origin = {2, -42}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
  equation
    connect(gain1.y, feedback1.u2) annotation(
      Line(points = {{-4, -42}, {-14, -42}, {-14, -28}, {-40, -28}, {-40, -38}, {-40, -38}}, color = {0, 0, 127}));
    connect(speedSensor1.v, gain1.u) annotation(
      Line(points = {{20, -42}, {10, -42}, {10, -42}, {10, -42}, {10, -42}}, color = {0, 0, 127}));
    connect(PID.y, electricMotor.u) annotation(
      Line(points = {{-40, 0}, {-22, 0}, {-22, 0}, {-22, 0}}, color = {0, 0, 127}));
    connect(feedback1.y, PID.u) annotation(
      Line(points = {{-30, -46}, {-26, -46}, {-26, -16}, {-80, -16}, {-80, 0}, {-64, 0}, {-64, 0}}, color = {0, 0, 127}));
    connect(InputSignal.y, feedback1.u1) annotation(
      Line(points = {{-74, -46}, {-48, -46}, {-48, -46}, {-48, -46}}, color = {0, 0, 127}));
    connect(speedSensor1.flange, Hull.flange_prop) annotation(
      Line(points = {{42, -42}, {56, -42}, {56, 0}, {56, 0}}, color = {0, 127, 0}));
    connect(Hull.flange_water, fixed1.flange) annotation(
      Line(points = {{76, 0}, {86, 0}, {86, -54}, {86, -54}}, color = {0, 127, 0}));
    connect(Propeller.flangeT, Hull.flange_prop) annotation(
      Line(points = {{36, 0}, {56, 0}, {56, 0}, {56, 0}}, color = {0, 127, 0}));
    connect(Propeller.flangeR, electricMotor.flange_a) annotation(
      Line(points = {{16, 0}, {-4, 0}, {-4, 0}, {-4, 0}}));
    annotation(
      uses(Modelica(version = "3.2.2")));
  end SpeedControlTest;

  model EngineControlTest
    Modelica.Blocks.Sources.Ramp InputSignal(duration = 200, height = 6) annotation(
      Placement(visible = true, transformation(origin = {-86, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Components.Fixed fixed1 annotation(
      Placement(visible = true, transformation(origin = {86, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShipEnergySystemsModelicaLibrary.Propeller.FixedEfficiency Propeller(K_Q = 0.035, K_T = 0.25, diameter = 2, useSupportR = false, useSupportT = false) annotation(
      Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShipHull.FixedC Hull(C1 = 100) annotation(
      Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor1 annotation(
      Placement(visible = true, transformation(origin = {32, -42}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Continuous.PID PID(k = 0.5) annotation(
      Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Feedback feedback1 annotation(
      Placement(visible = true, transformation(origin = {-40, -46}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
    Modelica.Blocks.Math.Gain gain1(k = 1) annotation(
      Placement(visible = true, transformation(origin = {2, -42}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
  equation
    connect(gain1.y, feedback1.u2) annotation(
      Line(points = {{-4, -42}, {-14, -42}, {-14, -28}, {-40, -28}, {-40, -38}, {-40, -38}}, color = {0, 0, 127}));
    connect(speedSensor1.v, gain1.u) annotation(
      Line(points = {{20, -42}, {10, -42}, {10, -42}, {10, -42}, {10, -42}}, color = {0, 0, 127}));
    connect(feedback1.y, PID.u) annotation(
      Line(points = {{-30, -46}, {-26, -46}, {-26, -16}, {-80, -16}, {-80, 0}, {-64, 0}, {-64, 0}}, color = {0, 0, 127}));
    connect(InputSignal.y, feedback1.u1) annotation(
      Line(points = {{-74, -46}, {-48, -46}, {-48, -46}, {-48, -46}}, color = {0, 0, 127}));
    connect(speedSensor1.flange, Hull.flange_prop) annotation(
      Line(points = {{42, -42}, {56, -42}, {56, 0}, {56, 0}}, color = {0, 127, 0}));
    connect(Hull.flange_water, fixed1.flange) annotation(
      Line(points = {{76, 0}, {86, 0}, {86, -54}, {86, -54}}, color = {0, 127, 0}));
    connect(Propeller.flangeT, Hull.flange_prop) annotation(
      Line(points = {{36, 0}, {56, 0}, {56, 0}, {56, 0}}, color = {0, 127, 0}));
    annotation(
      uses(Modelica(version = "3.2.2")));
  end EngineControlTest;

  model testHENmonophase
    package secondaryMedium = Modelica.Media.Water.StandardWater;
    Modelica.Blocks.Sources.Constant Massflow(k = 2) annotation(
      Placement(visible = true, transformation(origin = {68, 84}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant cp(k = 1018) annotation(
      Placement(visible = true, transformation(origin = {68, 56}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant Density(k = 1) annotation(
      Placement(visible = true, transformation(origin = {68, 24}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Step Temperature(height = 300, offset = 400, startTime = 7000) annotation(
      Placement(visible = true, transformation(origin = {68, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = secondaryMedium, T = 300, m_flow = 2, nPorts = 1) annotation(
      Placement(visible = true, transformation(origin = {-60, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT ambient(nPorts = 1, p = 150000, redeclare package Medium = secondaryMedium) annotation(
      Placement(visible = true, transformation(origin = {46, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    ShipEnergySystemsModelicaLibrary.Other.inputForFlangeCdot inputForFlangeCdot1 annotation(
      Placement(visible = true, transformation(origin = {24, 54}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    ShipEnergySystemsModelicaLibrary.HeatExchanger.HEX_1f HEX(A_cs_pf = 0.5, A_cs_sf = 1, A_ht = 20, L = 0.5, M_wall = 1, Tsf_SU_start = 873.15, U_pf = 600, U_sf = 50, c_wall = 450, dTsf_start(displayUnit = "K") = 20, eta_fin = 1, k_wall = 60, rho_wall(displayUnit = "kg/m3") = 5000) annotation(
      Placement(visible = true, transformation(origin = {-14, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(HEX.OutFlow, ambient.ports[1]) annotation(
      Line(points = {{-4, -2}, {20, -2}, {20, -38}, {36, -38}, {36, -38}, {36, -38}}, color = {0, 127, 255}));
    connect(boundary.ports[1], HEX.InFlow) annotation(
      Line(points = {{-50, -4}, {-24, -4}, {-24, -2}, {-24, -2}}, color = {0, 127, 255}, thickness = 0.5));
    connect(inputForFlangeCdot1.OutFlow_sf, HEX.InFlow_sf) annotation(
      Line(points = {{14, 60}, {6, 60}, {6, 10}, {-4, 10}, {-4, 10}}, color = {255, 0, 0}));
    connect(Temperature.y, inputForFlangeCdot1.u4[1]) annotation(
      Line(points = {{58, -6}, {44, -6}, {44, 46}, {36, 46}, {36, 46}}, color = {0, 0, 127}));
    connect(Density.y, inputForFlangeCdot1.u3[1]) annotation(
      Line(points = {{58, 24}, {50, 24}, {50, 52}, {36, 52}, {36, 52}}, color = {0, 0, 127}));
    connect(cp.y, inputForFlangeCdot1.u2[1]) annotation(
      Line(points = {{58, 56}, {36, 56}, {36, 58}, {36, 58}}, color = {0, 0, 127}));
    connect(Massflow.y, inputForFlangeCdot1.u1[1]) annotation(
      Line(points = {{58, 84}, {48, 84}, {48, 62}, {36, 62}, {36, 64}}, color = {0, 0, 127}));
    annotation(
      uses(Modelica(version = "3.2.2")));
  end testHENmonophase;


  model testHENtwophase
    package secondaryMedium = Modelica.Media.Water.StandardWater;
    Modelica.Blocks.Sources.Constant Massflow(k = 10) annotation(
      Placement(visible = true, transformation(origin = {68, 84}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant cp(k = 1018) annotation(
      Placement(visible = true, transformation(origin = {68, 56}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant Density(k = 1) annotation(
      Placement(visible = true, transformation(origin = {68, 24}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Step Temperature(height = 300, offset = 400, startTime = 7000) annotation(
      Placement(visible = true, transformation(origin = {68, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = secondaryMedium, T = 300, m_flow = 1, nPorts = 1) annotation(
      Placement(visible = true, transformation(origin = {-60, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT ambient(nPorts = 1, p = 150000, redeclare package Medium = secondaryMedium) annotation(
      Placement(visible = true, transformation(origin = {46, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    ShipEnergySystemsModelicaLibrary.Other.inputForFlangeCdot inputForFlangeCdot1 annotation(
      Placement(visible = true, transformation(origin = {24, 54}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    ShipEnergySystemsModelicaLibrary.HeatExchanger.HEX_2f HEX(A = 50, L = 10, M_wall = 5, U_pf_ev = 2600, U_pf_sl = 450, U_sf = 50, c_wall = 450, k_wall = 45, rho_wall(displayUnit = "kg/m3") = 5000, L_ev(start = 0.01), mdot_pf_sl2ev(start = 0.2), T_w_sf_in(start = 350), T_w_pf_in(start = 350), T_w_sl2ev(start = 350)) annotation(
      Placement(visible = true, transformation(origin = {-14, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(HEX.OutFlow, ambient.ports[1]) annotation(
      Line(points = {{-4, -2}, {18, -2}, {18, -38}, {36, -38}, {36, -38}}, color = {0, 127, 255}));
    connect(HEX.InFlow_sf, inputForFlangeCdot1.OutFlow_sf) annotation(
      Line(points = {{-4, 10}, {8, 10}, {8, 60}, {12, 60}, {12, 60}, {14, 60}}, color = {255, 0, 0}));
    connect(boundary.ports[1], HEX.InFlow) annotation(
      Line(points = {{-50, -4}, {-32, -4}, {-32, -2}, {-24, -2}, {-24, -2}}, color = {0, 127, 255}, thickness = 0.5));
    connect(Temperature.y, inputForFlangeCdot1.u4[1]) annotation(
      Line(points = {{58, -6}, {44, -6}, {44, 46}, {36, 46}, {36, 46}}, color = {0, 0, 127}));
    connect(Density.y, inputForFlangeCdot1.u3[1]) annotation(
      Line(points = {{58, 24}, {50, 24}, {50, 52}, {36, 52}, {36, 52}}, color = {0, 0, 127}));
    connect(cp.y, inputForFlangeCdot1.u2[1]) annotation(
      Line(points = {{58, 56}, {36, 56}, {36, 58}, {36, 58}}, color = {0, 0, 127}));
    connect(Massflow.y, inputForFlangeCdot1.u1[1]) annotation(
      Line(points = {{58, 84}, {48, 84}, {48, 62}, {36, 62}, {36, 64}}, color = {0, 0, 127}));
    annotation(
      uses(Modelica(version = "3.2.2")));
  end testHENtwophase;
  
  model testHENthreePhase
  package secondaryMedium = Modelica.Media.Water.StandardWater;
  Modelica.Blocks.Sources.Constant Massflow(k = 10) annotation(
    Placement(visible = true, transformation(origin = {68, 84}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant cp(k = 1018) annotation(
    Placement(visible = true, transformation(origin = {68, 56}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant Density(k = 1) annotation(
    Placement(visible = true, transformation(origin = {68, 24}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step Temperature(height = 300, offset = 500, startTime = 7000) annotation(
    Placement(visible = true, transformation(origin = {68, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = secondaryMedium, T = 300, m_flow = 1, nPorts = 1) annotation(
    Placement(visible = true, transformation(origin = {-60, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT ambient(nPorts = 1, p = 150000, redeclare package Medium = secondaryMedium) annotation(
    Placement(visible = true, transformation(origin = {46, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  ShipEnergySystemsModelicaLibrary.Other.inputForFlangeCdot inputForFlangeCdot1 annotation(
    Placement(visible = true, transformation(origin = {24, 54}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  ShipEnergySystemsModelicaLibrary.HeatExchanger.HEX_3f hEX_3f1(A = 10, ETA = 1, L = 2, L_SB(start = 1.98), L_TP(start = 0.01), M_tot = 5, Tsf_SU_start = 623.15, TwSB(start = 373.15), TwSH(start = 393.15), TwTP(start = 383.15), U_SB = 160, U_SH = 80, U_TP = 500, Usf = 60, Void = 0.5, c_wall = 500, dTsf_start(displayUnit = "K") = 20, dVoid_dh = 0, dVoid_dp = 0, h_EX(start = 100000), p(start = 100000), rho_wall(displayUnit = "kg/m3") = 6000)  annotation(
      Placement(visible = true, transformation(origin = {-10, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(hEX_3f1.InFlow, boundary.ports[1]) annotation(
      Line(points = {{-20, -2}, {-42, -2}, {-42, -4}, {-50, -4}, {-50, -4}}, color = {0, 127, 255}));
    connect(hEX_3f1.OutFlow, ambient.ports[1]) annotation(
      Line(points = {{0, -2}, {24, -2}, {24, -38}, {36, -38}, {36, -38}}, color = {0, 127, 255}));
    connect(hEX_3f1.InFlow_sf, inputForFlangeCdot1.OutFlow_sf) annotation(
      Line(points = {{0, 10}, {6, 10}, {6, 60}, {14, 60}, {14, 60}}, color = {255, 0, 0}));
    connect(Temperature.y, inputForFlangeCdot1.u4[1]) annotation(
      Line(points = {{58, -6}, {44, -6}, {44, 46}, {36, 46}, {36, 46}}, color = {0, 0, 127}));
    connect(Density.y, inputForFlangeCdot1.u3[1]) annotation(
      Line(points = {{58, 24}, {50, 24}, {50, 52}, {36, 52}, {36, 52}}, color = {0, 0, 127}));
    connect(cp.y, inputForFlangeCdot1.u2[1]) annotation(
      Line(points = {{58, 56}, {36, 56}, {36, 58}, {36, 58}}, color = {0, 0, 127}));
    connect(Massflow.y, inputForFlangeCdot1.u1[1]) annotation(
      Line(points = {{58, 84}, {48, 84}, {48, 62}, {36, 62}, {36, 64}}, color = {0, 0, 127}));
    annotation(
    uses(Modelica(version = "3.2.2")));
end testHENthreePhase;

model testHENthreePhaseSS

package secondaryMedium = Modelica.Media.Water.StandardWater;
Modelica.Blocks.Sources.Constant Massflow(k = 10) annotation(
  Placement(visible = true, transformation(origin = {68, 84}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
Modelica.Blocks.Sources.Constant cp(k = 1018) annotation(
  Placement(visible = true, transformation(origin = {68, 56}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
Modelica.Blocks.Sources.Constant Density(k = 1) annotation(
  Placement(visible = true, transformation(origin = {68, 24}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
Modelica.Blocks.Sources.Step Temperature(height = 300, offset = 500, startTime = 7000) annotation(
  Placement(visible = true, transformation(origin = {68, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = secondaryMedium, T = 300, m_flow = 1, nPorts = 1) annotation(
  Placement(visible = true, transformation(origin = {-60, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
Modelica.Fluid.Sources.Boundary_pT ambient(nPorts = 1, p = 150000, redeclare package Medium = secondaryMedium) annotation(
  Placement(visible = true, transformation(origin = {46, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
ShipEnergySystemsModelicaLibrary.Other.inputForFlangeCdot inputForFlangeCdot1 annotation(
  Placement(visible = true, transformation(origin = {24, 54}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
ShipEnergySystemsModelicaLibrary.HeatExchanger.HEX_3f_ss hEX_3f1(A = 10, ETA = 1, L = 2, L_SB(start = 1.98), L_TP(start = 0.01), M_tot = 5, Tsf_SU_start = 623.15, U_wf_SB = 160, U_wf_SH = 80, U_wf_TP = 500, U_sf = 60, h_EX(start = 100000), p(start = 100000))  annotation(
    Placement(visible = true, transformation(origin = {-10, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(hEX_3f1.InFlow, boundary.ports[1]) annotation(
    Line(points = {{-20, -2}, {-42, -2}, {-42, -4}, {-50, -4}, {-50, -4}}, color = {0, 127, 255}));
  connect(hEX_3f1.OutFlow, ambient.ports[1]) annotation(
    Line(points = {{0, -2}, {24, -2}, {24, -38}, {36, -38}, {36, -38}}, color = {0, 127, 255}));
  connect(hEX_3f1.InFlow_sf, inputForFlangeCdot1.OutFlow_sf) annotation(
    Line(points = {{0, 10}, {6, 10}, {6, 60}, {14, 60}, {14, 60}}, color = {255, 0, 0}));
  connect(Temperature.y, inputForFlangeCdot1.u4[1]) annotation(
    Line(points = {{58, -6}, {44, -6}, {44, 46}, {36, 46}, {36, 46}}, color = {0, 0, 127}));
  connect(Density.y, inputForFlangeCdot1.u3[1]) annotation(
    Line(points = {{58, 24}, {50, 24}, {50, 52}, {36, 52}, {36, 52}}, color = {0, 0, 127}));
  connect(cp.y, inputForFlangeCdot1.u2[1]) annotation(
    Line(points = {{58, 56}, {36, 56}, {36, 58}, {36, 58}}, color = {0, 0, 127}));
  connect(Massflow.y, inputForFlangeCdot1.u1[1]) annotation(
    Line(points = {{58, 84}, {48, 84}, {48, 62}, {36, 62}, {36, 64}}, color = {0, 0, 127}));
  annotation(
  uses(Modelica(version = "3.2.2")));

end testHENthreePhaseSS;




  
end Examples;
