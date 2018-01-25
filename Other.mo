within ShipEnergySystems;

package Other

connector FlangeCdot_A "Flange for a Cdot-type heat source. It can be used for all fluids with a known and fixed cp and density"
  input Modelica.SIunits.MassFlowRate Mdot "Mass flow rate";
  input Modelica.SIunits.SpecificHeatCapacity cp "Specific Heat capacity";
  input Modelica.SIunits.Density rho "Density of entering fluid";
  input Modelica.SIunits.Temperature T "Temperature of entering fluid";
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={Ellipse(extent={{-100,100},{100,-100}},
            lineColor={255,0,0})}));
end FlangeCdot_A;

connector FlangeCdot_B "B-type flange connector for Cdot-type heat source"
  output Modelica.SIunits.MassFlowRate Mdot "Mass flow rate";
  output Modelica.SIunits.SpecificHeatCapacity cp "Specific Heat capacity";
  output Modelica.SIunits.Density rho "Density of the fluid";
  output Modelica.SIunits.Temperature T "Temperature of the fluid";
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={Ellipse(extent={{-100,100},{100,-100}},
            lineColor={255,0,0}), Ellipse(
          extent={{-42,44},{44,-40}},
          lineColor={255,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end FlangeCdot_B;

model inputForFlangeCdot "This block takes four separate inputs and merges them in the output information required by the flange of Cdot type"
  ShipEnergySystems.Other.FlangeCdot_A OutFlow_sf annotation (Placement(transformation(extent={{110,50},{90,70}})));
  Modelica.Blocks.Interfaces.RealInput u1[1]
    "Connector of Real input signals 1" annotation (Placement(transformation(
          extent={{-140,70},{-100,110}})));
  Modelica.Blocks.Interfaces.RealInput u2[1]
    "Connector of Real input signals 2" annotation (Placement(transformation(
          extent={{-140,10},{-100,50}})));
  Modelica.Blocks.Interfaces.RealInput u3[1]
    "Connector of Real input signals 3" annotation (Placement(transformation(
          extent={{-140,-50},{-100,-10}})));
  Modelica.Blocks.Interfaces.RealInput u4[1]
    "Connector of Real input signals 4" annotation (Placement(transformation(
          extent={{-140,-110},{-100,-70}})));
  
equation
  OutFlow_sf.Mdot = u1[1];
  OutFlow_sf.cp = u2[1];
  OutFlow_sf.rho = u3[1];
  OutFlow_sf.T = u4[1];
end inputForFlangeCdot;

end Other;
