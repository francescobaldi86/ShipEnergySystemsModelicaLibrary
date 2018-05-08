within ShipEnergySystemsModelicaLibrary;

package ElectricMachinery

model DCmotor
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
    Placement(visible = true, transformation(origin = {-48, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Inductor inductor1(L = 5)  annotation(
    Placement(visible = true, transformation(origin = {20, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Resistor resistor1(R = 10)  annotation(
    Placement(visible = true, transformation(origin = {-22, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.EMF emf(k = -1)  annotation(
    Placement(visible = true, transformation(origin = {48, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage1 annotation(
    Placement(visible = true, transformation(origin = {-48, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation(
    Placement(visible = true, transformation(origin = {80, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u annotation(
    Placement(visible = true, transformation(origin = {-100, 6}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 6}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
equation
  connect(signalVoltage1.v, u) annotation(
    Line(points = {{-56, 6}, {-82, 6}, {-82, 6}, {-100, 6}}, color = {0, 0, 127}));
  connect(emf.flange, flange_a) annotation(
    Line(points = {{58, 2}, {68, 2}, {68, 2}, {78, 2}, {78, 2}, {80, 2}}));
  connect(ground1.p, signalVoltage1.p) annotation(
    Line(points = {{-48, -38}, {-48, -4}}, color = {0, 0, 255}));
  connect(signalVoltage1.n, resistor1.p) annotation(
    Line(points = {{-48, 16}, {-48, 40}, {-32, 40}}, color = {0, 0, 255}));
  connect(emf.n, ground1.p) annotation(
    Line(points = {{48, -8}, {-48, -8}, {-48, -38}}, color = {0, 0, 255}));
  connect(inductor1.n, emf.p) annotation(
    Line(points = {{30, 40}, {48, 40}, {48, 12}, {48, 12}, {48, 12}, {48, 12}, {48, 12}}, color = {0, 0, 255}));
  connect(resistor1.n, inductor1.p) annotation(
    Line(points = {{-12, 40}, {8, 40}, {8, 40}, {10, 40}, {10, 40}}, color = {0, 0, 255}));
  annotation(
    uses(Modelica(version = "3.2.2")), Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          lineColor={82,0,2},
          fillColor={252,37,57},
          fillPattern=FillPattern.HorizontalCylinder,
          extent={{-100.0,-50.0},{30.0,50.0}},
          radius=10.0),
        Polygon(
          fillColor={64,64,64},
          fillPattern=FillPattern.Solid,
          points={{-100.0,-90.0},{-90.0,-90.0},{-60.0,-20.0},{-10.0,-20.0},{20.0,-90.0},{30.0,-90.0},{30.0,-100.0},{-100.0,-100.0},{-100.0,-90.0}}),
        Rectangle(
          lineColor={64,64,64},
          fillColor={255,255,255},
          fillPattern=FillPattern.HorizontalCylinder,
          extent={{30.0,-10.0},{90.0,10.0}})}));
end DCmotor;


end ElectricMachinery;
