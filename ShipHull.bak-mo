within ShipEnergySystems;

package ShipHull "Package for ship hull models"
  // Importing some useful packages //
  import SI = Modelica.SIunits;

  partial model BasicHull "Basic hull model"
    Modelica.Mechanics.Translational.Interfaces.Flange_a flange_prop "Hull connection to the propeller" annotation (
        Placement(transformation(extent={{-110,-10},{-90,10}})));
    Modelica.Mechanics.Translational.Interfaces.Flange_b flange_water "Hull connection to the water" annotation (
        Placement(transformation(extent={{90,-10},{110,10}})));
    // Model components
    Modelica.Mechanics.Translational.Components.Mass Displacement annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShipResistance Resistance annotation(
      Placement(visible = true, transformation(origin = {12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // Model-specific parameters
    // Model-specific variables
    equation
    connect(Resistance.flange_b, flange_water) annotation(
      Line(points = {{22, 0}, {98, 0}, {98, 0}, {100, 0}}, color = {0, 127, 0}));
    connect(Displacement.flange_b, Resistance.flange_a) annotation(
      Line(points = {{-32, 0}, {2, 0}, {2, 0}, {2, 0}}, color = {0, 127, 0}));
    connect(flange_prop, Displacement.flange_a) annotation(
      Line(points = {{-100, 0}, {-52, 0}, {-52, 0}, {-52, 0}}, color = {0, 127, 0}));
    annotation(
      Documentation(info = "<html>
<p>
We assume that the basic hull model is simply an inertia, with no friction
</p>

</html>"));
  end BasicHull;





























  model FixedRt "Basic hull model with externally provided ship resistance"
    extends BasicHull;
    // Model-specific parameters
    parameter SI.Force R_T "Ship resistance";
    // Model-specific equations
  equation
    f = R_T;
    annotation(
      Documentation(info = "<html>
<p>
In this case, we assume RT (TowingResistance) to fixed and provided by the user
</p>

</html>"));
  end FixedRt;

  model FixedC "Basic hull model where the C1 coefficient for the towing resistance is provided"
    extends BasicHull;
    // Model-specific parameters
    parameter Real C1(quantity = "Resistance Coefficient", unit = "N/m^2");
    // Model-specific equations
  equation
    Resistance.f = C1 * Resistance.v ^ 2;//
    annotation(
      Documentation(info = "<html>
<p>
In this case, we assume RT (TowingResistance) to be equal to R_T = C1 * v^2. C1 is provided externally by the user
</p>

</html>"));
  end FixedC;






  model ShipResistance
    extends Modelica.Mechanics.Translational.Interfaces.PartialCompliant;
    // Model-specific components
    // Model-specific parameters
    // Model-specific variables
    SI.Velocity v;
    // Model-specific equations
  equation
    v = der(s_rel);
  end ShipResistance;















end ShipHull;
