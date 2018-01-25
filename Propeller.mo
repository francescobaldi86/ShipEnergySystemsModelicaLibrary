within ShipEnergySystems;

package Propeller "Package for propeller models"
  // Importing some useful packages //
  import SI = Modelica.SIunits;

  partial model BasicPropeller "Basic propeller model"
    extends Modelica.Mechanics.Rotational.Interfaces.PartialElementaryRotationalToTranslational;
    annotation(Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-100,15},{20,-15}},
          lineColor={0,24,48},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={0,127,255},
          radius=10),
        Text(
          extent={{60,0},{100,40}},
          lineColor={0,0,0},
          textString="T"),
        Text(
          extent={{-100,40},{-60,80}},
          lineColor={0,0,0},
          textString="R"),
        Ellipse(
          extent={{20,-90},{60,-10}},
          lineColor={0,24,48},
          fillPattern=FillPattern.Sphere,
          fillColor={0,127,255}),
        Ellipse(
          extent={{20,10},{60,90}},
          lineColor={0,24,48},
          fillPattern=FillPattern.Sphere,
          fillColor={0,127,255}),
        Ellipse(
          extent={{30,-10},{50,10}},
          lineColor={0,24,48},
          fillPattern=FillPattern.Sphere,
          fillColor={0,127,255})}),
      Documentation(info = "<html>
<p>
A propeller is substantially nothing more than an element that converts rotational motion into translational
</p>

</html>"));
    // Model-specific parameters
    // Note: the model is partial because the efficiency is missing. This will be defined by the specific models
    // Model-specific variables
    SI.Angle phi(start=0) "Absolute rotation angle of component" ;
    // annotation (Dialog(group="Initialization", showStartAttribute=true)); //
    SI.AngularVelocity omega(start = 0) "Absolute angular velocity of component (= der(phi))" ;
    // annotation (Dialog(group="Initialization", showStartAttribute=true)); //
    SI.Torque Q "Rotational torque (= tau on flangeR)"; 
    // annotation (Dialog(group="Initialization", showStartAttribute=true)); //
    SI.Position s(start=0) "Position of the propeller" ;
    // annotation (Dialog(group="Initialization", showStartAttribute=true)); //
    SI.Velocity v(start=0) "Speed of the water through the propeller " ;
    // annotation (Dialog(group="Initialization", showStartAttribute=true)); //
    SI.Force T "Forward trhust delivered by the propeller" ;
    // annotation (Dialog(group="Initialization", showStartAttribute=true)); //
    // Model-specific equations
    equation
      // Propeller shaft side //
      phi = flangeR.phi;
      omega = der(phi); // Propeller shaft speed //
      Q = flangeR.tau; // Propeller torque (shaft side)
      // Water side //
      s = flangeT.s;
      v = der(s); // Propeller speed through water //
      T = flangeT.f ; // Propeller thrust //
      
  end BasicPropeller;






  
  model FixedEfficiency "Propeller with fixed, user-provided efficiency"
    extends BasicPropeller;
      // Model-specific parameters
    parameter Real K_T(start = 0.25);
    parameter Real K_Q(start = 0.035);
    parameter SI.Length diameter(start = 5);
      // Model-specific variables
    Real n(unit="round/s");
      // Model-specific equations
    equation
      n = omega / 2 / Modelica.Constants.pi ;
      T = K_T * 1000 * n^2 * diameter^4 ;
      Q = K_Q * 1000 * n^3 * diameter^4 ;
  
  end FixedEfficiency;













  
  
  
end Propeller;
