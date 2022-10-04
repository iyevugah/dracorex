#include "dracorexApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
dracorexApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

dracorexApp::dracorexApp(InputParameters parameters) : MooseApp(parameters)
{
  dracorexApp::registerAll(_factory, _action_factory, _syntax);
}

dracorexApp::~dracorexApp() {}

void
dracorexApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"dracorexApp"});
  Registry::registerActionsTo(af, {"dracorexApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
dracorexApp::registerApps()
{
  registerApp(dracorexApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
dracorexApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  dracorexApp::registerAll(f, af, s);
}
extern "C" void
dracorexApp__registerApps()
{
  dracorexApp::registerApps();
}
