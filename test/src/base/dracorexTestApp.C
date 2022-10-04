//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "dracorexTestApp.h"
#include "dracorexApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
dracorexTestApp::validParams()
{
  InputParameters params = dracorexApp::validParams();
  return params;
}

dracorexTestApp::dracorexTestApp(InputParameters parameters) : MooseApp(parameters)
{
  dracorexTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

dracorexTestApp::~dracorexTestApp() {}

void
dracorexTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  dracorexApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"dracorexTestApp"});
    Registry::registerActionsTo(af, {"dracorexTestApp"});
  }
}

void
dracorexTestApp::registerApps()
{
  registerApp(dracorexApp);
  registerApp(dracorexTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
dracorexTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  dracorexTestApp::registerAll(f, af, s);
}
extern "C" void
dracorexTestApp__registerApps()
{
  dracorexTestApp::registerApps();
}
