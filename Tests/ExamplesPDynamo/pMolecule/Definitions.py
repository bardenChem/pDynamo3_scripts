"""Some definitions."""

import os, os.path

from pCore import TestScript_InputDataPath  , \
                  TestScript_OutputDataPath

# . Local name.
_name = "pMolecule"

# . Paths.
dataPath       = TestScript_InputDataPath  ( _name )
outPath        = TestScript_OutputDataPath ( _name )
structuresPath = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "structures" )
yamlPath       = os.path.join ( dataPath, "yaml" )

# . Other options.
_FullVerificationSummary = False
