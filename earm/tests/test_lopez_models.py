"""
Checks the MOMP-only versions of the Lopez models, that is

- :py:mod:`earm.mito.lopez_embedded`
- :py:mod:`earm.mito.lopez_direct`
- :py:mod:`earm.mito.lopez_indirect`

against a previously validated and serialized state. Each test in the module
uses the function pysb.testing.check_model_against_component_list to perform
the comparison.
"""

from earm.mito import lopez_embedded
from earm.mito import lopez_direct
from earm.mito import lopez_indirect

from pysb.testing import *
import pickle
import os

def test_lopez_embedded():
    """Test the earm.mito.lopez_embedded model against a serialized state."""
    full_pickle_file_path = os.path.join(os.path.dirname(__file__),
                    'earm.mito.lopez_embedded_validated.pickle')
    f = open(full_pickle_file_path, 'r')
    component_list = pickle.load(f)
    check_model_against_component_list(lopez_embedded.model, component_list)

def test_lopez_direct():
    """Test the earm.mito.lopez_direct model against a serialized state."""
    full_pickle_file_path = os.path.join(os.path.dirname(__file__),
                    'earm.mito.lopez_direct_validated.pickle')
    f = open(full_pickle_file_path, 'r')
    component_list = pickle.load(f)
    check_model_against_component_list(lopez_direct.model, component_list)

def test_lopez_indirect():
    """Test the earm.mito.lopez_indirect model against a serialized state."""
    full_pickle_file_path = os.path.join(os.path.dirname(__file__),
                    'earm.mito.lopez_indirect_validated.pickle')
    f = open(full_pickle_file_path, 'r')
    component_list = pickle.load(f)
    check_model_against_component_list(lopez_indirect.model, component_list)

def pickle_lopez_models():
    """The pickling procedure that was used to serialize the components
    of the Lopez models as a record of a validated state.
    """
    serialize_component_list(lopez_embedded.model,
                             'lopez_embedded_validated.pickle')
    serialize_component_list(lopez_direct.model,
                             'lopez_direct_validated.pickle')
    serialize_component_list(lopez_indirect.model,
                             'lopez_indirect_validated.pickle')
